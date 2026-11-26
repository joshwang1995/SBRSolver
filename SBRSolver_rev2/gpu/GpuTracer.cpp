#include "gpu/GpuTracer.h"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "CpuTracer.h"
#include "common/Constants.h"
#include "common/Profiling.h"
#include "gpu/GpuRecords.h"
#include "trace_double.spv.h"
#include "trace_single.spv.h"

namespace sbr
{
namespace gpu
{
namespace
{
	// The trace shader walks the BVH with a stack of this many entries.
	const int kMaxBvhDepth = 63;

	struct TracePush
	{
		DispatchBase base;
		uint32_t pendingFirst;
		uint32_t pendingCount;
		int32_t maxReflections;
		int32_t maxTransmissions;
	};

	// Segments a launch ray can produce with the children rule of RayLaunch:
	// one per sequence of r reflections and t transmissions, r <= R, t <= T.
	uint64_t TreeBound(int maxReflections, int maxTransmissions)
	{
		// Pascal's triangle by rows, so no factorial overflows.
		std::vector<uint64_t> row(size_t(maxTransmissions) + 1, 1);
		uint64_t total = row.size();
		for (int r = 1; r <= maxReflections; r++)
		{
			for (int t = 1; t <= maxTransmissions; t++) row[t] += row[t - 1];
			for (uint64_t value : row) total += value;
		}
		return total;
	}

	int Depth(const BVH<Triangle>& bvh, int node)
	{
		if (bvh.nodes[node].children[0] < 0) return 1;
		return 1 + std::max(Depth(bvh, bvh.nodes[node].children[0]), Depth(bvh, bvh.nodes[node].children[1]));
	}

	template <typename Real>
	struct FlatBvh
	{
		std::vector<GpuNode<Real>> nodes;
		std::vector<GpuTriangle<Real>> triangles;
		std::vector<uint32_t> leafTriangles;
	};

	// The CPU's BVH as index arrays, in the same node order, with the leaf
	// triangles listed in the leaves' stored order.
	template <typename Real>
	FlatBvh<Real> Flatten(const BVH<Triangle>& bvh, const std::vector<Triangle*>& mesh)
	{
		const int depth = Depth(bvh, 0);
		if (depth > kMaxBvhDepth)
		{
			throw std::runtime_error("the BVH is " + std::to_string(depth) + " levels deep; the GPU tracer handles "
				+ std::to_string(kMaxBvhDepth) + " (run with --device cpu)");
		}
		FlatBvh<Real> flat;
		flat.nodes.reserve(bvh.nodes.size());
		for (const BVHNode<Triangle>& node : bvh.nodes)
		{
			GpuNode<Real> record{};
			for (int axis = 0; axis < 3; axis++)
			{
				record.bmin[axis] = Real(node.bbox(0, axis));
				record.bmax[axis] = Real(node.bbox(1, axis));
			}
			if (node.children[0] >= 0)
			{
				record.left = node.children[0];
				record.right = node.children[1];
			}
			else
			{
				record.first = int32_t(flat.leafTriangles.size());
				record.count = int32_t(node.objects.size());
				for (const Triangle* triangle : node.objects) flat.leafTriangles.push_back(uint32_t(triangle->triangleId));
			}
			flat.nodes.push_back(record);
		}
		flat.triangles.reserve(mesh.size());
		for (const Triangle* triangle : mesh)
		{
			// The edges are the CPU's doubles, rounded once for the single build.
			const Vec3 e1 = triangle->v2 - triangle->v1;
			const Vec3 e2 = triangle->v3 - triangle->v1;
			GpuTriangle<Real> record{};
			for (int axis = 0; axis < 3; axis++)
			{
				record.v1[axis] = Real(triangle->v1[axis]);
				record.e1[axis] = Real(e1[axis]);
				record.e2[axis] = Real(e2[axis]);
				record.n[axis] = Real(triangle->norm[axis]);
			}
			record.id = triangle->triangleId;
			record.coplanar = triangle->coplanarId;
			record.material = triangle->materialId;
			flat.triangles.push_back(record);
		}
		return flat;
	}

	template <typename T>
	void UploadVector(GpuContext& gpu, Buffer& buffer, const std::vector<T>& data)
	{
		gpu.Upload(buffer, data.data(), uint64_t(data.size()) * sizeof(T));
	}

	// Turns the records of the device, in append order, into the CPU's
	// segments: launch by launch, breadth first, the transmit child before
	// the reflect child, with the Ray fields as RayLaunch sets them.
	template <typename Real>
	void Renumber(const std::vector<GpuSegment<Real>>& records, int launches, TraceResult& result)
	{
		std::vector<std::vector<int>> children(records.size());
		for (size_t s = size_t(launches); s < records.size(); s++)
		{
			children[size_t(records[s].parent)].push_back(int(s));
		}
		for (std::vector<int>& list : children)
		{
			std::sort(list.begin(), list.end(), [&records](int a, int b) { return records[size_t(a)].kind > records[size_t(b)].kind; });
		}

		std::vector<int> newIndex(records.size(), -1);
		result.segments.clear();
		result.segments.reserve(records.size());
		result.launchFirst.assign(size_t(launches) + 1, 0);
		std::deque<int> queue;
		for (int launch = 0; launch < launches; launch++)
		{
			result.launchFirst[size_t(launch)] = int(result.segments.size());
			queue.push_back(launch);
			while (!queue.empty())
			{
				const int old = queue.front();
				queue.pop_front();
				const GpuSegment<Real>& record = records[size_t(old)];
				newIndex[size_t(old)] = int(result.segments.size());

				Segment segment;
				segment.launch = record.launch;
				segment.parent = (record.parent >= 0) ? newIndex[size_t(record.parent)] : -1;
				Ray& ray = segment.ray;
				ray.sourcePoint = Vec3(double(record.ox), double(record.oy), double(record.oz));
				const Vec3 direction(double(record.dx), double(record.dy), double(record.dz));
				ray.targetPoint = ray.sourcePoint + direction * double(record.t);
				if (record.parent >= 0)
				{
					const int32_t parentMaterial = records[size_t(record.parent)].material;
					ray.reflectionMaterialId = (record.kind == kKindReflect) ? parentMaterial : -1;
					ray.penetrationMaterialId = (record.kind == kKindTransmit) ? parentMaterial : -1;
				}
				ray.hitSurfaceID = record.triangle;
				ray.hitCoplanarId = record.coplanar;
				ray.pathLength = (record.triangle >= 0) ? double(record.unfoldedStart) + (ray.targetPoint - ray.sourcePoint).norm() : INF;
				result.segments.push_back(segment);
				for (int child : children[size_t(old)]) queue.push_back(child);
			}
		}
		result.launchFirst[size_t(launches)] = int(result.segments.size());
	}

	template <typename Real>
	void TraceOnGpu(const TracerInput& input, GpuContext& gpu, const uint32_t* spirv, size_t words, TraceResult& result)
	{
		SBR_ZONE_NAMED("GpuTrace");
		const FlatBvh<Real> flat = Flatten<Real>(*input.bvh, *input.triangles);
		const int launches = int(input.launches->size());
		const uint64_t capacity = uint64_t(launches) * TreeBound(input.maxReflections, input.maxTransmissions);

		Buffer nodes = gpu.CreateBuffer(uint64_t(flat.nodes.size()) * sizeof(GpuNode<Real>));
		Buffer triangles = gpu.CreateBuffer(uint64_t(flat.triangles.size()) * sizeof(GpuTriangle<Real>));
		Buffer leaves = gpu.CreateBuffer(uint64_t(std::max<size_t>(flat.leafTriangles.size(), 1)) * sizeof(uint32_t));
		Buffer segments = gpu.CreateBuffer(capacity * sizeof(GpuSegment<Real>));
		Buffer counter = gpu.CreateBuffer(sizeof(uint32_t));
		UploadVector(gpu, nodes, flat.nodes);
		UploadVector(gpu, triangles, flat.triangles);
		if (!flat.leafTriangles.empty()) UploadVector(gpu, leaves, flat.leafTriangles);

		// The launch rays are the roots; the shader appends every child
		// behind them, so one level is one contiguous range.
		std::vector<GpuSegment<Real>> roots(static_cast<size_t>(launches));
		for (int i = 0; i < launches; i++)
		{
			const Vec3& origin = (*input.transmitters)[size_t((*input.launchTx)[size_t(i)])].position;
			const Vec3& direction = (*input.launches)[size_t(i)];
			GpuSegment<Real>& root = roots[size_t(i)];
			root.ox = Real(origin.x());
			root.oy = Real(origin.y());
			root.oz = Real(origin.z());
			root.unfoldedStart = Real(0.0);
			root.dx = Real(direction.x());
			root.dy = Real(direction.y());
			root.dz = Real(direction.z());
			root.t = Real(0.0);
			root.parent = -1;
			root.launch = i;
			root.triangle = -1;
			root.coplanar = -1;
			root.material = -1;
			root.kind = kKindRoot;
			root.reflections = 0;
			root.transmissions = 0;
		}
		UploadVector(gpu, segments, roots);
		uint32_t count = uint32_t(launches);
		gpu.Upload(counter, &count, sizeof(count));

		Pipeline trace = gpu.CreateCompute(spirv, words, 5, sizeof(TracePush));
		uint32_t levelStart = 0;
		uint32_t levelEnd = count;
		for (int level = 0; levelEnd > levelStart; level++)
		{
			TracePush push{};
			push.pendingFirst = levelStart;
			push.pendingCount = levelEnd - levelStart;
			push.maxReflections = input.maxReflections;
			push.maxTransmissions = input.maxTransmissions;
			gpu.Dispatch(trace, { &nodes, &triangles, &leaves, &segments, &counter }, &push, sizeof(push), push.pendingCount);
			std::cout << "\tGPU trace level " << level << " -> " << push.pendingCount << " segments" << std::endl;
			gpu.Download(counter, &count, sizeof(count));
			levelStart = levelEnd;
			levelEnd = count;
		}

		std::vector<GpuSegment<Real>> records(count);
		gpu.Download(segments, records.data(), uint64_t(count) * sizeof(GpuSegment<Real>));
		Renumber(records, launches, result);
	}
}

void TraceAndCapture(const TracerInput& input, GpuContext& gpu, GpuPrecision precision, TraceResult& result)
{
	if (precision == GpuPrecision::Double)
	{
		TraceOnGpu<double>(input, gpu, kSpirv_trace_double, kSpirvWords_trace_double, result);
	}
	else
	{
		TraceOnGpu<float>(input, gpu, kSpirv_trace_single, kSpirvWords_trace_single, result);
	}
	std::cout << "\tnote: GPU trace, CPU capture" << std::endl;
	CaptureCandidates(input, result);
}

} // namespace gpu
} // namespace sbr
