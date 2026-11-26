#include "gpu/GpuTracer.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "common/Constants.h"
#include "common/Profiling.h"
#include "gpu/GpuRecords.h"
#include "capture_double.spv.h"
#include "capture_single.spv.h"
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
	// A batch of launches may use this much device memory for its segments,
	// unless the device limits allow less; SBR_GPU_BATCH_BYTES lowers it.
	const uint64_t kBatchBytes = uint64_t(1) << 30;
	// Candidate pairs the buffer holds at first; SBR_GPU_CANDIDATE_CAPACITY
	// sets it. A batch that captures more reruns with the count it needs.
	const uint64_t kCandidateCapacity = 4194304;

	struct TracePush
	{
		DispatchBase base;
		uint32_t pendingFirst;
		uint32_t pendingCount;
		int32_t maxReflections;
		int32_t maxTransmissions;
	};

	struct CapturePush
	{
		DispatchBase base;
		uint32_t segmentCount;
		uint32_t receiverFirst;
		uint32_t receiverCount;
		uint32_t capacity;
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
			for (int t = 1; t <= maxTransmissions; t++) row[size_t(t)] += row[size_t(t) - 1];
			for (uint64_t value : row) total += value;
		}
		return total;
	}

	int Depth(const BVH<Triangle>& bvh, int node)
	{
		if (bvh.nodes[size_t(node)].children[0] < 0) return 1;
		return 1 + std::max(Depth(bvh, bvh.nodes[size_t(node)].children[0]), Depth(bvh, bvh.nodes[size_t(node)].children[1]));
	}

	uint64_t EnvironmentValue(const char* name, uint64_t fallback)
	{
		const char* text = std::getenv(name);
		if (text == nullptr || *text == '\0') return fallback;
		char* end = nullptr;
		const unsigned long long value = std::strtoull(text, &end, 10);
		if (end == text || *end != '\0' || value == 0)
		{
			throw std::runtime_error(std::string(name) + " must be a positive integer, got '" + text + "'");
		}
		return value;
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
		if (!data.empty()) gpu.Upload(buffer, data.data(), uint64_t(data.size()) * sizeof(T));
	}

	template <typename T>
	Buffer CreateVectorBuffer(GpuContext& gpu, const std::vector<T>& data)
	{
		Buffer buffer = gpu.CreateBuffer(uint64_t(std::max<size_t>(data.size(), 1)) * sizeof(T));
		UploadVector(gpu, buffer, data);
		return buffer;
	}

	// One run of the GPU tracer in one precision: the scene on the device,
	// the two pipelines, and the batches of launches.
	template <typename Real>
	class Run
	{
	public:
		Run(const TracerInput& input, GpuContext& gpu, const uint32_t* traceSpirv, size_t traceWords, const uint32_t* captureSpirv, size_t captureWords)
			: _input(input), _gpu(gpu)
		{
			const FlatBvh<Real> flat = Flatten<Real>(*input.bvh, *input.triangles);
			_nodes = CreateVectorBuffer(gpu, flat.nodes);
			_triangles = CreateVectorBuffer(gpu, flat.triangles);
			_leaves = CreateVectorBuffer(gpu, flat.leafTriangles);

			std::vector<GpuReceiver<Real>> receivers;
			for (const Vec3& r : *input.receivers) receivers.push_back({ Real(r.x()), Real(r.y()), Real(r.z()) });
			_receivers = CreateVectorBuffer(gpu, receivers);
			std::vector<Real> slopes;
			for (int tx : *input.launchTx) slopes.push_back(Real(input.tanBeta[size_t(tx)]));
			_slopes = CreateVectorBuffer(gpu, slopes);
			_counter = gpu.CreateBuffer(sizeof(uint32_t));

			// A batch holds as many launches as its segment buffer allows.
			const DeviceInfo& info = gpu.Info();
			const uint64_t limit = std::min({ kBatchBytes, info.heapBytes / 4, info.storageRange, info.allocationLimit });
			const uint64_t batchBytes = std::min(limit, EnvironmentValue("SBR_GPU_BATCH_BYTES", limit));
			_bound = TreeBound(input.maxReflections, input.maxTransmissions);
			const uint64_t recordBytes = _bound * sizeof(GpuSegment<Real>);
			if (recordBytes > batchBytes)
			{
				std::ostringstream message;
				message << "one launch ray can produce " << _bound << " segments (" << recordBytes << " bytes), more than the "
					<< batchBytes << " bytes a batch may use; lower max_reflections or max_transmissions, or raise SBR_GPU_BATCH_BYTES";
				throw std::runtime_error(message.str());
			}
			const int launches = int(input.launches->size());
			_launchesPerBatch = int(std::min<uint64_t>(uint64_t(launches), batchBytes / recordBytes));
			_segments = gpu.CreateBuffer(uint64_t(_launchesPerBatch) * recordBytes);
			_capacity = std::min<uint64_t>(EnvironmentValue("SBR_GPU_CANDIDATE_CAPACITY", kCandidateCapacity), limit / (2 * sizeof(uint32_t)));
			_pairs = gpu.CreateBuffer(_capacity * 2 * sizeof(uint32_t));

			_trace = gpu.CreateCompute(traceSpirv, traceWords, 5, sizeof(TracePush));
			_capture = gpu.CreateCompute(captureSpirv, captureWords, 5, sizeof(CapturePush));
		}

		void Execute(TraceResult& result)
		{
			const int launches = int(_input.launches->size());
			const int batches = (launches + _launchesPerBatch - 1) / _launchesPerBatch;
			result.segments.clear();
			result.launchFirst.clear();
			result.candidates.clear();
			for (int batch = 0; batch < batches; batch++)
			{
				const int first = batch * _launchesPerBatch;
				const int end = std::min(launches, first + _launchesPerBatch);
				std::cout << "\tGPU batch " << batch + 1 << "/" << batches << " (launches " << first << ".." << end - 1 << ")" << std::endl;
				std::vector<GpuSegment<Real>> records;
				TraceBatch(first, end, records);
				std::vector<uint32_t> pairs;
				CaptureBatch(uint32_t(records.size()), 0, uint32_t(_input.receivers->size()), pairs);
				Renumber(records, first, end, pairs, result);
			}
			result.launchFirst.push_back(int(result.segments.size()));
			std::sort(result.candidates.begin(), result.candidates.end());
		}

	private:
		// Traces the launches [first, end) level by level into `records`, in
		// the device's append order: the roots first, then every child
		// behind its level.
		void TraceBatch(int first, int end, std::vector<GpuSegment<Real>>& records)
		{
			SBR_ZONE_NAMED("GpuTrace");
			const int launches = end - first;
			std::vector<GpuSegment<Real>> roots(static_cast<size_t>(launches));
			for (int i = 0; i < launches; i++)
			{
				const int launch = first + i;
				const Vec3& origin = (*_input.transmitters)[size_t((*_input.launchTx)[size_t(launch)])].position;
				const Vec3& direction = (*_input.launches)[size_t(launch)];
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
				root.launch = launch;
				root.triangle = -1;
				root.coplanar = -1;
				root.material = -1;
				root.kind = kKindRoot;
				root.reflections = 0;
				root.transmissions = 0;
			}
			UploadVector(_gpu, _segments, roots);
			uint32_t count = uint32_t(launches);
			_gpu.Upload(_counter, &count, sizeof(count));

			uint32_t levelStart = 0;
			uint32_t levelEnd = count;
			for (int level = 0; levelEnd > levelStart; level++)
			{
				TracePush push{};
				push.pendingFirst = levelStart;
				push.pendingCount = levelEnd - levelStart;
				push.maxReflections = _input.maxReflections;
				push.maxTransmissions = _input.maxTransmissions;
				_gpu.Dispatch(_trace, { &_nodes, &_triangles, &_leaves, &_segments, &_counter }, &push, sizeof(push), push.pendingCount);
				std::cout << "\tGPU trace level " << level << " -> " << push.pendingCount << " segments" << std::endl;
				_gpu.Download(_counter, &count, sizeof(count));
				levelStart = levelEnd;
				levelEnd = count;
			}
			records.resize(count);
			_gpu.Download(_segments, records.data(), uint64_t(count) * sizeof(GpuSegment<Real>));
		}

		// Captures the receivers [first, end) against the batch's segments
		// and appends (receiver, local segment) pairs. A count above the
		// capacity reruns with the count, or splits the receivers when the
		// buffer may not grow that far.
		void CaptureBatch(uint32_t segmentCount, uint32_t first, uint32_t end, std::vector<uint32_t>& pairs)
		{
			SBR_ZONE_NAMED("GpuCapture");
			int reruns = 0;
			for (;;)
			{
				uint32_t zero = 0;
				_gpu.Upload(_counter, &zero, sizeof(zero));
				CapturePush push{};
				push.segmentCount = segmentCount;
				push.receiverFirst = first;
				push.receiverCount = end - first;
				push.capacity = uint32_t(_capacity);
				_gpu.Dispatch(_capture, { &_segments, &_receivers, &_slopes, &_pairs, &_counter }, &push, sizeof(push), segmentCount, end - first);
				uint32_t count = 0;
				_gpu.Download(_counter, &count, sizeof(count));
				if (count <= _capacity)
				{
					std::cout << "\tGPU capture -> " << count << " candidates (capacity " << _capacity << ", reruns " << reruns << ")" << std::endl;
					const size_t before = pairs.size();
					pairs.resize(before + size_t(count) * 2);
					if (count > 0) _gpu.Download(_pairs, pairs.data() + before, uint64_t(count) * 2 * sizeof(uint32_t));
					return;
				}
				const uint64_t needed = uint64_t(count) * 2 * sizeof(uint32_t);
				const DeviceInfo& info = _gpu.Info();
				if (needed <= std::min({ kBatchBytes, info.heapBytes / 4, info.storageRange, info.allocationLimit }))
				{
					_capacity = count;
					_pairs = Buffer();
					_pairs = _gpu.CreateBuffer(needed);
					reruns++;
					continue;
				}
				if (end - first < 2)
				{
					std::ostringstream message;
					message << "receiver " << first << " alone captures " << count << " segments, more than the device can hold; reduce the tessellation";
					throw std::runtime_error(message.str());
				}
				const uint32_t middle = first + (end - first) / 2;
				CaptureBatch(segmentCount, first, middle, pairs);
				CaptureBatch(segmentCount, middle, end, pairs);
				return;
			}
		}

		// Turns a batch's records, in append order, into the CPU's segments:
		// launch by launch, breadth first, the transmit child before the
		// reflect child, with the Ray fields as RayLaunch sets them. The
		// candidate pairs follow their segments to the new indices.
		void Renumber(const std::vector<GpuSegment<Real>>& records, int first, int end, const std::vector<uint32_t>& pairs, TraceResult& result)
		{
			const int launches = end - first;
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
			std::deque<int> queue;
			for (int i = 0; i < launches; i++)
			{
				result.launchFirst.push_back(int(result.segments.size()));
				queue.push_back(i);
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
			for (size_t p = 0; p + 1 < pairs.size(); p += 2)
			{
				result.candidates.push_back(std::make_pair(int(pairs[p]), newIndex[pairs[p + 1]]));
			}
		}

		const TracerInput& _input;
		GpuContext& _gpu;
		Buffer _nodes;
		Buffer _triangles;
		Buffer _leaves;
		Buffer _receivers;
		Buffer _slopes;
		Buffer _segments;
		Buffer _pairs;
		Buffer _counter;
		Pipeline _trace;
		Pipeline _capture;
		uint64_t _bound = 0;
		int _launchesPerBatch = 1;
		uint64_t _capacity = 0;
	};
}

void TraceAndCapture(const TracerInput& input, GpuContext& gpu, GpuPrecision precision, TraceResult& result)
{
	if (precision == GpuPrecision::Double)
	{
		Run<double> run(input, gpu, kSpirv_trace_double, kSpirvWords_trace_double, kSpirv_capture_double, kSpirvWords_capture_double);
		run.Execute(result);
	}
	else
	{
		Run<float> run(input, gpu, kSpirv_trace_single, kSpirvWords_trace_single, kSpirv_capture_single, kSpirvWords_capture_single);
		run.Execute(result);
	}
}

} // namespace gpu
} // namespace sbr
