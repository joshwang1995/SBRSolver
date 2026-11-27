#include "CpuTracer.h"

#include <algorithm>
#include <cmath>
#include <deque>

#include "ReceiverGrid.h"
#include "common/Constants.h"
#include "common/ParallelFor.h"
#include "common/Profiling.h"

namespace
{
	// A ray that waits for its intersection test: the arguments of one call
	// of the old recursive RayLaunch.
	struct PendingRay
	{
		Vec3 source;
		Vec3 direction;
		int transMaterialId;
		int refMaterialId;
		int reflectionCount;
		int transmissionCount;
		double totalPathLength;
		double lastAngleFromNormal;
		int parent;
	};

	// A ray that leaves the scene keeps going this far, so that receivers far
	// away can still capture it.
	const double kEscapeLength = 1E10;

	// Segments per chunk of the capture loop; every chunk collects its own
	// candidates, which keeps the workers off a shared list.
	const int kCaptureChunk = 1024;

	// Traces one launch into `out`, in breadth-first order with the transmit
	// child before the reflect child. Parent indices are local to `out`.
	void TraceLaunch(const TracerInput& input, int launch, std::vector<Segment>& out)
	{
		const TransmitterSource& tx = (*input.transmitters)[(*input.launchTx)[launch]];
		std::deque<PendingRay> pending;
		pending.push_back({ tx.position, (*input.launches)[launch], -1, -1, 0, 0, 0.0, INF, -1 });
		while (!pending.empty())
		{
			const PendingRay p = pending.front();
			pending.pop_front();

			Segment segment;
			segment.parent = p.parent;
			segment.launch = launch;
			Ray& ray = segment.ray;
			HitInfo hit;
			if (input.bvh->RayIntersects(0, p.source, p.direction, hit))
			{
				const double angleFromNormal = AngleBetween(hit.normal, p.direction);
				ray.sourcePoint = p.source;
				ray.targetPoint = hit.pointIntersect;
				ray.penetrationMaterialId = p.transMaterialId;
				ray.reflectionMaterialId = p.refMaterialId;
				ray.hitSurfaceID = hit.triangelId;
				ray.hitCoplanarId = hit.coplanarId;
				ray.pathLength = p.totalPathLength + hit.distance;
				ray.angleFromSurfaceNormal = p.lastAngleFromNormal;
				const int index = int(out.size());
				out.push_back(segment);

				if (p.transmissionCount < input.maxTransmissions)
				{
					pending.push_back({ hit.pointIntersect, p.direction, hit.materialId, -1, p.reflectionCount, p.transmissionCount + 1, ray.pathLength, angleFromNormal, index });
				}
				if (p.reflectionCount < input.maxReflections)
				{
					const Vec3 reflected = Reflect(p.direction, hit.normal).normalized();
					pending.push_back({ hit.pointIntersect, reflected, -1, hit.materialId, p.reflectionCount + 1, p.transmissionCount, ray.pathLength, angleFromNormal, index });
				}
			}
			else
			{
				ray.sourcePoint = p.source;
				ray.targetPoint = p.source + (p.direction * kEscapeLength);
				ray.penetrationMaterialId = p.transMaterialId;
				ray.reflectionMaterialId = p.refMaterialId;
				ray.hitSurfaceID = -1;
				ray.hitCoplanarId = -1;
				ray.pathLength = INF;
				ray.angleFromSurfaceNormal = p.lastAngleFromNormal;
				out.push_back(segment);
			}
		}
	}
}

void TraceSegments(const TracerInput& input, TraceResult& result)
{
	const int launches = int(input.launches->size());
	std::vector<std::vector<Segment>> perLaunch(launches);
	ParallelFor(launches, input.threads, [&](int i)
	{
		SBR_ZONE_NAMED("RayLaunch");
		TraceLaunch(input, i, perLaunch[i]);
	});

	result.launchFirst.assign(launches + 1, 0);
	for (int i = 0; i < launches; i++)
	{
		result.launchFirst[i + 1] = result.launchFirst[i] + int(perLaunch[i].size());
	}
	result.segments.clear();
	result.segments.reserve(result.launchFirst[launches]);
	for (int i = 0; i < launches; i++)
	{
		for (Segment& segment : perLaunch[i])
		{
			if (segment.parent >= 0) segment.parent += result.launchFirst[i];
			result.segments.push_back(segment);
		}
	}
}

void CaptureCandidates(const TracerInput& input, TraceResult& result)
{
	SBR_ZONE_NAMED("CaptureCandidates");
	const VecVec3& receivers = *input.receivers;
	result.candidates.clear();
	if (receivers.empty() || result.segments.empty()) return;

	const ReceiverGrid grid(receivers);
	Vec3 low = receivers[0];
	Vec3 high = receivers[0];
	for (const Vec3& r : receivers)
	{
		low = low.cwiseMin(r);
		high = high.cwiseMax(r);
	}
	const Vec3 centre = 0.5 * (low + high);
	const double rho = 0.5 * (high - low).norm();

	const int segmentCount = int(result.segments.size());
	const int chunks = (segmentCount + kCaptureChunk - 1) / kCaptureChunk;
	std::vector<std::vector<std::pair<int, int>>> perChunk(chunks);
	ParallelFor(chunks, input.threads, [&](int chunk)
	{
		std::vector<int> near;
		const int end = std::min(segmentCount, (chunk + 1) * kCaptureChunk);
		for (int s = chunk * kCaptureChunk; s < end; s++)
		{
			const Segment& segment = result.segments[s];
			const Ray& ray = segment.ray;
			const double unfoldedStart = (segment.parent < 0) ? 0.0 : result.segments[segment.parent].ray.pathLength;
			const double tanBeta = input.tanBeta[(*input.launchTx)[segment.launch]];

			// A point at distance s along the segment can capture a receiver
			// only while the sphere (u0 + s) tan(beta) still reaches the
			// receivers' bounding sphere: |c - P(s)| >= s - |S - c| - rho.
			const Vec3 direction = ray.targetPoint - ray.sourcePoint;
			const double length = direction.norm();
			const double sMax = ((ray.sourcePoint - centre).norm() + rho + unfoldedStart * tanBeta) / (1.0 - tanBeta);
			const double sEnd = std::min(length, sMax);
			const Vec3 end = ray.sourcePoint + direction * (sEnd / length);
			const double radiusEnd = (unfoldedStart + sEnd) * tanBeta;
			Mat23 box;
			box.row(0) = (ray.sourcePoint.cwiseMin(end).array() - radiusEnd).transpose();
			box.row(1) = (ray.sourcePoint.cwiseMax(end).array() + radiusEnd).transpose();

			near.clear();
			grid.Query(box, near);
			for (int r : near)
			{
				Vec3 capturePoint;
				if (HitReceptionSphere(ray.sourcePoint, ray.targetPoint, receivers[r], unfoldedStart, tanBeta, capturePoint))
				{
					perChunk[chunk].push_back(std::make_pair(r, s));
				}
			}
		}
	});

	for (const std::vector<std::pair<int, int>>& list : perChunk)
	{
		result.candidates.insert(result.candidates.end(), list.begin(), list.end());
	}
	std::sort(result.candidates.begin(), result.candidates.end());
}

void TraceAndCapture(const TracerInput& input, TraceResult& result)
{
	TraceSegments(input, result);
	CaptureCandidates(input, result);
}

double ReceptionSlope(int tessellation)
{
	// Reference: 191115_Ray_Launching.pdf
	const double beta = 1.0 / (3.0 * tessellation) * acos(-1.0 / sqrt(5.0));
	return tan(beta);
}

bool HitReceptionSphere
(
	const Vec3& sourcePoint,
	const Vec3& targetPoint,
	const Vec3& sphereCenter,
	double pathLength,
	double tanBeta,
	Vec3& capturePoint
)
{
	// Find the closest distance between the ray and the sphere center
	// Reference: https://www.geometrictools.com/Documentation/DistancePointLine.pdf
	Vec3 rayDir = targetPoint - sourcePoint; // Do not normalize, this is a line segment
	double t0 = rayDir.dot(sphereCenter - sourcePoint) / rayDir.squaredNorm();

	if (t0 < 0)
	{
		// Projection is negative, receiver is in the opposite direction of the ray direction
		return false;
	}
	else if (t0 >= 1)
	{
		// We can choose to clamp t0 to 1 here. This means that the closest point is the target point
		// However, I choose to ignore it because the reception sphere algorithm checks to see if the capture
		// point is in LOS of the receiver.
		t0 = 1; // Clamp t0 to 1 because the closest point is the target point
	}

	// Takahiro's method: the sphere grows with the unfolded path length from
	// the transmitter to the closest point of this segment.
	// Note: the unfolded path length is the path length from previous ray segments plus the distance
	//  the point of intersection and the receiver
	double unfoldedLength = pathLength + (t0 * rayDir).norm();
	double sphereRadius = unfoldedLength * tanBeta;

	// Rappaport's approximation
	// double unfoldedLength = pathLength + (sourcePoint + t0 * rayDir).norm();
	// double alpha = (69.0 * PI / 180.0) / txTesslation;
	// double sphereRadius = alpha * unfoldedLength / sqrt(3.0);

	// double unfoldedLength = pathLength + (t0 * rayDir).norm();
	// double phi = (sqrt(5) + 1) / 2;
	// double alpha = acos(phi / (phi * phi + 1)) / (pow(2, txTesslation - 1)) / sqrt(3);
	// double sphereRadius = alpha * unfoldedLength;


	double distance = (sphereCenter - (sourcePoint + t0 * rayDir)).norm();
	if (distance <= sphereRadius)
	{
		// Distance from closest point on the ray to sphere center is smaller
		// than the sphere radius, therefore it is not captured
		capturePoint = sourcePoint + t0 * rayDir;
		return true;
	}
	return false;
}
