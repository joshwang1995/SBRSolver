#pragma once

#include <utility>
#include <vector>

#include "BVH/BVH.h"
#include "Ray.h"
#include "Triangle.h"
#include "common/DataStructures.h"
#include "common/VecMatDef.h"

// One traced ray segment of a run. The segments of a launch are stored in
// the order of the old tree walk (breadth first, transmit child before
// reflect child), launch after launch, so that a receiver's candidate paths
// come out in the order the duplicate removal depends on.
struct Segment
{
	Ray ray;          // points, materials, hit ids and unfolded length as RayLaunch set them
	int parent = -1;  // index of the parent segment, -1 for the root of a launch
	int launch = 0;   // index into the launch list
};

// Everything a tracer reads. The pointers stay valid for the run.
struct TracerInput
{
	BVH<Triangle>* bvh = nullptr;
	const std::vector<TransmitterSource>* transmitters = nullptr;
	const std::vector<Vec3>* launches = nullptr;  // unit direction of every launch
	const std::vector<int>* launchTx = nullptr;   // transmitter of every launch
	const VecVec3* receivers = nullptr;
	std::vector<double> tanBeta;                  // reception sphere slope per transmitter
	int maxReflections = 0;
	int maxTransmissions = 0;
	int threads = 1;
};

// What a tracer produces: every segment, and every (receiver, segment) pair
// whose reception sphere the segment crosses, sorted by receiver, then by
// segment index.
struct TraceResult
{
	std::vector<Segment> segments;
	std::vector<int> launchFirst;  // first segment of every launch, plus the end
	std::vector<std::pair<int, int>> candidates;
};
