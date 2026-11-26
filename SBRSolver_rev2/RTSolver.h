#pragma once

#include "Ray.h"
#include "Tracer.h"
#include "timer/Timer.h"
#include "common/VecMatDef.h"
#include "common/Constants.h"
#include "common/DataStructures.h"
#include "Icosahedron.h"
#include "Triangle.h"
#include "Preprocessor.h"
#include "FieldCompute.h"
#include "viewer/SceneData.h"
#include <mutex>

class RTSolver
{
public:
	RTSolver();
	virtual ~RTSolver();
	int Init(MaterialProperties materialProperties[], int materialPropertiesCount, BVH<Triangle>& triangles);
	// Threads used by ExecuteRayTracing (0 = one per hardware thread).
	void SetThreadCount(int threads) { _threadCount = threads; }
	void Cleanup();
	// Launches the rays of every transmitter together, finds the paths to
	// every receiver and sums the fields of all transmitters at each one,
	// as phasors or, with powerSum, as powers. Returns the number of
	// candidate paths before the image-method correction.
	size_t ExecuteRayTracing
	(
		const std::vector<TransmitterSource>& transmitters,
		int maxReflectionCount,
		int maxTransmissionCount,
		double frequency,
		bool useFresnelCoeff,
		bool powerSum,
		VecVec3& receivers,
		std::vector<Triangle*>& triangleMesh
	);
	// Results of the last ExecuteRayTracing for the viewer: the total field
	// per receiver and every image-corrected path with its own field and
	// its transmitter.
	const std::vector<Vec3c>& Fields() const;
	std::vector<sbr::viewer::ScenePath> CollectPaths() const;
	bool SavePathsAsVtk(std::string fname);
	bool SaveReceiversAsVtk(std::string fname);
	bool SaveFieldAsCsv(std::string fname);
	// The launch directions of every transmitter, as lines from its position.
	bool SaveIcosahedronAsVtk(std::string fname);
	void CmdLineDebug();
protected:
	MaterialProperties* _materialProperties = nullptr;
	// The paths of every receiver, in launch order and then in the order of
	// the ray tree; the image-method correction and the duplicate removal
	// prune them in place.
	std::vector<std::vector<ReceiverPath>> _receiverPaths;
	BVH<Triangle>* _bvh = nullptr;
	std::vector<Vec3>* _shootRayList = nullptr;
	VecVec3* _receivers = nullptr;
	std::vector<Vec3c>* _efield = nullptr;
	// Field of every captured path, per receiver, in the order CollectPaths
	// visits them (launch, then path). Grows with the captured paths only.
	std::vector<std::vector<Vec3c>> _pathFields;

	// The transmitters of the last run and the transmitter of every launch
	// direction in _shootRayList.
	std::vector<TransmitterSource> _transmitters;
	std::vector<int> _launchTx;

	bool _isSceneInitialized = false;
	int _maxTransmissionCount = 0;
	int _maxReflectionCount = 0;
	int _receiverCount = 0;
	int _threadCount = 0;
	// Triangles of every coplanar surface, indexed by Triangle::coplanarId.
	std::vector<std::vector<Triangle*>> _coplanarGroups;
private:
	// Replaces every captured path of a receiver by the exact path with the
	// same interaction sequence (image method) from its own transmitter and
	// drops the ones that do not exist: no intersection with the surface, or
	// a segment blocked by another surface.
	void ImagePathCorrection(int receiverId, const std::vector<Triangle*>& triangleMesh);
	bool MultiPathCorrection(std::vector<Ray>& multiPath, int receiverId, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh);
	bool DirectPathCorrection(std::vector<Ray>& directPath, int receiverId, const Vec3& txPoint);
	// True when the open segment from `from` to `to` crosses a triangle that
	// belongs to neither of the two coplanar surfaces (-1: none) at its ends.
	bool SegmentBlocked(const Vec3& from, const Vec3& to, int skipCoplanarA, int skipCoplanarB) const;
	// Intersection of the ray from `origin` towards `target` with a coplanar
	// surface; returns the triangle hit (or nullptr) and the point.
	Triangle* IntersectSurface(int coplanarId, const Vec3& origin, const Vec3& target, Vec3& point) const;
	void RemoveDuplicatePath(int receiverId);
	// One path of _receiverPaths, as the outputs address it.
	struct PathIndex
	{
		int launch;
		int receiver;
		size_t position;
	};
	// Every path of the run in the order the outputs list them: by launch,
	// then by receiver, then in the receiver's own order.
	std::vector<PathIndex> OutputOrder() const;
};

void PathToVector
(
	const std::vector<Ray>& path,
	std::vector<Vec3>& vertices,
	std::vector<std::vector<int>>& lineIndex,
	std::vector<int>& launchIds,
	std::vector<int>& receiverIds,
	int launchId,
	int receiverId,
	int& totalLines
);

std::vector<int> GetHitSurfaceIds(const std::vector<Ray>& rayPaths);
std::vector<int> GetHitCoplanarIds(const std::vector<Ray>& rayPaths);
double GetTotalRayLength(const std::vector<Ray>& rayPaths);
double DistanceToReceiver(const std::vector<Ray>& rayPaths, const Vec3& receiver);
Vec3 ImagePoint(const Triangle& triangle, const Vec3& src);
std::vector<std::pair<int, Vec3>> GetImageSources
(
	const std::vector<Ray>& path, 
	const Vec3& txPoint, 
	const std::vector<Triangle*>& triangleMesh
);