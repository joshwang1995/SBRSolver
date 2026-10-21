#pragma once

#include "Ray.h"
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
	int ExecuteRayTracing
	(
		Vec3 sourcePoint,
		const Mat3& txCoordSys,
		int maxReflectionCount,
		int maxTransmissionCount,
		double frequency,
		double pt,
		bool useFresnelCoeff,
		int antennaType,
		VecVec3& receivers,
		std::vector<Triangle*>& triangleMesh,
		int txTesslation
	);
	// Results of the last ExecuteRayTracing for the viewer: the total field
	// per receiver and every image-corrected path with its own field.
	const std::vector<Vec3c>& Fields() const;
	std::vector<sbr::viewer::ScenePath> CollectPaths() const;
	bool SavePathsAsVtk(std::string fname);
	bool SaveReceiversAsVtk(std::string fname);
	bool SaveFieldAsCsv(std::string fname);
	bool SaveIcosahedronAsVtk(std::string fname, Vec3 rayOrg, int tessellation);
	void CmdLineDebug();
protected:
	MaterialProperties* _materialProperties = nullptr;
	Paths** _rayPaths = nullptr;
	BVH<Triangle>* _bvh = nullptr;
	std::vector<Vec3>* _shootRayList = nullptr;
	VecVec3* _receivers = nullptr;
	std::vector<Vec3c>* _efield = nullptr;

	bool _isSceneInitialized = false;
	int _txTesslation = 0;
	int _maxTransmissionCount = 0;
	int _maxReflectionCount = 0;
	int _pathsCount = 0;
	int _receiverCount = 0;
	int _threadCount = 0;
	// Triangles of every coplanar surface, indexed by Triangle::coplanarId.
	std::vector<std::vector<Triangle*>> _coplanarGroups;
private:
	void RayLaunch
	(
		PathTreeNode* rayTreeNode, 
		Vec3& sourcePoint, 
		Vec3& directionPoint, 
		int transMaterialID,
		int refMaterialID,
		int reflectionCnt,
		int transmissionCnt,
		double totalPathLength,
		bool isRoot, 
		double lastAnglefromN
	);
	void RayCapture
	(
		PathTreeNode* rayTreeNode,
		const Vec3& receiver, 
		double totalPathLength
	);
	// Replaces every captured path of a receiver by the exact path with the
	// same interaction sequence (image method) and drops the ones that do not
	// exist: no intersection with the surface, or a segment blocked by another
	// surface.
	void ImagePathCorrection(int receiverId, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh);
	bool MultiPathCorrection(std::vector<Ray>& multiPath, int receiverId, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh);
	bool DirectPathCorrection(std::vector<Ray>& directPath, int receiverId, const Vec3& txPoint);
	// True when the open segment from `from` to `to` crosses a triangle that
	// belongs to neither of the two coplanar surfaces (-1: none) at its ends.
	bool SegmentBlocked(const Vec3& from, const Vec3& to, int skipCoplanarA, int skipCoplanarB) const;
	// Intersection of the ray from `origin` towards `target` with a coplanar
	// surface; returns the triangle hit (or nullptr) and the point.
	Triangle* IntersectSurface(int coplanarId, const Vec3& origin, const Vec3& target, Vec3& point) const;
	void RemoveDuplicatePath(int receiverId);
	void InitRayPaths();
	void DeleteRayPaths();
};

void PathsToVector
(
	const std::vector<std::vector<Ray>>& rayPaths,
	std::vector<Vec3>& vertices,
	std::vector<std::vector<int>>& lineIndex,
	std::vector<int>& launchIds,
	std::vector<int>& receiverIds,
	int launchId,
	int receiverId,
	int& totalLines
);

bool HitReceptionSphere
(
	const Vec3& sourcePoint,
	const Vec3& targetPoint,
	const Vec3& sphereCenter,
	double pathLength,
	int txTesslation,
	Vec3& capturePoint
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