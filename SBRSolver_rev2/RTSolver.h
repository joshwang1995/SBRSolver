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
#include <mutex>


class RTSolver
{
public:
	RTSolver();
	virtual ~RTSolver();
	int Init(MaterialProperties materialProperties[], int materialPropertiesCount, BVH<Triangle>& triangles);
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
		VecVec3& receivers,
		std::vector<Triangle*>& triangleMesh,
		int txTesslation
	);
	bool SavePathsAsVtk(std::string fname);
	bool SaveReceiversAsVtk(std::string fname);
	bool SaveFieldAsCsv(std::string fname);
	bool SaveIcosahedronAsVtk(std::string fname, Vec3 rayOrg, int tessellation);
	void CmdLineDebug();
protected:
	MaterialProperties* _materialProperties = nullptr;
	Paths** _rayPaths;
	BVH<Triangle>* _bvh;
	std::vector<Vec3>* _shootRayList;
	VecVec3* _receivers;
	std::vector<Vec3c>* _efield;

	bool _isSceneInitialized;
	int _txTesslation;
	int _maxTransmissionCount;
	int _maxReflectionCount;
	int _pathsCount;
	int _receiverCount;
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
	void ImagePathCorrection(int receiverId, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh);
	void MultiPathCorrection(std::vector<Ray>& multiPath, int receiverId, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh);
	void DirectPathCorrection(std::vector<Ray>& directPath, int receiverId, const Vec3& txPoint);
	void RemoveDuplicatePath(const Vec3& receiver, int receiverId);
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
double GetTotalRayLength(const std::vector<Ray>& rayPaths);
double DistanceToReceiver(const std::vector<Ray>& rayPaths, const Vec3& receiver);
bool sortbysec(const std::pair<int, int>& a, const std::pair<int, int>& b);
Vec3 ImagePoint(const Triangle& triangle, const Vec3& src);
std::vector<std::pair<int, Vec3>> GetImageSources
(
	const std::vector<Ray>& path, 
	const Vec3& txPoint, 
	const std::vector<Triangle*>& triangleMesh
);