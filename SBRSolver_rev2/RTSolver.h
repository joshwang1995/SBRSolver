#pragma once

#include "Ray.h"
#include "timer/Timer.h"
#include "common/VecMatDef.h"
#include "common/Constants.h"
#include "common/DataStructures.h"
#include "Icosahedron.h"
#include "Triangle.h"
#include "Preprocessor.h"
#include <mutex>


class RTSolver
{
public:
	RTSolver();
	virtual ~RTSolver();
	int Init(MaterialProperties materialProperties[], int materialPropertiesCount, BVH<Triangle>& triangles);
	void Cleanup();
	int ExecuteRayTracing(Vec3 sourcePoint, int maxBounceCount, double maxPathLoss, int txTesslation = 0);
	bool SavePathsAsVtk(std::string fname);
	void DebugFunc();
protected:
	bool _isSceneInitialized;
	MaterialProperties* _materialProperties = nullptr;
	int _maxPenetrationCount;
	int _maxBounceCount;
	double _maxPathLoss;
	Paths* _rayPaths;
	int _pathsCount;
	BVH<Triangle>* _bvh;
private:
	std::mutex _mutex;
	void RayLaunch(PathTreeNode* rayTreeNode, Vec3& sourcePoint, Vec3& directionPoint, int transmissionLossId, int reflectionMaterialId, int bounceCnt, int penetrationCnt, double totalLoss, bool isRoot, float lastAnglefromN);
	bool RayLaunch2(Vec3& rayOrig, Vec3& rayDir, double totalPathLen, int bounceCnt, int penetrationCnt);
};