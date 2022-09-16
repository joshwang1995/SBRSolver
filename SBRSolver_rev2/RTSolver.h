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
	int ExecuteRayTracing
	(
		Vec3 sourcePoint, 
		int maxReflectionCount, 
		int maxTransmissionCount, 
		VecVec3 receivers,
		int txTesslation = 0
	);
	bool SavePathsAsVtk(std::string fname);
	bool SaveIcosahedronAsVtk(std::string fname, Vec3 rayOrg, int tessellation);
	void CmdLineDebug();
protected:
	bool _isSceneInitialized;
	MaterialProperties* _materialProperties = nullptr;
	int _txTesslation;
	int _maxTransmissionCount;
	int _maxReflectionCount;
	Paths* _capturedPaths;
	Paths* _rayPaths;
	int _pathsCount;
	BVH<Triangle>* _bvh;
	std::vector<Vec3>* _shootRayList;
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

	void RayLaunchAndCapture
	(
		PathTreeNode* rayTreeNode, 
		Vec3& sourcePoint, 
		Vec3& directionPoint, 
		const Vec3& receiver, 
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
};

bool HitReceptionSphere
(
	const Vec3& sourcePoint,
	const Vec3& targetPoint,
	const Vec3& sphereCenter,
	double pathLength,
	int txTesslation,
	Vec3& capturePoint
);