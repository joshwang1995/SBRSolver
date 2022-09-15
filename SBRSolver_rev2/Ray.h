#pragma once

#include "common/VecMatDef.h"
#include "common/DataStructures.h"
#include "Triangle.h"
#include "BVH/BVH.h"
#include <vector>
#include <queue>
#include <tuple>
#include <list>
#include <atomic>

extern std::atomic<int> RayGlobalId;

struct Ray
{
	Vec3 sourcePoint;
	Vec3 targetPoint;
	int penetrationMaterialId;
	int reflectionMaterialId;
	double angleFromSurfaceNormal;
	double pathLength;
	bool captured = false;
	int id;
};

class TreeNode
{
public:
	Ray ray;
	std::list<TreeNode> child;
	TreeNode();
	~TreeNode();
};

struct PathTreeNode
{
	Ray ray;
	PathTreeNode* childDirect;
	PathTreeNode* childReflect;
};
PathTreeNode* newPathTreeNode
(
	Vec3 SourcePoint,
	Vec3 TargetPoint,
	int PenetrationMaterialId,
	int ReflectionMaterialId,
	double PathLength,
	double AngleFromSurfaceNormal
);

class Paths
{
public:
	std::vector<std::vector<Ray>> rayPaths;
	Paths();
	Paths(TreeNode* rootNode);
	Paths(PathTreeNode* rootNode);
	~Paths();
};