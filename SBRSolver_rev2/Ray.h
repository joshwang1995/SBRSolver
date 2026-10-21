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
	Vec3 sourcePoint = Vec3(INF,INF,INF);
	Vec3 targetPoint = Vec3(INF,INF,INF);
	int penetrationMaterialId = -1;
	int reflectionMaterialId = -1;
	int hitSurfaceID = -1;
	int hitCoplanarId = -1;
	double angleFromSurfaceNormal = INF;
	double pathLength = INF;
	bool captured = false;
	int id = -1;
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
	PathTreeNode* childTransmit = nullptr;
	PathTreeNode* childReflect = nullptr;
};
PathTreeNode* newPathTreeNode
(
	Vec3 sourcePoint,
	Vec3 targetPoint,
	int penetrationMaterialId,
	int reflectionMaterialId,
	int surfaceId,
	int coplanarId,
	double pathLength,
	double angleFromSurfaceNormal
);
PathTreeNode* DeleteChildNodes(PathTreeNode* node);
PathTreeNode* CloneNode(PathTreeNode* node);
void CloneTree(PathTreeNode* orgTree, PathTreeNode* cloneTree);

class Paths
{
public:
	std::vector<std::vector<Ray>> rayPaths;
	Paths();
	Paths(TreeNode* rootNode);
	Paths(PathTreeNode* rootNode);
	~Paths();
};