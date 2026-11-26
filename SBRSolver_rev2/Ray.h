#pragma once

#include "common/VecMatDef.h"
#include "common/DataStructures.h"
#include "Triangle.h"
#include "BVH/BVH.h"
#include <vector>

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

class Paths
{
public:
	std::vector<std::vector<Ray>> rayPaths;
	Paths();
	~Paths();
};
