#pragma once

#include "VecMatDef.h"
#include "Constants.h"
#include <map>
#include <vector>

// Information to stroe after an intersection
struct HitInfo
{
	Vec3 normal;
	Vec3 pointIntersect;
	int materialId = -1;
	int triangelId = -1;
	int coplanarId = -1;
	double distance = INF;
};

// Material property for material at single frequency
struct MaterialProperties
{
	double frequency = 0.0;
	double transmissionLoss = 0.0;
	double reflectionLoss = 0.0;
	double relPermittivityRe = 0.0;
	double relPermittivityIm = 0.0;
	double relConductivity = 0.0;
	double width = 0.0;
};
using Materials = std::vector<MaterialProperties>;

// Antenna gain information
struct GainInfo
{
    float data[4]; //theta gain, theta phase, phi gain, phi phase
};
using ThetaPhiAngle = std::pair<float, float>; 
using GainMap = std::map<ThetaPhiAngle, GainInfo>; // map theta phi angles to gain