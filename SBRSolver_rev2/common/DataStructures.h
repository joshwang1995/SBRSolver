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
	double distance = INF;
};

// Material property for material at single frequency
struct MaterialProperties
{
	double frequency;
	double transmissionLoss;
	double reflectionLoss;
	double relPermittivityRe;
	double relPermittivityIm;
	double relConductivity;
};
using Materials = std::vector<MaterialProperties>;

// Antenna gain information
struct GainInfo
{
    float data[4]; //theta gain, theta phase, phi gain, phi phase
};
using ThetaPhiAngle = std::pair<float, float>; 
using GainMap = std::map<ThetaPhiAngle, GainInfo>; // map theta phi angles to gain