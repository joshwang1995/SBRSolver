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
	double frequency;
	double transmissionLoss;
	double reflectionLoss;
	double relPermittivityRe;
	double relPermittivityIm;
	double relConductivity;
};
using Materials = std::vector<MaterialProperties>;

// A transmitter as the ray tracer sees it: where its rays start, how the
// antenna is oriented and how many rays are launched (10 * 4^tessellation + 2).
struct TransmitterSource
{
	Vec3 position = Vec3::Zero();
	Mat3 coordSys = Mat3::Identity(); // rows: antenna x, y, z axes in global coordinates
	double powerWatts = 1.0;
	int antennaType = 1;              // see AntennaType in SimulationConfig.h
	int tessellation = 2;
};

// Antenna gain information
struct GainInfo
{
    float data[4]; //theta gain, theta phase, phi gain, phi phase
};
using ThetaPhiAngle = std::pair<float, float>; 
using GainMap = std::map<ThetaPhiAngle, GainInfo>; // map theta phi angles to gain