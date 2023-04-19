#pragma once
#include <complex>
#include "Ray.h"
#include "common/Global.h"
#include "common/DataStructures.h"
#include "common/VecMatDef.h"
#include "RTSolver.h"

class FieldCompute
{
public:
	FieldCompute();
	FieldCompute
	(
		Paths** rayPaths,
		std::vector<Triangle*>& triangleMesh,
		VecVec3& receivers,
		int pathsCount,
		MaterialProperties* materials,
		Mat3 txCoordSys,
		double frequency,
		double pt,
		bool useFresnelCoeff
	)
		:
		_rayPaths(rayPaths),
		_triangleMesh(&triangleMesh),
		_receivers(&receivers),
		_pathsCount(pathsCount),
		_materials(materials),
		txCoordSys(txCoordSys),
		frequency(frequency),
		txPower(pt),
		useFresnelCoeff(useFresnelCoeff)
	{
		lamda = SPEED_OF_LIGHT / frequency;
		k = (2 * PI) / lamda;
	};
	~FieldCompute();

	const Mat3 globalCoordSys = Mat3::Identity();

	double lamda = 0.0;
	double k = 0.0;
	Mat3 txCoordSys = Mat3();
	double frequency = 0.0;
	double txPower = 0.0;
	bool useFresnelCoeff = false;
	

	Vec3c FieldAtReceiver(int receiverId);
	void RefCoeffTest(int numPts, double freq, cdouble epsilon1, cdouble epsilon2, std::string fname);
	void TransCoeffTest(int numPts, double freq, cdouble epsilon1, cdouble epsilon2, std::string fname);

	Paths** _rayPaths;
	std::vector<Triangle*>* _triangleMesh;
	int _pathsCount;
	VecVec3* _receivers;
	MaterialProperties* _materials;
	Vec3c FieldForPath(const std::vector<Ray>& path);

	Vec3c ComputeRefcField
	(
		const Vec3& vecGlobal,
		const Vec3c& efield_i_sph,
		int materialId,
		const Mat3& currentCoordSys,
		const Mat3& nextCoordSys
	);

	Vec3c ComputeTransField
	(
		const Vec3& vecGlobal,
		const Vec3c& efield_i_sph,
		int materialId,
		const Mat3& currentCoordSys,
		const Mat3& nextCoordSys
	);

	void GetTECoeff
	(
		cdouble theta_i, cdouble theta_t,
		cdouble rel_perm_i, cdouble rel_perm_t,
		double width, bool inf_wall,
		cdouble& ref_coeff, cdouble& tran_coeff
	);

	void GetTMCoeff
	(
		cdouble theta_i, cdouble theta_t,
		cdouble rel_perm_i, cdouble rel_perm_t,
		double width, bool inf_wall,
		cdouble& ref_coeff, cdouble& tran_coeff
	);
	cdouble GetTransAngle(double thetaIncident, cdouble epsilonIncident, cdouble epsilonTransmit);
	Mat3 GetSurfCoordSys(const int& hitSurfId, const Ray& rayIncident);
	Vec3c GetAnalyticEfieldPattern(int antennaType, double theta, double phi, double pt);
};