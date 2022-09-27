#pragma once
#include <complex>
#include "Ray.h"
#include "common/Constants.h"
#include "common/DataStructures.h"
#include "common/VecMatDef.h"

class FieldCompute
{
public:
	FieldCompute();
	FieldCompute(Paths** rayPaths, std::vector<Triangle*>& triangleMesh);
	~FieldCompute();
private:
	Vec3c ComputeRefcField
	(
		const Vec3c& efield_i_sph,
		double rel_perm,
		double sigma,
		double freq,
		double theta_i,
		double width,
		bool inf_wall
	);

	Vec3c ComputeTransField
	(
		const Vec3c& efield_i_sph,
		double rel_perm,
		double sigma,
		double freq,
		double theta_i,
		double width,
		bool inf_wall
	);

	void GetTECoeff
	(
		cdouble theta_i, cdouble theta_t,
		cdouble rel_perm_i, cdouble rel_perm_t,
		double lamda, double width, bool inf_wall,
		cdouble& ref_coeff, cdouble& tran_coeff
	);

	void GetTMCoeff
	(
		cdouble theta_i, cdouble theta_t,
		cdouble rel_perm_i, cdouble rel_perm_t,
		double lamda, double width, bool inf_wall,
		cdouble& ref_coeff, cdouble& tran_coeff
	);
	cdouble GetTransAngle(double thetaIncident, cdouble epsilonIncident, cdouble epsilonTransmit);
};

