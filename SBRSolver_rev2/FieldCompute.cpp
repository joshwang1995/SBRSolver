#include "FieldCompute.h"

FieldCompute::FieldCompute()
{
}

FieldCompute::FieldCompute(Paths** rayPaths, std::vector<Triangle*>& triangleMesh)
{
}

FieldCompute::~FieldCompute()
{
}

Vec3c FieldCompute::ComputeRefcField(const Vec3c& efield_i_sph, double rel_perm, double sigma, double freq, double theta_i, double width, bool inf_wall)
{
	cdouble epsilonTransmit;
	epsilonTransmit.real(rel_perm);
	epsilonTransmit.imag(sigma / (E0 * 2 * PI * freq));

	double lamda = SPEED_OF_LIGHT / freq;

	// Permittivity of air
	cdouble epsilonIncident(1, 0);

	cdouble theta_t = trans_angle(theta_i, epsilonIncident, epsilonTransmit);

	cdouble te_refc, tm_refc, te_tranc, tm_tranc;
	GetTMCoeff(theta_i, theta_t, epsilonIncident, epsilonTransmit, lamda, width, inf_wall, tm_refc, tm_tranc);
	GetTECoeff(theta_i, theta_t, epsilonIncident, epsilonTransmit, lamda, width, inf_wall, te_refc, te_tranc);

	// x is r, y is theta, z is phi
	cdouble Er = efield_i_sph[0];
	cdouble Etheta = efield_i_sph[1];
	cdouble Ephi = efield_i_sph[2];

	Vec3c refc_field (Er, Etheta * efield_i_sph.y() * te_refc); 
	refc_field = ;
	refc_field.phi = ;
	refc_field.theta = efield_i_sph.theta * tm_refc;

	return refc_field;
}

Vec3c FieldCompute::ComputeTransField(const Vec3c& efield_i_sph, double rel_perm, double sigma, double freq, double theta_i, double width, bool inf_wall)
{
	complex<double> e_t;
	e_t.real(rel_perm);
	e_t.imag(sigma / (E0 * 2 * PI * freq));

	double lamda = SPEED_OF_LIGHT / freq;

	// Permittivity of air
	complex<double> e_i(1, 0);

	complex<double> theta_t = trans_angle(theta_i, e_i, e_t);

	complex<double> te_refc, tm_refc, te_tranc, tm_tranc;
	tm_coeff(theta_i, theta_t, e_i, e_t, lamda, width, inf_wall, tm_refc, tm_tranc);
	te_coeff(theta_i, theta_t, e_i, e_t, lamda, width, inf_wall, te_refc, te_tranc);

	Cvect3dsph trans_field;
	trans_field.r = efield_i_sph.r;
	trans_field.phi = efield_i_sph.phi * te_tranc;
	trans_field.theta = efield_i_sph.theta * tm_tranc;

	return trans_field;
}

void FieldCompute::GetTECoeff(cdouble theta_i, cdouble theta_t, cdouble rel_perm_i, cdouble rel_perm_t, double lamda, double width, bool inf_wall, cdouble& ref_coeff, cdouble& tran_coeff)
{
}

void FieldCompute::GetTMCoeff(cdouble theta_i, cdouble theta_t, cdouble rel_perm_i, cdouble rel_perm_t, double lamda, double width, bool inf_wall, cdouble& ref_coeff, cdouble& tran_coeff)
{
}

cdouble FieldCompute::GetTransAngle(double thetaIncident, cdouble efieldIncident, cdouble efieldTransmit)
{
	return cdouble();
}
