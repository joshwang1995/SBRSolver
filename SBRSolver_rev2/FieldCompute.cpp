#include "FieldCompute.h"

FieldCompute::FieldCompute()
{
}

FieldCompute::~FieldCompute()
{
}

Vec3c FieldCompute::FieldAtReceiver(int receiverId)
{
	Vec3c totalField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };

	// Get all paths for the receiver
	for (int i = 0; i < _pathsCount; i++)
	{
		int numRayPaths = int(_rayPaths[i][receiverId].rayPaths.size());
		for (int j = 0; j < numRayPaths; j++)
		{
			totalField += FieldForPath(_rayPaths[i][receiverId].rayPaths[j]);
		}
	}
	return totalField;
}

Vec3c FieldCompute::FieldForPath(const std::vector<Ray>& path)
{
	double lamda = SPEED_OF_LIGHT / _frequency;
	double k = (2 * PI) / lamda;

	Mat3 globalCoordSys = Mat3::Identity();
	Mat3 currentCoordSys = globalCoordSys; // Global coordinate system
	Mat3 nextCoordSys = _txCoordSys;

	Vec3c incidentField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };
	Vec3c totalField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };
	cdouble propagation_term{ 0,0 };

	for (int i = 0; i < int(path.size()); i++)
	{
		const Ray& ray = path[i];
		Vec3 vecGlobal = ray.targetPoint - ray.sourcePoint;
		Vec3 vecLocal = RotateToNewCoordSys(vecGlobal, globalCoordSys, nextCoordSys);
		Vec3 vecSph = CartesianToSpherical(vecLocal);

		currentCoordSys = nextCoordSys;
		if (ray.hitSurfaceID == -1)
		{
			nextCoordSys = globalCoordSys;
		}
		else
		{
			nextCoordSys = _triangleMesh->at(ray.hitSurfaceID)->coordSys;
		}
		

		double r = vecSph(0);
		double theta = vecSph(1);
		double phi = vecSph(2);
		propagation_term.real(cos(-1 * k * r) / r);
		propagation_term.imag(sin(-1 * k * r) / r);

		if (i == 0)
		{
			// The first ray is always from TX to the next facet or receiver
			// Get either the gain from the TX or an analytical pattern
			// E(r,theta,phi) = E(0) * exp(-jkr)/r
			// E(0) need to be calculated from the field coefficient sqrt(eta * Pr * Gain / 2Pi)
			incidentField = propagation_term * GetAnalyticEfieldPattern(0, theta, phi, _txPower); // local spherical
		}
		else
		{
			incidentField = RotateToNewCoordSys(totalField, globalCoordSys, currentCoordSys);
			incidentField = CartesianToSphericalVector(incidentField, theta, phi);
			incidentField = propagation_term * incidentField;
		}

		// Flip the theta phi components accroding to Catedra
		incidentField(1) = -incidentField(1);
		incidentField(2) = -incidentField(2);

		if (ray.reflectionMaterialId >= 0)
		{
			double incidentAngle = AngleBetween(vecGlobal, currentCoordSys(2,Eigen::all), false, true);
			double relPerm = _materials[ray.reflectionMaterialId].relPermittivityRe;
			double sigma = _materials[ray.reflectionMaterialId].relConductivity;
			totalField = ComputeRefcField(incidentField, relPerm, sigma, _frequency, incidentAngle, 0, _useFresnelCoeff); // local spherical
		}
		else if (ray.penetrationMaterialId >= 0)
		{
			double incidentAngle = AngleBetween(vecGlobal, currentCoordSys(2, Eigen::all), false, true);
			double relPerm = _materials[ray.penetrationMaterialId].relPermittivityRe;
			double sigma = _materials[ray.penetrationMaterialId].relConductivity;
			totalField = ComputeTransField(incidentField, relPerm, sigma, _frequency, incidentAngle, 0, _useFresnelCoeff);
		}
		else
		{
			totalField = incidentField;
		}
		// Need to do: add a case where receiver gain can be applied here
		totalField = SphericalToCartesianVector(totalField, theta, phi);
		totalField = RotateToNewCoordSys(totalField, currentCoordSys, globalCoordSys);
	}

	return totalField;
}

Vec3c FieldCompute::GetAnalyticEfieldPattern(int antennaType, double theta, double phi, double pt)
{
	Vec3c eFieldSph;
	if (antennaType == 0)
	{
		// Isotropic antenna
		eFieldSph(0) = cdouble(sqrt(60 * pt), 0);
		eFieldSph(1) = 0;
		eFieldSph(2) = 0;
	}

	if (antennaType == 1)
	{
		// Field of Hertzian Dipole in TX Spherical Coordinate System
		eFieldSph(0) = 0;
		eFieldSph(1) = cdouble(0, sqrt(60*pt*1.5) * sin(theta));
		eFieldSph(2) = 0;
	}
	else if (antennaType == 2)
	{
		// Field of half wave dipole in TX Spherical Coordinate System
		eFieldSph(0) = 0;
		eFieldSph(1) = cdouble(0, sqrt(60 * pt) * (cos(PI * cos(theta) / 2.0) / sin(theta)));
		eFieldSph(2) = 0;
	}
	else if (antennaType == 3)
	{
		// Field of antenna array in TX Spherical Coordinate System
	}
	return eFieldSph;
}

Vec3c FieldCompute::ComputeRefcField(const Vec3c& efield_i_sph, double rel_perm, double sigma, double freq, double thetaIncident, double width, bool inf_wall)
{
	cdouble epsilonTransmit;
	epsilonTransmit.real(rel_perm);
	epsilonTransmit.imag(-1 * sigma / (E0 * 2 * PI * freq));

	double lamda = SPEED_OF_LIGHT / freq;

	// Permittivity of air
	cdouble epsilonIncident(1, 0);

	cdouble theta_t = GetTransAngle(thetaIncident, epsilonIncident, epsilonTransmit);

	cdouble refTE, refTM, transTE, transTM;
	GetTMCoeff(thetaIncident, theta_t, epsilonIncident, epsilonTransmit, lamda, width, inf_wall, refTM, transTM);
	GetTECoeff(thetaIncident, theta_t, epsilonIncident, epsilonTransmit, lamda, width, inf_wall, refTE, transTE);

	// x is r, y is theta, z is phi
	cdouble Er = efield_i_sph[0];
	cdouble Etheta = efield_i_sph[1];
	cdouble Ephi = efield_i_sph[2];

	Vec3c refcField (Er, Etheta * refTM, Ephi * refTE); 
	return refcField;
}

Vec3c FieldCompute::ComputeTransField(const Vec3c& efield_i_sph, double rel_perm, double sigma, double freq, double thetaIncident, double width, bool inf_wall)
{
	cdouble epsilonTransmit;
	epsilonTransmit.real(rel_perm);
	epsilonTransmit.imag(-1 * sigma / (E0 * 2 * PI * freq));

	double lamda = SPEED_OF_LIGHT / freq;

	// Permittivity of air
	cdouble epsilonIncident(1, 0);

	cdouble theta_t = GetTransAngle(thetaIncident, epsilonIncident, epsilonTransmit);

	cdouble refTE, refTM, transTE, transTM;
	GetTMCoeff(thetaIncident, theta_t, epsilonIncident, epsilonTransmit, lamda, width, inf_wall, refTM, transTM);
	GetTECoeff(thetaIncident, theta_t, epsilonIncident, epsilonTransmit, lamda, width, inf_wall, refTE, transTE);

	// x is r, y is theta, z is phi
	cdouble Er = efield_i_sph[0];
	cdouble Etheta = efield_i_sph[1];
	cdouble Ephi = efield_i_sph[2];

	Vec3c transField(Er, Etheta * transTM, Ephi * transTE);
	return transField;
}

void FieldCompute::GetTECoeff(cdouble theta_i, cdouble theta_t, cdouble rel_perm_i, cdouble rel_perm_t, double lamda, double width, bool inf_wall, cdouble& ref_coeff, cdouble& tran_coeff)
{
	cdouble n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);
	cdouble cos_i = cos(theta_i);
	cdouble cos_t = cos(theta_t);
	cdouble sin_i = sin(theta_i);

	cdouble gamma_te = ((n_i * cos_i) - (n_t * cos_t)) / ((n_i * cos_i) + (n_t * cos_t));
	if (inf_wall)
	{
		ref_coeff = gamma_te;
		tran_coeff = 0;
	}
	else
	{
		cdouble q = (TWOPI * width / lamda) * sqrt(n_t - (sin_i * sin_i));
		cdouble j(0.0, 1.0);
		cdouble phi_factor = exp(-j * 2.0 * q);
		cdouble denom = 1.0 - (gamma_te * gamma_te * phi_factor);
		ref_coeff = (gamma_te * (1.0 - phi_factor)) / denom;
		tran_coeff = ((1.0 - gamma_te * gamma_te) * phi_factor) / denom;
	}
}

void FieldCompute::GetTMCoeff(cdouble theta_i, cdouble theta_t, cdouble rel_perm_i, cdouble rel_perm_t, double lamda, double width, bool inf_wall, cdouble& ref_coeff, cdouble& tran_coeff)
{
	cdouble n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);
	cdouble cos_i = cos(theta_i);
	cdouble cos_t = cos(theta_t);
	cdouble sin_i = sin(theta_i);

	cdouble gamma_tm = (n_t * cos_i - n_i * cos_t) / (n_i * cos_t + n_t * cos_i);

	if (inf_wall)
	{
		ref_coeff = gamma_tm;
		tran_coeff = 0;
	}
	else
	{
		cdouble q = (TWOPI * width / lamda) * sqrt(n_t - (sin_i * sin_i));
		cdouble phi_factor = exp(-j * 2.0 * q);
		cdouble denom = 1.0 - (gamma_tm * gamma_tm * phi_factor);
		ref_coeff = (gamma_tm * (1.0 - phi_factor)) / denom;
		tran_coeff = ((1.0 - gamma_tm * gamma_tm) * phi_factor) / denom;
	}
}

cdouble FieldCompute::GetTransAngle(double thetaIncident, cdouble epsilonIncident, cdouble epsilonTransmit)
{
	// Snell's Law
	cdouble theta_i(thetaIncident, 0);
	

	cdouble n_i = sqrt(epsilonIncident);
	cdouble n_t = sqrt(epsilonTransmit);
	cdouble theta_t = asin(n_i * sin(theta_i) / n_t);
	return theta_t;
}

