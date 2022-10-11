#include "FieldCompute.h"

FieldCompute::FieldCompute()
{
}

FieldCompute::~FieldCompute()
{
}



Vec3c FieldCompute::FieldAtReceiver(int receiverId)
{
	// Get all paths for the receiver
	for (int i = 0; i < _pathsCount; i++)
	{
		//_rayPaths[i][receiverId].rayPaths
	}
	return Vec3c();
}

Vec3c FieldCompute::FieldForPath(const std::vector<Ray>& path, const Mat3& txCoordSys)
{
	/* This is for testing purposes only */
	double lamda = SPEED_OF_LIGHT / 1e9;
	double k = (2 * PI) / lamda;

	Mat3 currentCoordSys = Mat3::Identity(); // Global coordinate system
	Mat3 nextCoordSys = txCoordSys;


	// Step 1: convert ray from global coordinate system to local coordinate system
	// Step 2: once in local coordinates, convert to spherical coordinates
	// Step 3: find field contribution in spherical coordinate system
	// Step 4:  

	Vec3c totalField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };
	cdouble propagation_term{ 0,0 };

	for (int i = 0; i < int(path.size()); i++)
	{
		const Ray& ray = path[i];
		Vec3 pGlobal = ray.targetPoint - ray.sourcePoint; // This is just a translation to global origin
		Vec3 pLocal = PointInNewCoordSys(pGlobal, currentCoordSys, nextCoordSys);
		Vec3 pSph = CartesianToSpherical(pLocal);

		double r = pSph(0);
		double theta = pSph(1);
		double phi = pSph(2);

		propagation_term.real(cos(-1 * k * r) / r);
		propagation_term.imag(sin(-1 * k * r) / r);

		if (i == 0)
		{
			// The first ray is always from TX to the next facet or receiver
			// Get either the gain from the TX or an analytical pattern
			// E(r,theta,phi) = E(0) * exp(-jkr)/r
			totalField += propagation_term * GetAnalyticEfieldPattern(1, theta, phi, 1);
		}


		
	}

	/*
	for (const Ray& r : path)
	{


		
		cdouble propagation_term;
		propagation_term.real(cos(-1 * k * pSph(0) / pSph(0)));
		propagation_term.imag(sin(-1 * k * pSph(0) / pSph(0)));

		Vec3 surfaceNormal = _triangleMesh[r.hitSurfaceID]->norm;
		MaterialProperties m = _materials[_triangleMesh[r.hitSurfaceID]->materialId];

		
	}
	*/
	return Vec3c();
}

Vec3c FieldCompute::GetAnalyticEfieldPattern(int antenna_type, double theta, double phi, double pt)
{
	Vec3c e_field_sph;

	if (antenna_type == 1)
	{
		// Field of Hertzian Dipole in TX Spherical Coordinate System
		e_field_sph(0) = 0;
		e_field_sph(1) = sin(theta);
		e_field_sph(2) = 0;
	}
	else if (antenna_type == 2)
	{
		// Field of half wave dipole in TX Spherical Coordinate System
		e_field_sph(0) = 0;
		e_field_sph(1) = sqrt(60 * pt) * (cos(PI * cos(theta) / 2.0) / sin(theta));
		e_field_sph(2) = 0;
	}
	else if (antenna_type == 3)
	{
		// Field of antenna array in TX Spherical Coordinate System
	}
	return e_field_sph;
}

Vec3c FieldCompute::ComputeRefcField(const Vec3c& efield_i_sph, double rel_perm, double sigma, double freq, double thetaIncident, double width, bool inf_wall)
{
	cdouble epsilonTransmit;
	epsilonTransmit.real(rel_perm);
	epsilonTransmit.imag(sigma / (E0 * 2 * PI * freq));

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
	epsilonTransmit.imag(sigma / (E0 * 2 * PI * freq));

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

	cdouble gamma_te = ((n_i * cos_i) - (n_t * cos_t)) - ((n_i * cos_i) + (n_t * cos_t));

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
	cdouble t_i(thetaIncident, 0);
	cdouble t_t;

	cdouble n_i = sqrt(epsilonIncident);
	cdouble n_t = sqrt(epsilonTransmit);
	t_t = asin(n_i * sin(t_i) / n_t);
	return t_t;
}

Mat3 FieldCompute::GetSurfCoordSys(const Vec3& n, const Ray& rayIncident)
{
	// Need to fix: if the ray is normal incident, then the code return nans
	
	// Fail-safe procedure: ensure both ray direction and wall norm is unit vector
	Vec3 normal = n.normalized();
	Vec3 rayDir = (rayIncident.targetPoint - rayIncident.sourcePoint).normalized();

	// This would fail if the normal and rayDir are not unit vectors
	double theta_i = acos(rayDir.dot(normal));

	//double theta_i = acos((zw*inc_ray.dir)/(mag(zw)*mag(inc_ray.dir)));
	if (theta_i > PI / 2)
	{
		theta_i = PI - theta_i;

		// From Neeraj's ray tracer code:
		// If the angle between incident ray and the surface normal is > 90 degrees,
		//	this indicates that the surface normal is inward pointing. 
		// 	Reverse the surface normal and also reverse orientation of vertices

		//wall.unit_norm_ = -1*wall.unit_norm_;
		//wall.d = - wall.d;
	}

	Vec3 zw = normal;
	Vec3 yw = rayDir.cross(normal) / sin(theta_i);
	Vec3 xw = yw.cross(zw);

	Mat3 surfaceCoord;
	surfaceCoord(0, Eigen::all) = xw;
	surfaceCoord(1, Eigen::all) = yw;
	surfaceCoord(2, Eigen::all) = zw;

	return surfaceCoord;
}

