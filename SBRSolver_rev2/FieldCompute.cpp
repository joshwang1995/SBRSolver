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

Vec3c FieldCompute::FieldForPath(const std::vector<Ray>& path)
{
	for (const Ray& r : path)
	{
		Vec3 rayDir = r.targetPoint - r.sourcePoint;
		Vec3 rayOrig = r.sourcePoint;
		Vec3 surfaceNormal = _triangleMesh[r.hitSurfaceID]->norm;
		MaterialProperties m = _materials[_triangleMesh[r.hitSurfaceID]->materialId];
	}
	return Vec3c();
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