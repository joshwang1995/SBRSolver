#include "FieldCompute.h"

FieldCompute::FieldCompute()
{
}

FieldCompute::~FieldCompute()
{
}

Vec3c FieldCompute::FieldAtReceiver(int receiverId)
{
#if DEBUG_LEVEL > 1
	std::cout << "\n[Entering] FieldAtReceiver..." << std::endl;
#endif
	Vec3c totalField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };

	// Get all paths for the receiver
	for (int i = 0; i < _pathsCount; i++)
	{
		int numRayPaths = int(_rayPaths[i][receiverId].rayPaths.size());
		for (int j = 0; j < numRayPaths; j++)
		{
#if DEBUG_LEVEL > 1
			std::cout << "\nStart Field Computation for path " << i << " at receiver " << receiverId << std::endl;
#endif
			totalField += FieldForPath(_rayPaths[i][receiverId].rayPaths[j]);
		}
	}
#if DEBUG_LEVEL > 1
	std::cout << std::endl;
	std::cout << "\t\tCombined Field Components (Ex, Ey, Ez) = " << totalField.x() << ", " << totalField.y() << ", " << totalField.z() << std::endl;
	std::cout << "\t\tCombined Field Components (|Ex|,|Ey|,|Ez|) = " << abs(totalField.x()) << ", " << abs(totalField.y()) << ", " << abs(totalField.z()) << std::endl;
	std::cout << "\t\tCombined Field Components arg(Ex, Ey, Ez) = " << Rad2Deg(std::arg(totalField.x())) << ", " << Rad2Deg(std::arg(totalField.y())) << ", " << Rad2Deg(std::arg(totalField.z())) << std::endl;
	std::cout << "\t\tCombined Field Strength (dBuVm): " << 20.0 * log10(totalField.norm() * 1e6 / sqrt(2)) << std::endl;
	std::cout << "[Leaving] FieldAtReceiver...\n" << std::endl;
#endif
	return totalField;
}

void FieldCompute::RefCoeffTest(int numPts, double freq, cdouble epsilon1, cdouble epsilon2, std::string fname)
{
	double lamda = SPEED_OF_LIGHT / freq;
	Vec thetaArray = Vec::LinSpaced(numPts, 0.0, PI / 2);

	std::ofstream ofs;
	ofs.open(fname);

	cdouble refTE, refTM, transTE, transTM, theta_t;
	ofs << "theta_i,refTM_mag,refTM_phase,refTE_mag, refTE_phase" << std::endl;
	for (auto theta_i : thetaArray)
	{
		theta_t = GetTransAngle(theta_i, epsilon1, epsilon2);
		GetTMCoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTM, transTM);
		GetTECoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTE, transTE);
		ofs << Rad2Deg(theta_i) << "," << std::abs(refTM) << "," << std::arg(refTM) << ",";
		ofs << std::abs(refTE) << "," << std::arg(refTE) << std::endl;
	}

	ofs.close();
}

void FieldCompute::TransCoeffTest(int numPts, double freq, cdouble epsilon1, cdouble epsilon2, std::string fname)
{
	double lamda = SPEED_OF_LIGHT / freq;
	Vec thetaArray = Vec::LinSpaced(numPts, 0.0, PI / 2);

	std::ofstream ofs;
	ofs.open(fname);

	cdouble refTE, refTM, transTE, transTM, theta_t;
	ofs << "theta_i,tranTM_mag,tranTM_phase, tranTE_mag, tranTE_phase" << std::endl;
	for (auto theta_i : thetaArray)
	{
		theta_t = GetTransAngle(theta_i, epsilon1, epsilon2);
		GetTMCoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTM, transTM);
		GetTECoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTE, transTE);
		ofs << Rad2Deg(theta_i) << "," << std::abs(transTM) << "," << std::arg(transTM) << ",";
		ofs << std::abs(transTE) << "," << std::arg(transTE) << std::endl;
	}

	ofs.close();
}

Vec3c FieldCompute::FieldForPath(const std::vector<Ray>& path)
{
#if DEBUG_LEVEL > 1
	using namespace std;
	Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", " ; ", "", "");
	cout << "[Entering] FieldForPath..." << endl;
#endif

	Vec3c totalField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };
	double totalPathLength = 0.0;
	int numRef = 0;

	for (int i = 0; i < int(path.size()); i++)
	{
		//cout << "\tComputing field for Ray" << i << "..." << endl;

		// Extract information from current ray
		const Ray& ray = path[i];
		Vec3 vecGlobal = ray.targetPoint - ray.sourcePoint;
		totalPathLength += vecGlobal.norm();

		//cout << "\t\tSource Point: " << ray.sourcePoint.transpose().format(CleanFmt) << endl;
		//cout << "\t\tTarget Point: " << ray.targetPoint.transpose().format(CleanFmt) << endl;
		//cout << "\t\tNext Hit Surface ID: " << ray.hitSurfaceID << endl;
		//if (ray.hitSurfaceID > -1)
		//{
			//cout << "\t\tNext Hit Surface Normal: " << _triangleMesh->at(ray.hitSurfaceID)->norm.transpose().format(CleanFmt) << endl;
		//}
		//cout << "\t\tRay vector global cartesian: " << vecGlobal.transpose().format(CleanFmt) << endl;

		/* Updating the incident field*/
		if (i == 0)
		{
			Vec3 vecSph = CartesianToSpherical(RotateToNewCoordSys(vecGlobal, globalCoordSys, txCoordSys));
			totalField = GetAnalyticEfieldPattern(0, vecSph[1], vecSph[2], txPower); // local spherical
			//cout << "\t\tTX Coordinate System: " << txCoordSys.format(CleanFmt) << endl;
			//cout << "\t\t\tTX theta: " << Rad2Deg(vecSph[1]) << endl;
			//cout << "\t\t\tTX phi: " << Rad2Deg(vecSph[2]) << endl;
			//cout << "\t\tIncident Field Spherical: " << incidentField.transpose().format(CleanFmt) << endl;
			totalField = RotateToNewCoordSys(SphericalToCartesianVector(totalField, vecSph[1], vecSph[2]), txCoordSys, globalCoordSys);
		}

		/* Updating the reflected or transmitted field*/
		if (ray.reflectionMaterialId >= 0)
		{
			numRef++;
			Vec3 kIncident = path[i-1].targetPoint - path[i-1].sourcePoint;
			totalField = ComputeRefcField(kIncident, vecGlobal, totalField, ray.reflectionMaterialId, GetSurfCoordSys(path[i-1].hitSurfaceID, ray));
		}
		else if (ray.penetrationMaterialId >= 0)
		{
			Vec3 kIncident = path[i - 1].targetPoint - path[i - 1].sourcePoint;
			totalField = ComputeTransField(kIncident, vecGlobal, totalField, ray.penetrationMaterialId, GetSurfCoordSys(path[i - 1].hitSurfaceID, ray)); // local spherical
		}
		//cout << "\t\tTotal Field: " << totalField.transpose().format(CleanFmt) << std::endl;
	}
	totalField = totalField * exp(-j * k * totalPathLength) / totalPathLength;

#if DEBUG_LEVEL > 1
	//cout << endl;
	cout << "\tDelay [ns]: " << (totalPathLength/ SPEED_OF_LIGHT)*1e9 << endl;
	cout << "\tPhase [deg]: " << Rad2Deg(WrapToTwoPi(k * totalPathLength + numRef * PI)) << endl;
	cout << "\tField Strength [dBuvm]: " << 20.0 * log10(totalField.norm() * 1e6 / sqrt(2)) << endl;
	cout << "\tField Components [Ex, Ey, Ez]: " << totalField.transpose().format(CleanFmt) << endl;
	cout << "[Leaving] FieldForPath..." << endl;
#endif

	return totalField;
}

Vec3c FieldCompute::GetAnalyticEfieldPattern(int antennaType, double theta, double phi, double pt)
{
	Vec3c eFieldSph;
	if (antennaType == 0)
	{
		// Isotropic antenna
		eFieldSph(0) = 0;
		eFieldSph(1) = cdouble(0, sqrt(60 * pt));
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

Vec3c FieldCompute::ComputeRefcField(const Vec3& kIncident, const Vec3& kReflect, const Vec3c& efieldGlobal, int materialId, const Mat3& surfCoordSys)
{
	// Permittivity of air
	cdouble epsilon_i(1, 0);

	// Material property of the impinging medium
	cdouble epsilon_t;
	epsilon_t.real(_materials[materialId].relPermittivityRe);
	epsilon_t.imag(-1 * _materials[materialId].relConductivity/ (E0 * 2 * PI * frequency));
	
	// Get theta_i and theta_t
	double theta_i = GetIncidentAngle(kIncident, surfCoordSys(2, Eigen::indexing::all));
	cdouble theta_t = GetTransAngle(theta_i, epsilon_i, epsilon_t);

	cdouble refTE, refTM, transTE, transTM;
	GetTECoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTE, transTE);
	GetTMCoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTM, transTM);

	// Rotate field from global coordinate system to current coordinate system
	Vec3c efield_local_cart = RotateToNewCoordSys(efieldGlobal, globalCoordSys, surfCoordSys);
	Vec3c efield_local_sph = CartesianToSphericalVector(efield_local_cart, theta_i, 0);
	Vec3c refcField (efield_local_sph[0], efield_local_sph[1] * refTM, efield_local_sph[2] * refTE);
	Vec3c refcField_cart = RotateToNewCoordSys(SphericalToCartesianVector(refcField, theta_i, PI), surfCoordSys, globalCoordSys);

#if DEBUG_LEVEL > 2
	using namespace std;
	Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", " ; ", "", "");
	cout << "\t\tSurface Coordinate System: " << surfCoordSys.format(CleanFmt) << endl;
	cout << "\t\tIncident wave vector: " << kIncident.transpose().format(CleanFmt) << endl;
	cout << "\t\t\tAngle between k_i and normal: " << Rad2Deg(theta_i) << endl;
	cout << "\t\t\t|RefTM| = " << abs(refTM) << ", arg(RefTM) = " << Rad2Deg(std::arg(refTM)) << endl;
	cout << "\t\t\t|RefTE| = " << abs(refTE) << ", arg(RefTE) = " << Rad2Deg(std::arg(refTE)) << endl;
	cout << "\t\tEfield local cartesian: " << efield_local_cart.transpose().format(CleanFmt) << endl;
	cout << "\t\tEfield local spherical: " << efield_local_sph.transpose().format(CleanFmt) << endl;
	cout << "\t\tReflected field local spherical: " << refcField.transpose().format(CleanFmt) << endl;
	cout << "\t\tReflected field global: " << refcField_cart.transpose().format(CleanFmt) << endl;
#endif

	return refcField_cart;
}

Vec3c FieldCompute::ComputeTransField(const Vec3& kIncident, const Vec3& kReflect, const Vec3c& efieldGlobal, int materialId, const Mat3& surfCoordSys)
{
	// Permittivity of air
	cdouble epsilon_i(1, 0);

	// Material property of the impinging medium
	cdouble epsilon_t;
	epsilon_t.real(_materials[materialId].relPermittivityRe);
	epsilon_t.imag(-1 * _materials[materialId].relConductivity / (E0 * 2 * PI * frequency));

	// Get theta_i and theta_t
	double theta_i = GetIncidentAngle(kIncident, surfCoordSys(2, Eigen::indexing::all));
	cdouble theta_t = GetTransAngle(theta_i, epsilon_i, epsilon_t);

	cdouble refTE, refTM, transTE, transTM;
	GetTECoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTE, transTE);
	GetTMCoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTM, transTM);

	// Rotate field from global coordinate system to current coordinate system
	Vec3 k_i_local_sph = CartesianToSpherical(RotateToNewCoordSys(kIncident, globalCoordSys, surfCoordSys));
	Vec3c efield_local_cart = RotateToNewCoordSys(efieldGlobal, globalCoordSys, surfCoordSys);
	Vec3c efield_local_sph = CartesianToSphericalVector(efield_local_cart, k_i_local_sph[1], k_i_local_sph[2]);
	Vec3c transField(efield_local_sph[0], efield_local_sph[1] * transTM, efield_local_sph[2] * transTE);
	Vec3c transField_cart = RotateToNewCoordSys(SphericalToCartesianVector(transField, k_i_local_sph[1], k_i_local_sph[2] * 0), surfCoordSys, globalCoordSys);

#if DEBUG_LEVEL > 2
	using namespace std;
	Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", " ; ", "", "");
	cout << "\t\tSurface Coordinate System: " << surfCoordSys.format(CleanFmt) << endl;
	cout << "\t\tIncident wave vector: " << kIncident.transpose().format(CleanFmt) << endl;
	cout << "\t\tIncident wave vector Spherical: " << k_i_local_sph[0] << ", " << Rad2Deg(k_i_local_sph[1]) << ", " << Rad2Deg(k_i_local_sph[2]) << endl;
	cout << "\t\t\tAngle between k_i and normal: " << Rad2Deg(theta_i) << endl;
	cout << "\t\t\tTransTM: " << transTM << endl;
	cout << "\t\t\tTransTE: " << transTE << endl;
	cout << "\t\tEfield local cartesian: " << efield_local_cart.transpose().format(CleanFmt) << endl;
	cout << "\t\tEfield local spherical: " << efield_local_sph.transpose().format(CleanFmt) << endl;
	cout << "\t\tTransmitted field local spherical: " << transField.transpose().format(CleanFmt) << endl;
	cout << "\t\tTransmitted field global: " << transField_cart.transpose().format(CleanFmt) << endl;
#endif

	return transField_cart;
}

void FieldCompute::GetTECoeff(cdouble theta_i, cdouble theta_t, cdouble rel_perm_i, cdouble rel_perm_t, double width, bool inf_wall, cdouble& ref_coeff, cdouble& tran_coeff)
{
	cdouble n_i, n_t, eta_i, eta_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);
	eta_i = ETA0 / n_i;
	eta_t = ETA0 / n_t;
	cdouble cos_i = cos(theta_i);
	cdouble cos_t = cos(theta_t);
	cdouble sin_i = sin(theta_i);

	cdouble gamma_te = (-eta_t * cos_i + eta_i * cos_t) / (eta_t * cos_i + eta_i * cos_t);
	cdouble tau_te = (2.0 * eta_t * cos_i) / (eta_t * cos_i + eta_i * cos_t);
	if (inf_wall)
	{
		ref_coeff = gamma_te;
		tran_coeff = tau_te;
	}
	else
	{
		cdouble q = (TWOPI * width / lamda) * sqrt(n_t - (sin_i * sin_i));
		cdouble phi_factor = exp(-j * 2.0 * q);
		cdouble denom = 1.0 - (gamma_te * gamma_te * phi_factor);
		ref_coeff = (gamma_te * (1.0 - phi_factor)) / denom;
		tran_coeff = ((1.0 - gamma_te * gamma_te) * phi_factor) / denom;
	}
}

void FieldCompute::GetTMCoeff(cdouble theta_i, cdouble theta_t, cdouble rel_perm_i, cdouble rel_perm_t, double width, bool inf_wall, cdouble& ref_coeff, cdouble& tran_coeff)
{
	cdouble n_i, n_t, eta_i, eta_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);
	eta_i = ETA0 / n_i;
	eta_t = ETA0 / n_t;
	cdouble cos_i = cos(theta_i);
	cdouble cos_t = cos(theta_t);
	cdouble sin_i = sin(theta_i);

	cdouble gamma_tm = (-eta_t * cos_t + eta_i * cos_i) / (eta_t * cos_t + eta_i * cos_i);
	cdouble tau_tm = (2.0 * eta_t * cos_i) / (eta_t * cos_t + eta_i * cos_i);
	if (inf_wall)
	{
		ref_coeff = gamma_tm;
		tran_coeff = tau_tm;
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

inline double FieldCompute::GetIncidentAngle(const Vec3& v, const Vec3& normal)
{
	double angleBetween = AngleBetween(v, normal, false, true);
	return angleBetween >= PI / 2.0 ? PI - angleBetween : angleBetween;
}

inline cdouble FieldCompute::GetTransAngle(double thetaIncident, cdouble epsilon_i, cdouble epsilon_t)
{
	// Snell's Law
	cdouble theta_i(thetaIncident, 0);
	cdouble n_i = sqrt(epsilon_i);
	cdouble n_t = sqrt(epsilon_t);
	cdouble theta_t = asin((n_i/n_t) * sin(theta_i));
	return theta_t;
}

Mat3 FieldCompute::GetSurfCoordSys(const int& hitSurfId, const Ray& rayIncident)
{
	// Handle the case where the next hit Surfac is a receiver
	if (hitSurfId == -1)
	{
		return globalCoordSys;
	}

	// Fail-safe procedure: ensure both ray direction and wall norm is unit vector
	Vec3 normal = _triangleMesh->at(hitSurfId)->norm;
	Vec3 rayDir = (rayIncident.targetPoint - rayIncident.sourcePoint).normalized();

	// This would fail if the normal and rayDir are not unit vectors
	double theta_i = acos(rayDir.dot(normal));

	//double theta_i = acos((zw*inc_ray.dir)/(mag(zw)*mag(inc_ray.dir)));
	if (theta_i > PI / 2)
	{
		theta_i = PI - theta_i;
	}

	if (theta_i == 0.0)
	{
		return globalCoordSys;
	}

	Vec3 zw = normal;
	Vec3 yw = rayDir.cross(normal) / sin(theta_i);
	Vec3 xw = yw.cross(zw);

	Mat3 surfaceCoord;
	surfaceCoord(0, Eigen::indexing::all) = xw;
	surfaceCoord(1, Eigen::indexing::all) = yw;
	surfaceCoord(2, Eigen::indexing::all) = zw;

	return surfaceCoord;
}