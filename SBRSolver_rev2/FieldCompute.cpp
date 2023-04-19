#include "FieldCompute.h"

FieldCompute::FieldCompute()
{
}

FieldCompute::~FieldCompute()
{
}

Vec3c FieldCompute::FieldAtReceiver(int receiverId)
{
	std::cout << "\n[Entering] FieldAtReceiver..." << std::endl;

	Vec3c totalField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };

	// Get all paths for the receiver
	for (int i = 0; i < _pathsCount; i++)
	{
		int numRayPaths = int(_rayPaths[i][receiverId].rayPaths.size());
		for (int j = 0; j < numRayPaths; j++)
		{
			std::cout << "\nStart Field Computation for path " << i << " at receiver " << receiverId << std::endl;
			totalField += FieldForPath(_rayPaths[i][receiverId].rayPaths[j]);
		}
	}

	std::cout << "[Leaving] FieldAtReceiver...\n" << std::endl;
	return totalField;
}

void FieldCompute::RefCoeffTest(int numPts, double freq, cdouble epsilon1, cdouble epsilon2, std::string fname)
{
	double lamda = SPEED_OF_LIGHT / freq;
	Vec thetaArray = Vec::LinSpaced(numPts, 0.0, EIGEN_PI / 2);

	std::ofstream ofs;
	ofs.open(fname);

	cdouble refTE, refTM, transTE, transTM, theta_t;
	ofs << "theta_i,refTM_mag,refTM_phase,refTE_mag, refTE_phase" << std::endl;
	for (auto theta_i : thetaArray)
	{
		theta_t = GetTransAngle(theta_i, epsilon1, epsilon2);
		GetTMCoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTM, transTM);
		GetTECoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTE, transTE);
		ofs << theta_i * (180.0/EIGEN_PI) << "," << std::abs(refTM) << "," << std::arg(refTM) << ",";
		ofs << std::abs(refTE) << "," << std::arg(refTE) << std::endl;
	}

	ofs.close();
}

void FieldCompute::TransCoeffTest(int numPts, double freq, cdouble epsilon1, cdouble epsilon2, std::string fname)
{
	double lamda = SPEED_OF_LIGHT / freq;
	Vec thetaArray = Vec::LinSpaced(numPts, 0.0, EIGEN_PI / 2);

	std::ofstream ofs;
	ofs.open(fname);

	cdouble refTE, refTM, transTE, transTM, theta_t;
	ofs << "theta_i,tranTM_mag,tranTM_phase, tranTE_mag, tranTE_phase" << std::endl;
	for (auto theta_i : thetaArray)
	{
		theta_t = GetTransAngle(theta_i, epsilon1, epsilon2);
		GetTMCoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTM, transTM);
		GetTECoeff(theta_i, theta_t, epsilon1, epsilon2, 0.0, true, refTE, transTE);
		ofs << theta_i * (180.0 / EIGEN_PI) << "," << std::abs(transTM) << "," << std::arg(transTM) << ",";
		ofs << std::abs(transTE) << "," << std::arg(transTE) << std::endl;
	}

	ofs.close();
}

Vec3c FieldCompute::FieldForPath(const std::vector<Ray>& path)
{
	using namespace std;
	Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", " ; ", "", "");
	cout << "[Entering] FieldForPath..." << endl;

	Mat3 currentCoordSys = globalCoordSys;
	Mat3 nextCoordSys = txCoordSys;

	Vec3c incidentField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };
	Vec3c totalField{ cdouble(0,0), cdouble(0,0), cdouble(0,0) };
	double totalPathLength = 0.0;
	
	for (int i = 0; i < int(path.size()); i++)
	{
		cout << "\tComputing field for Ray" << i << "..." << endl;

		// Extract information from current ray
		const Ray& ray = path[i];
		Vec3 vecGlobal = ray.targetPoint - ray.sourcePoint;
		Vec3 vecLocal = RotateToNewCoordSys(vecGlobal, globalCoordSys, nextCoordSys);
		Vec3 vecSph = CartesianToSpherical(vecLocal);
		
		double r = vecSph(0);
		double theta = vecSph(1);
		double phi = vecSph(2);
		totalPathLength += r;

		cout << "\t\tRay vector global cartesian: " << vecGlobal.transpose().format(CleanFmt) << endl;
		cout << "\t\tRotate ray from global to this coordinate system: " << nextCoordSys.format(CleanFmt) << endl;
		cout << "\t\tRay vector rotated to local: " << vecLocal.transpose().format(CleanFmt) << endl;
		cout << "\t\tRay vector in local spherical coodinates: " << r << ", " << theta * (180.0 / PI) << ", " << phi * (180.0 / PI) << endl;

		/* Updating the coordinate systems*/
		currentCoordSys = nextCoordSys;
		nextCoordSys = GetSurfCoordSys(ray.hitSurfaceID, ray);
		cout << "\t\tUpdated Current Coordinate System: " << currentCoordSys.format(CleanFmt) << endl;
		cout << "\t\tUpdated Next Coordinate System: " << nextCoordSys.format(CleanFmt) << endl;

		/* Updating the incident field*/
		if (i == 0)
		{
			// The first ray is always from TX to the next facet or receiver
			// Get either the gain from the TX or an analytical pattern
			// E(r,theta,phi) = E(0) * exp(-jkr)/r
			// E(0) need to be calculated from the field coefficient sqrt(eta * Pr * Gain / 2Pi)
			incidentField = GetAnalyticEfieldPattern(0, theta, phi, txPower); // local spherical
		}
		else
		{
			incidentField = RotateToNewCoordSys(totalField, globalCoordSys,currentCoordSys);
			incidentField = CartesianToSphericalVector(incidentField, theta, phi);
		}

		/* Updating the reflected or transmitted field*/
		if (ray.reflectionMaterialId >= 0)
		{
			//double incidentAngle = AngleBetween(vecGlobal, currentCoordSys(2,Eigen::indexing::all), false, true);
			//double relPerm = _materials[ray.reflectionMaterialId].relPermittivityRe;
			//double sigma = _materials[ray.reflectionMaterialId].relConductivity;
			totalField = ComputeRefcField(vecGlobal, incidentField, ray.reflectionMaterialId, currentCoordSys, nextCoordSys); // local spherical

			//totalField = ComputeRefcField(incidentField, relPerm, sigma, _frequency, incidentAngle, 0, _useFresnelCoeff); // local spherical
		}
		else if (ray.penetrationMaterialId >= 0)
		{
			// totalField = ComputeTransField(incidentField, relPerm, sigma, frequency, incidentAngle, 0, useFresnelCoeff);
			totalField = ComputeTransField(vecGlobal, incidentField, ray.penetrationMaterialId, currentCoordSys, nextCoordSys); // local spherical
		}
		else
		{
			totalField = incidentField;
		}
		//std::cout << "Total Field Before Rotation: " << totalField << std::endl;
		// Need to do: add a case where receiver gain can be applied here
		totalField = SphericalToCartesianVector(totalField, theta, phi);
		totalField = RotateToNewCoordSys(totalField, currentCoordSys, globalCoordSys);
		//std::cout << "Total Field After Rotation: " << totalField << std::endl;

	}
	totalField = totalField * exp(-j * k * totalPathLength) / totalPathLength;

	cout << endl;
	cout << "\tDelay [ns]: " << (totalPathLength/ SPEED_OF_LIGHT)*1e9 << endl;
	cout << "\tField Components [Ex, Ey, Ez]: " << totalField.transpose().format(CleanFmt) << endl;
	cout << "\tField Strength [dBuvm]: " << 20 * log10(totalField.norm() * 1e6) << endl;
	cout << "[Leaving] FieldForPath..." << endl;

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

Vec3c FieldCompute::ComputeRefcField(const Vec3& vecGlobal, const Vec3c& efield_i_sph, int materialId, const Mat3& currentCoordSys, const Mat3& nextCoordSys)
{
	// Material property of the impinging medium
	cdouble epsilon_t;
	epsilon_t.real(_materials[materialId].relPermittivityRe);
	epsilon_t.imag(-1 * _materials[materialId].relConductivity/ (E0 * 2 * PI * frequency));

	// Permittivity of air
	cdouble epsilon_i(1, 0);
	
	// Get theta_i and theta_t
	double theta_i = ConstrainAngleTo90(AngleBetween(vecGlobal, currentCoordSys(2, Eigen::indexing::all), false, true));
	cdouble theta_t = GetTransAngle(theta_i, epsilon_i, epsilon_t);

	cdouble refTE, refTM, transTE, transTM;
	GetTECoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTE, transTE);
	GetTMCoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTM, transTM);

	// x is r, y is theta, z is phi
	cdouble Er = efield_i_sph[0];
	cdouble Etheta = efield_i_sph[1];
	cdouble Ephi = efield_i_sph[2];

	Vec3c refcField (Er, Etheta * refTM, Ephi * refTE); 
	return refcField;
}

Vec3c FieldCompute::ComputeTransField(const Vec3& vecGlobal, const Vec3c& efield_i_sph, int materialId, const Mat3& currentCoordSys, const Mat3& nextCoordSys)
{
	// Material property of the impinging medium
	cdouble epsilon_t;
	epsilon_t.real(_materials[materialId].relPermittivityRe);
	epsilon_t.imag(-1 * _materials[materialId].relConductivity / (E0 * 2 * PI * frequency));

	// Permittivity of air
	cdouble epsilon_i(1, 0);

	// Get theta_i and theta_t
	double theta_i = ConstrainAngleTo90(AngleBetween(vecGlobal, currentCoordSys(2, Eigen::indexing::all), false, true));
	cdouble theta_t = GetTransAngle(theta_i, epsilon_i, epsilon_t);

	cdouble refTE, refTM, transTE, transTM;
	GetTECoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTE, transTE);
	GetTMCoeff(theta_i, theta_t, epsilon_i, epsilon_t, _materials[materialId].width, useFresnelCoeff, refTM, transTM);

	// x is r, y is theta, z is phi
	cdouble Er = efield_i_sph[0];
	cdouble Etheta = efield_i_sph[1];
	cdouble Ephi = efield_i_sph[2];

	Vec3c transField(Er, Etheta * transTM, Ephi * transTE);
	return transField;
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

	//  cdouble gamma_te = ((n_i * cos_i) - (n_t * cos_t)) / ((n_i * cos_i) + (n_t * cos_t));
	cdouble gamma_te = -((eta_t / cos_t) - (eta_i / cos_i)) / ((eta_t / cos_t) + (eta_i / cos_i));
	cdouble tau_te = (2.0 * eta_t * cos_i) / (eta_t * cos_i + eta_i * cos_t);
	if (inf_wall)
	{
		ref_coeff = gamma_te;
		// tran_coeff = 1.0 - gamma_te;
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

	// cdouble gamma_tm = (n_t * cos_i - n_i * cos_t) / (n_i * cos_t + n_t * cos_i);
	cdouble gamma_tm = (-eta_t * cos_t + eta_i * cos_i) / (eta_t * cos_t + eta_i * cos_i);
	cdouble tau_tm = (2.0 * eta_t * cos_i) / (eta_t * cos_t + eta_i * cos_i);
	if (inf_wall)
	{
		ref_coeff = gamma_tm;
		// tran_coeff = 1.0 - ref_coeff;
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

	// DEBUGGING
	if (hitSurfId == 0 || hitSurfId == 1)
	{
		return Mat3{ {0,0,1},{1,0,0},{0,1,0} };
	}
	else if (hitSurfId == 2 || hitSurfId == 3)
	{
		return Mat3{ {1,0,0},{0,0,1},{0,-1,0} };
	}
	else
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

	Vec3 zw = normal;
	Vec3 yw = rayDir.cross(normal) / sin(theta_i);
	Vec3 xw = yw.cross(zw);

	Mat3 surfaceCoord;
	surfaceCoord(0, Eigen::indexing::all) = xw;
	surfaceCoord(1, Eigen::indexing::all) = yw;
	surfaceCoord(2, Eigen::indexing::all) = zw;

	return surfaceCoord;
}