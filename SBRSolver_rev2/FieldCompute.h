#pragma once
#include <complex>
#include "Ray.h"
#include "common/Constants.h"
#include "common/DataStructures.h"
#include "common/VecMatDef.h"
#include "RTSolver.h"

// Computes the electric field that every captured (image-corrected) path
// delivers to a receiver, with the conventions of the RayTubeTracing
// reference solver:
//
//  * the source radiates E = sqrt(P eta0 / 4 pi) * pattern(theta) * theta_hat
//    at unit distance, so |E|^2 / eta0 is the power density (RMS phasors);
//  * every path contributes exp(-j k L) / L with L the unfolded path length
//    (point source, planar surfaces: the astigmatic spreading factors of the
//    reference multiply to 1 / L);
//  * reflections and transmissions act on the TE/TM components of the field
//    in the plane of incidence, with the Fresnel coefficients of a half space
//    (use_fresnel_coefficients = true) or of a perfect conductor (false).
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
		bool useFresnelCoeff,
		int antennaType
	)
		:
		_rayPaths(rayPaths),
		_triangleMesh(&triangleMesh),
		_receivers(&receivers),
		_pathsCount(pathsCount),
		_materials(materials),
		_txCoordSys(txCoordSys),
		_frequency(frequency),
		_txPower(pt),
		_useFresnelCoeff(useFresnelCoeff),
		_antennaType(antennaType)
	{};
	~FieldCompute();
	Vec3c FieldAtReceiver(int receiverId);

	// Field of one path (the sum over a receiver's paths is FieldAtReceiver).
	Vec3c FieldForPath(const std::vector<Ray>& path) const;

	// Amplitude pattern sqrt(gain(theta)) of the analytic antennas (see
	// AntennaType in SimulationConfig.h); the field is theta-polarised.
	static double PatternAmplitude(int antennaType, double theta);

	// Complex relative permittivity eps_r - j sigma / (omega eps0) of a material.
	static cdouble RelativePermittivity(const MaterialProperties& material, double frequency);

	// Fresnel reflection coefficients of a half space for the perpendicular
	// (TE) and parallel (TM) polarisations, incidence angle theta from the
	// normal, in the reference solver's sign convention (both -1 at grazing
	// incidence; at normal incidence rPerp = (1 - sqrt(eps)) / (1 + sqrt(eps))
	// and rPar = -rPerp).
	static void FresnelReflection(double theta, cdouble eps, cdouble& rPerp, cdouble& rPar);
	static void FresnelTransmission(double theta, cdouble eps, cdouble& tPerp, cdouble& tPar);

protected:
	Paths** _rayPaths;
	std::vector<Triangle*>* _triangleMesh;
	VecVec3* _receivers;
	int _pathsCount;
	MaterialProperties* _materials;
	Mat3 _txCoordSys;
	double _frequency;
	double _txPower;
	bool _useFresnelCoeff;
	int _antennaType; // see AntennaType in SimulationConfig.h

private:
	// Source field (global coordinates) at unit distance in direction `dir`.
	Vec3c SourceField(const Vec3& dir) const;
	// Field after a specular reflection / a transmission at a surface with
	// normal `n`; `dirBack` points from the interaction point back to the
	// previous point of the path.
	Vec3c Reflect(const Vec3c& field, const MaterialProperties& material, const Vec3& dirBack, const Vec3& n) const;
	Vec3c Transmit(const Vec3c& field, const MaterialProperties& material, const Vec3& dirBack, const Vec3& n) const;
};
