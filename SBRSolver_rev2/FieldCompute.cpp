#include "FieldCompute.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
	// Directivities of the analytic antennas. The half-wave value is
	// 4 pi / integral of |cos(pi/2 cos theta) / sin theta|^2 over the sphere.
	const double kHertzianDirectivity = 1.5;
	const double kHalfWaveDirectivity = 1.6409223769845667;

	// Unit vector perpendicular to `d` when the cross product that normally
	// provides it vanished (normal incidence), as in the reference solver.
	Vec3 AnyPerpendicular(const Vec3& d)
	{
		if (d.x() == 0.0) return Vec3(0.0, -d.z(), d.y()).normalized();
		if (d.y() == 0.0) return Vec3(d.z(), 0.0, -d.x()).normalized();
		return Vec3(d.z() * d.y(), d.z() * d.x(), d.x() * -2.0 * d.y()).normalized();
	}

	// Component of a complex field along a real unit vector (no conjugation).
	cdouble Along(const Vec3c& field, const Vec3& axis)
	{
		return field.x() * axis.x() + field.y() * axis.y() + field.z() * axis.z();
	}

	double ClampedAcos(double c)
	{
		return std::acos(std::max(-1.0, std::min(1.0, c)));
	}
}

FieldCompute::FieldCompute()
	:
	_rayPaths(nullptr),
	_triangleMesh(nullptr),
	_receivers(nullptr),
	_pathsCount(0),
	_materials(nullptr),
	_txCoordSys(Mat3::Identity()),
	_frequency(0.0),
	_txPower(0.0),
	_useFresnelCoeff(true),
	_antennaType(0)
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

double FieldCompute::PatternAmplitude(int antennaType, double theta)
{
	switch (antennaType)
	{
	case 0: // isotropic
		return 1.0;
	case 1: // Hertzian dipole along the antenna z axis
		return std::sqrt(kHertzianDirectivity) * std::sin(theta);
	case 2: // half-wave dipole along the antenna z axis
	{
		const double s = std::sin(theta);
		if (s < 1e-12) return 0.0;
		return std::sqrt(kHalfWaveDirectivity) * std::cos(0.5 * PI * std::cos(theta)) / s;
	}
	default:
		throw std::invalid_argument("unknown antenna type " + std::to_string(antennaType));
	}
}

cdouble FieldCompute::RelativePermittivity(const MaterialProperties& material, double frequency)
{
	return cdouble(material.relPermittivityRe, -material.relConductivity / (E0 * TWOPI * frequency));
}

void FieldCompute::FresnelReflection(double theta, cdouble eps, cdouble& rPerp, cdouble& rPar)
{
	const double s = std::sin(theta);
	const double c = std::cos(theta);
	const cdouble q = std::sqrt(eps - s * s);
	rPerp = (c - q) / (c + q);
	rPar = (eps * c - q) / (eps * c + q);
}

void FieldCompute::FresnelTransmission(double theta, cdouble eps, cdouble& tPerp, cdouble& tPar)
{
	const double s = std::sin(theta);
	const double c = std::cos(theta);
	const cdouble q = std::sqrt(eps - s * s);
	tPerp = (c + c) / (c + q);
	tPar = (std::sqrt(eps) * 2.0) * c / (eps * c + q);
}

Vec3c FieldCompute::SourceField(const Vec3& dir) const
{
	// Direction in the antenna frame (rows of _txCoordSys are the antenna
	// axes in global coordinates).
	const Vec3 local = _txCoordSys * dir;
	const double theta = ClampedAcos(local.z());
	const double phi = std::atan2(local.y(), local.x());
	const double amplitude = PatternAmplitude(_antennaType, theta) * std::sqrt(_txPower * ETA0 / (4.0 * PI));
	const Vec3 thetaHat(std::cos(theta) * std::cos(phi), std::cos(theta) * std::sin(phi), -std::sin(theta));
	const Vec3 global = _txCoordSys.transpose() * thetaHat;
	return Vec3c(cdouble(amplitude * global.x(), 0.0), cdouble(amplitude * global.y(), 0.0), cdouble(amplitude * global.z(), 0.0));
}

Vec3c FieldCompute::Reflect(const Vec3c& field, const MaterialProperties& material, const Vec3& dirBack, const Vec3& n) const
{
	// TE direction (perpendicular to the plane of incidence), TM directions
	// of the incident and of the reflected ray.
	Vec3 perp = dirBack.cross(n);
	perp = (perp.norm() == 0.0) ? AnyPerpendicular(dirBack) : perp.normalized();
	const Vec3 parIn = perp.cross(dirBack).normalized();
	const Vec3 reflected = 2.0 * dirBack.dot(n) * n - dirBack;
	const Vec3 parOut = reflected.cross(perp).normalized();

	cdouble rPerp, rPar;
	if (_useFresnelCoeff)
	{
		FresnelReflection(ClampedAcos(dirBack.dot(n)), RelativePermittivity(material, _frequency), rPerp, rPar);
	}
	else
	{
		// perfect electric conductor
		rPerp = cdouble(-1.0, 0.0);
		rPar = cdouble(1.0, 0.0);
	}
	const cdouble te = rPerp * Along(field, perp);
	const cdouble tm = rPar * Along(field, parIn);
	return Vec3c(perp.x() * te + parOut.x() * tm, perp.y() * te + parOut.y() * tm, perp.z() * te + parOut.z() * tm);
}

Vec3c FieldCompute::Transmit(const Vec3c& field, const MaterialProperties& material, const Vec3& dirBack, const Vec3& n) const
{
	Vec3 perp = dirBack.cross(n);
	perp = (perp.norm() == 0.0) ? AnyPerpendicular(dirBack) : perp.normalized();
	const Vec3 par = perp.cross(dirBack).normalized();

	cdouble tPerp(0.0, 0.0), tPar(0.0, 0.0); // a perfect conductor transmits nothing
	if (_useFresnelCoeff)
	{
		FresnelTransmission(ClampedAcos(dirBack.dot(n)), RelativePermittivity(material, _frequency), tPerp, tPar);
	}
	const cdouble te = tPerp * Along(field, perp);
	const cdouble tm = tPar * Along(field, par);
	return Vec3c(perp.x() * te + par.x() * tm, perp.y() * te + par.y() * tm, perp.z() * te + par.z() * tm);
}

Vec3c FieldCompute::FieldForPath(const std::vector<Ray>& path) const
{
	const double k = TWOPI * _frequency / SPEED_OF_LIGHT;

	// path[i] runs from sourcePoint to targetPoint; the interaction at
	// path[i].targetPoint (on surface path[i].hitSurfaceID) is described by
	// path[i + 1]: reflectionMaterialId or penetrationMaterialId >= 0.
	Vec3c field = SourceField((path[0].targetPoint - path[0].sourcePoint).normalized());
	double totalLength = 0.0;
	for (size_t i = 0; i < path.size(); i++)
	{
		const Ray& ray = path[i];
		totalLength += (ray.targetPoint - ray.sourcePoint).norm();
		if (i + 1 == path.size()) break;

		const Ray& next = path[i + 1];
		const Vec3 dirBack = (ray.sourcePoint - ray.targetPoint).normalized();
		Vec3 n = _triangleMesh->at(ray.hitSurfaceID)->norm;
		if (n.dot(dirBack) < 0.0) n = -n; // normal towards the side the ray comes from
		if (next.reflectionMaterialId >= 0)
		{
			field = Reflect(field, _materials[next.reflectionMaterialId], dirBack, n);
		}
		else if (next.penetrationMaterialId >= 0)
		{
			field = Transmit(field, _materials[next.penetrationMaterialId], dirBack, n);
		}
		else
		{
			throw std::logic_error("ray path segment without an interaction type");
		}
	}
	return field * (std::exp(-j * k * totalLength) / totalLength);
}
