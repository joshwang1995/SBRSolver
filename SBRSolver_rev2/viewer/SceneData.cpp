#include "viewer/SceneData.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>

#include "Triangle.h"
#include "common/Constants.h"

namespace sbr
{
namespace viewer
{

namespace
{
	const double kDbFloor = -300.0;

	double Magnitude(const Vec3c& field)
	{
		return std::sqrt(std::norm(field.x()) + std::norm(field.y()) + std::norm(field.z()));
	}

	double PhaseDeg(const cdouble& value)
	{
		return (value == cdouble(0.0, 0.0)) ? 0.0 : std::arg(value) * 180.0 / PI;
	}
}

const char* FieldQuantityName(FieldQuantity quantity)
{
	switch (quantity)
	{
	case FieldQuantity::Magnitude: return "|E| [V/m]";
	case FieldQuantity::MagnitudeDb: return "|E| [dBuV/m]";
	case FieldQuantity::PowerDbm: return "Received power [dBm]";
	case FieldQuantity::PhaseX: return "Phase Ex [deg]";
	case FieldQuantity::PhaseY: return "Phase Ey [deg]";
	case FieldQuantity::PhaseZ: return "Phase Ez [deg]";
	case FieldQuantity::MagnitudeX: return "|Ex| [V/m]";
	case FieldQuantity::MagnitudeY: return "|Ey| [V/m]";
	case FieldQuantity::MagnitudeZ: return "|Ez| [V/m]";
	default: return "?";
	}
}

bool FieldQuantityIsPhase(FieldQuantity quantity)
{
	return quantity == FieldQuantity::PhaseX || quantity == FieldQuantity::PhaseY || quantity == FieldQuantity::PhaseZ;
}

double ReceivedPowerWatts(const Vec3c& field, double frequencyHz)
{
	if (frequencyHz <= 0.0) return 0.0;
	const double lambda = SPEED_OF_LIGHT / frequencyHz;
	const double e2 = std::norm(field.x()) + std::norm(field.y()) + std::norm(field.z());
	return e2 * lambda * lambda / (4.0 * PI * ETA0);
}

double EvaluateQuantity(const Vec3c& field, FieldQuantity quantity, double frequencyHz)
{
	switch (quantity)
	{
	case FieldQuantity::Magnitude:
		return Magnitude(field);
	case FieldQuantity::MagnitudeDb:
	{
		const double m = Magnitude(field);
		return (m > 0.0) ? std::max(kDbFloor, 20.0 * std::log10(m * 1e6)) : kDbFloor;
	}
	case FieldQuantity::PowerDbm:
	{
		const double p = ReceivedPowerWatts(field, frequencyHz);
		return (p > 0.0) ? std::max(kDbFloor, 10.0 * std::log10(p * 1e3)) : kDbFloor;
	}
	case FieldQuantity::PhaseX: return PhaseDeg(field.x());
	case FieldQuantity::PhaseY: return PhaseDeg(field.y());
	case FieldQuantity::PhaseZ: return PhaseDeg(field.z());
	case FieldQuantity::MagnitudeX: return std::abs(field.x());
	case FieldQuantity::MagnitudeY: return std::abs(field.y());
	case FieldQuantity::MagnitudeZ: return std::abs(field.z());
	default: return 0.0;
	}
}

std::vector<double> EvaluateQuantity(const std::vector<Vec3c>& fields, FieldQuantity quantity, double frequencyHz)
{
	std::vector<double> values;
	values.reserve(fields.size());
	for (const Vec3c& f : fields)
	{
		values.push_back(EvaluateQuantity(f, quantity, frequencyHz));
	}
	return values;
}

SceneMesh BuildSceneMesh(const std::vector<Triangle*>& triangles, const std::vector<std::string>& materialNames)
{
	SceneMesh mesh;
	mesh.materialNames = materialNames;
	// Weld vertices by exact coordinates: the STL reader already unified the
	// corners of each file (float32 values), so shared corners are identical.
	std::map<std::tuple<double, double, double>, int> index;
	auto vertexIndex = [&](const Vec3& p)
	{
		const std::tuple<double, double, double> key(p.x(), p.y(), p.z());
		std::map<std::tuple<double, double, double>, int>::const_iterator it = index.find(key);
		if (it != index.end()) return it->second;
		const int id = int(mesh.vertices.size());
		mesh.vertices.push_back(p);
		index.emplace(key, id);
		return id;
	};
	mesh.faces.reserve(triangles.size());
	mesh.coplanarId.reserve(triangles.size());
	mesh.materialId.reserve(triangles.size());
	for (const Triangle* t : triangles)
	{
		mesh.faces.push_back({ vertexIndex(t->v1), vertexIndex(t->v2), vertexIndex(t->v3) });
		mesh.coplanarId.push_back(t->coplanarId);
		mesh.materialId.push_back(t->materialId);
	}
	return mesh;
}

std::vector<SceneEdge> BuildSceneEdges(const SceneMesh& mesh, double coplanarToleranceDeg)
{
	std::vector<Vec3> normals(mesh.faces.size(), Vec3::Zero());
	for (size_t f = 0; f < mesh.faces.size(); f++)
	{
		const std::array<int, 3>& face = mesh.faces[f];
		const Vec3 n = (mesh.vertices[face[1]] - mesh.vertices[face[0]]).cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
		normals[f] = (n.norm() > 0.0) ? Vec3(n.normalized()) : Vec3::Zero();
	}

	std::map<std::pair<int, int>, SceneEdge> edges;
	for (size_t f = 0; f < mesh.faces.size(); f++)
	{
		// A zero-area face has no normal and no extent: it contributes no
		// edge, so a sliver on a real corner does not hide the wedge there.
		if (normals[f].squaredNorm() == 0.0) continue;
		const std::array<int, 3>& face = mesh.faces[f];
		for (int k = 0; k < 3; k++)
		{
			const int a = face[k];
			const int b = face[(k + 1) % 3];
			if (a == b) continue;
			const std::pair<int, int> key(std::min(a, b), std::max(a, b));
			SceneEdge& e = edges[key];
			if (e.faceCount == 0)
			{
				e.v0 = key.first;
				e.v1 = key.second;
				e.face0 = int(f);
			}
			else if (e.faceCount == 1)
			{
				e.face1 = int(f);
			}
			e.faceCount++;
		}
	}

	std::vector<SceneEdge> result;
	result.reserve(edges.size());
	for (std::map<std::pair<int, int>, SceneEdge>::value_type& kv : edges)
	{
		SceneEdge& e = kv.second;
		if (e.faceCount >= 2)
		{
			const double c = std::max(-1.0, std::min(1.0, normals[e.face0].dot(normals[e.face1])));
			e.dihedralAngleDeg = std::acos(c) * 180.0 / PI;
			e.wedge = (e.faceCount == 2) && e.dihedralAngleDeg > coplanarToleranceDeg
				&& e.dihedralAngleDeg < 180.0 - coplanarToleranceDeg;
		}
		result.push_back(e);
	}
	return result;
}

} // namespace viewer
} // namespace sbr
