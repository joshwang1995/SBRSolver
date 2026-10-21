#pragma once

// Interpolation of the complex field sampled at the receivers onto regular
// grids: a 2D grid in the receivers' plane and a 3D grid over their bounding
// box. The complex components are interpolated (never the derived
// quantities), so magnitudes, powers and phases are evaluated afterwards on
// the grid with EvaluateQuantity. No rendering library is involved, so this
// is unit-tested on its own.
#include <vector>

#include "common/VecMatDef.h"

namespace sbr
{
namespace viewer
{

enum class InterpolationMethod
{
	Nearest = 0,         // value of the closest sample
	InverseDistance = 1  // Shepard weights 1 / d^power over the nearest samples
};

struct InterpolationSettings
{
	InterpolationMethod method = InterpolationMethod::InverseDistance;
	int neighbours = 8;          // samples used per grid node (inverse distance)
	double power = 2.0;          // exponent of the inverse-distance weights
	double radiusFactor = 3.0;   // only samples within radiusFactor * SampleSpacing() count; 0 = no limit
};

// Median distance from a sample to its nearest neighbour (0 with fewer than
// two samples); the natural length scale of the receiver layout.
double SampleSpacing(const std::vector<Vec3>& samples);

// Least-squares plane through the samples (principal axes): u and v span the
// plane along its largest extents, `normal` completes the right-handed frame.
struct PlaneFit
{
	bool ok = false;
	Vec3 origin = Vec3::Zero();   // centroid
	Vec3 u = Vec3::UnitX();
	Vec3 v = Vec3::UnitY();
	Vec3 normal = Vec3::UnitZ();
	double extentU = 0.0;         // sample range along u
	double extentV = 0.0;         // sample range along v
	double thickness = 0.0;       // largest distance of a sample from the plane
};
PlaneFit FitPlane(const std::vector<Vec3>& samples);

// Interpolation at arbitrary query points. `valid[i]` is 0 when no sample lies
// within the search radius of query i (its value is then the nearest sample).
void Interpolate(const std::vector<Vec3>& samples, const std::vector<Vec3c>& values, const std::vector<Vec3>& queries,
                 const InterpolationSettings& settings, std::vector<Vec3c>& out, std::vector<unsigned char>& valid);

// Regular grid of (nu x nv) nodes in a plane; node (i, j) is at index j * nu + i.
struct PlaneGrid
{
	PlaneFit plane;
	int nu = 0;
	int nv = 0;
	std::vector<Vec3> nodes;
	std::vector<Vec3c> values;
	std::vector<unsigned char> valid;
	size_t Index(int i, int j) const { return size_t(j) * size_t(nu) + size_t(i); }
};

// `resolution` nodes along the longer side of the plane (at least 2); the
// grid covers the sample extents, padded by half a sample spacing.
PlaneGrid InterpolateOnPlane(const std::vector<Vec3>& samples, const std::vector<Vec3c>& values, const PlaneFit& plane,
                             int resolution, const InterpolationSettings& settings);

// Regular grid of (nx x ny x nz) nodes over the samples' bounding box; node
// (i, j, k) is at index (k * ny + j) * nx + i.
struct VolumeGrid
{
	Vec3 boundMin = Vec3::Zero();
	Vec3 boundMax = Vec3::Zero();
	int nx = 0;
	int ny = 0;
	int nz = 0;
	std::vector<Vec3> nodes;
	std::vector<Vec3c> values;
	std::vector<unsigned char> valid;
	size_t Index(int i, int j, int k) const { return (size_t(k) * size_t(ny) + size_t(j)) * size_t(nx) + size_t(i); }
};

// `resolution` nodes along the longest side of the box (at least 2); sides
// thinner than one sample spacing get a single node (flat receiver layouts
// still give a usable grid).
VolumeGrid InterpolateInVolume(const std::vector<Vec3>& samples, const std::vector<Vec3c>& values, int resolution,
                               const InterpolationSettings& settings);

} // namespace viewer
} // namespace sbr
