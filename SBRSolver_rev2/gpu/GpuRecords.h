#pragma once

#include <cstddef>
#include <cstdint>

// C++ mirrors of the std430 records in gpu/shaders/common.glsl. `Real` is
// float for the single-precision shaders and double for the double ones;
// the static_asserts hold both layouts to the GLSL ones.
namespace sbr
{
namespace gpu
{

template <typename Real>
struct GpuSegment
{
	Real ox, oy, oz, unfoldedStart;
	Real dx, dy, dz, t;
	int32_t parent, launch, triangle, coplanar, material, kind, reflections, transmissions;
};

template <typename Real>
struct GpuNode
{
	Real bmin[3];
	Real bmax[3];
	int32_t left, right, first, count;
};

template <typename Real>
struct GpuTriangle
{
	Real v1[3];
	Real e1[3];
	Real e2[3];
	Real n[3];
	int32_t id, coplanar, material, pad;
};

// The kinds of a segment, as the shaders write them.
const int32_t kKindRoot = 0;
const int32_t kKindReflect = 1;
const int32_t kKindTransmit = 2;

static_assert(sizeof(GpuSegment<float>) == 64, "single-precision segment record");
static_assert(sizeof(GpuSegment<double>) == 96, "double-precision segment record");
static_assert(offsetof(GpuSegment<float>, dx) == 16 && offsetof(GpuSegment<float>, parent) == 32, "segment offsets (float)");
static_assert(offsetof(GpuSegment<double>, dx) == 32 && offsetof(GpuSegment<double>, parent) == 64, "segment offsets (double)");
static_assert(offsetof(GpuSegment<float>, transmissions) == 60 && offsetof(GpuSegment<double>, transmissions) == 92, "segment offsets");
static_assert(sizeof(GpuNode<float>) == 40 && sizeof(GpuNode<double>) == 64, "BVH node record");
static_assert(offsetof(GpuNode<float>, left) == 24 && offsetof(GpuNode<double>, left) == 48, "BVH node offsets");
static_assert(sizeof(GpuTriangle<float>) == 64 && sizeof(GpuTriangle<double>) == 112, "triangle record");
static_assert(offsetof(GpuTriangle<float>, id) == 48 && offsetof(GpuTriangle<double>, id) == 96, "triangle offsets");

} // namespace gpu
} // namespace sbr
