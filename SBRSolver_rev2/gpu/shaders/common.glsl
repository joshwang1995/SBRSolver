// Definitions shared by every compute shader of the GPU tracer. Each shader
// is compiled twice: with SBR_DOUBLE=0 `real` is a 32-bit float, with
// SBR_DOUBLE=1 a 64-bit float. The host picks the build that matches the
// device (see GpuContext::Resolve). The C++ mirrors of the records below
// are in gpu/GpuRecords.h, with static_asserts on every offset.
#if SBR_DOUBLE
#extension GL_ARB_gpu_shader_fp64 : require
#define real double
#define rvec3 dvec3
#define SBR_INFINITY packDouble2x32(uvec2(0u, 0x7FF00000u))
#else
#define real float
#define rvec3 vec3
#define SBR_INFINITY uintBitsToFloat(0x7F800000u)
#endif

// Every shader starts its push constants with these two values. The host
// runs a large dispatch in chunks and writes the first item of each chunk
// here; a shader adds them to its invocation id (see GpuContext::Dispatch).
#define SBR_DISPATCH_BASE uint firstX; uint firstY;

// The CPU's EPSILON (Constants.h): the float epsilon, used in double too.
const real kEpsilon = real(1.1920928955078125e-07);
// A ray that leaves the scene keeps going this far (RayLaunch's 1E10).
const real kMissT = real(1.0e10);

// One traced ray segment: its origin, unit direction, the hit parameter t
// (kMissT for a miss) and the unfolded length from the transmitter to the
// origin. `kind` is 0 for a launch ray, 1 for a reflected and 2 for a
// transmitted child. `triangle`, `coplanar` and `material` describe the
// surface this segment hits, -1 before it is traced and for a miss.
struct Segment
{
	real ox, oy, oz, unfoldedStart;
	real dx, dy, dz, t;
	int parent, launch, triangle, coplanar, material, kind, reflections, transmissions;
};

// A node of the flat BVH (BVH.h): an inner node names its two children, a
// leaf (count > 0) a range of the leaf triangle list.
struct BvhNode
{
	real bmin[3];
	real bmax[3];
	int left;
	int right;
	int first;
	int count;
};

// A triangle as Triangle::RayIntersects uses it: a vertex, the two edges
// from it, the unit normal, and the ids the tracer records.
struct Triangle
{
	real v1[3];
	real e1[3];
	real e2[3];
	real n[3];
	int id;
	int coplanar;
	int material;
	int pad;
};
