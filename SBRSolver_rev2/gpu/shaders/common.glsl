// Definitions shared by every compute shader of the GPU tracer. Each shader
// is compiled twice: with SBR_DOUBLE=0 `real` is a 32-bit float, with
// SBR_DOUBLE=1 a 64-bit float. The host picks the build that matches the
// device (see GpuContext::Resolve).
#if SBR_DOUBLE
#extension GL_ARB_gpu_shader_fp64 : require
#define real double
#define rvec3 dvec3
#else
#define real float
#define rvec3 vec3
#endif

// Every shader starts its push constants with these two values. The host
// runs a large dispatch in chunks and writes the first item of each chunk
// here; a shader adds them to its invocation id (see GpuContext::Dispatch).
#define SBR_DISPATCH_BASE uint firstX; uint firstY;
