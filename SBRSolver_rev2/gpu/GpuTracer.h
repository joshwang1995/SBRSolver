#pragma once

#include "Tracer.h"
#include "gpu/GpuContext.h"
#include "gpu/GpuTypes.h"

namespace sbr
{
namespace gpu
{

/**
 * @brief The GPU tracer: traces every launch on the Vulkan device and fills
 * the same TraceResult as the CPU tracer, so the correction, the duplicate
 * removal and the field run unchanged. The capture still runs on the CPU.
 * @param precision Single or Double, as GpuContext::Resolve returned it.
 * @note Throws std::runtime_error when the scene exceeds a device limit.
 */
void TraceAndCapture(const TracerInput& input, GpuContext& gpu, GpuPrecision precision, TraceResult& result);

} // namespace gpu
} // namespace sbr
