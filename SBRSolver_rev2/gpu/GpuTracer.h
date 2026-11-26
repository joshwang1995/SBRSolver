#pragma once

#include "Tracer.h"
#include "gpu/GpuContext.h"
#include "gpu/GpuTypes.h"

namespace sbr
{
namespace gpu
{

/**
 * @brief The GPU tracer: traces every launch and captures every receiver on
 * the Vulkan device, in batches of launches that fit the device, and fills
 * the same TraceResult as the CPU tracer, so the correction, the duplicate
 * removal and the field run unchanged.
 * @param precision Single or Double, as GpuContext::Resolve returned it.
 * @note Throws std::runtime_error when the scene exceeds a device limit.
 * SBR_GPU_BATCH_BYTES and SBR_GPU_CANDIDATE_CAPACITY shrink the batches
 * and the candidate buffer, for tests.
 */
void TraceAndCapture(const TracerInput& input, GpuContext& gpu, GpuPrecision precision, TraceResult& result);

} // namespace gpu
} // namespace sbr
