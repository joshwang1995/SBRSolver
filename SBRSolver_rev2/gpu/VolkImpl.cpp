// The one translation unit that holds volk's function pointers. volk loads
// the Vulkan loader at run time, so the solver starts on a machine without
// Vulkan and falls back to the CPU tracer.
#define VOLK_IMPLEMENTATION
#include <volk.h>
