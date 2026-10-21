#pragma once

// Tracy instrumentation. With SBR_HAS_TRACY (CMake option SBR_ENABLE_TRACY)
// the macros below record zones that the Tracy profiler GUI shows live once
// it connects to the running solver; otherwise they compile to nothing.
//
//   SBR_ZONE();                scope named after the enclosing function
//   SBR_ZONE_NAMED("name");    scope with an explicit name
//   SBR_FRAME_MARK();          end of a frame (the viewer's render loop)
//   SBR_THREAD_NAME("name");   name the calling thread in the profiler
//   SBR_PLOT("name", value);   plot a value over time
// Several zones may be opened in the same scope (each gets its own variable).
#if defined(SBR_HAS_TRACY)
#include <tracy/Tracy.hpp>
#define SBR_CONCAT_(a, b) a##b
#define SBR_CONCAT(a, b) SBR_CONCAT_(a, b)
#define SBR_ZONE() ZoneNamed(SBR_CONCAT(sbr_zone_, __LINE__), true)
#define SBR_ZONE_NAMED(name) ZoneNamedN(SBR_CONCAT(sbr_zone_, __LINE__), name, true)
#define SBR_FRAME_MARK() FrameMark
#define SBR_FRAME_MARK_NAMED(name) FrameMarkNamed(name)
#define SBR_THREAD_NAME(name) tracy::SetThreadName(name)
#define SBR_PLOT(name, value) TracyPlot(name, value)
#else
#define SBR_ZONE()
#define SBR_ZONE_NAMED(name)
#define SBR_FRAME_MARK()
#define SBR_FRAME_MARK_NAMED(name)
#define SBR_THREAD_NAME(name)
#define SBR_PLOT(name, value)
#endif
