#pragma once

// Tracy instrumentation. With SBR_HAS_TRACY (CMake option SBR_ENABLE_TRACY)
// the macros below record zones that the Tracy profiler GUI shows live once
// it connects to the running solver; otherwise they compile to nothing.
//
//   SBR_ZONE();                scope named after the enclosing function
//   SBR_ZONE_NAMED("name");    scope with an explicit name
//   SBR_ZONE_VALUE(v);         attach an integer to the current zone
//   SBR_FRAME_MARK();          end of a frame (the viewer's render loop)
//   SBR_THREAD_NAME("name");   name the calling thread in the profiler
#if defined(SBR_HAS_TRACY)
#include <tracy/Tracy.hpp>
#define SBR_ZONE() ZoneScoped
#define SBR_ZONE_NAMED(name) ZoneScopedN(name)
#define SBR_ZONE_VALUE(value) ZoneValue(static_cast<uint64_t>(value))
#define SBR_ZONE_TEXT(text, size) ZoneText(text, size)
#define SBR_FRAME_MARK() FrameMark
#define SBR_FRAME_MARK_NAMED(name) FrameMarkNamed(name)
#define SBR_THREAD_NAME(name) tracy::SetThreadName(name)
#define SBR_PLOT(name, value) TracyPlot(name, value)
#else
#define SBR_ZONE()
#define SBR_ZONE_NAMED(name)
#define SBR_ZONE_VALUE(value)
#define SBR_ZONE_TEXT(text, size)
#define SBR_FRAME_MARK()
#define SBR_FRAME_MARK_NAMED(name)
#define SBR_THREAD_NAME(name)
#define SBR_PLOT(name, value)
#endif
