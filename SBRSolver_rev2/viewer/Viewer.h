#pragma once

// The built-in 3D viewer (Polyscope). Shows the scene with face, material and
// surface ids, the edges (wedges are the diffraction candidates), the
// transmitter, the receivers coloured by any field quantity, every ray path,
// a 2D interpolated field view in the receivers' plane and a 3D slice view of
// the interpolated field.
#include <string>

#include "viewer/SceneData.h"

namespace sbr
{
namespace viewer
{

struct ViewerOptions
{
	// Polyscope backend: "" (automatic: a window, or headless rendering
	// without a display when a screenshot or a frame count is requested),
	// "openGL3_glfw", "openGL3_egl" (headless rendering on Linux) or
	// "openGL_mock" (no rendering at all, exercises everything else; used by
	// the tests). The interactive viewer needs a window.
	std::string backend;
	// When set, the view is rendered once, saved to this file (.png, .jpg or
	// .jpeg) and the viewer exits. With the mock backend no file is written.
	std::string screenshotFile;
	// Exit after this many frames (at least 1; -1: until the window is closed).
	long frames = -1;
	// Nodes along the longer side of the 2D view and the longest side of the
	// 3D view (clamped to 2..1000 and 2..64).
	int planeResolution = 120;
	int volumeResolution = 32;
};

// Runs the viewer on a finished run; returns 0 when it exits normally.
int RunViewer(const SceneData& scene, const ViewerOptions& options);

} // namespace viewer
} // namespace sbr
