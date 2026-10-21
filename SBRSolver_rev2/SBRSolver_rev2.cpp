// SBRSolver_rev2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <new>
#include <string>
#include <vector>
#include "BVH/BVH.h"
#include "common/VecMatDef.h"
#include "common/DataStructures.h"
#include "Preprocessor.h"
#include "RTSolver.h"
#include "SimulationConfig.h"
#include "timer/Timer.h"
#include "Triangle.h"
#include "common/Profiling.h"
#include "viewer/SceneData.h"
#if defined(SBR_HAS_VIEWER)
#include "viewer/Viewer.h"
#endif

// Directory that holds the reference configuration files, scenes and receiver
// files. The build system points this at <repo>/SBRSolver_rev2/data so the
// executable can be started from any working directory.
#ifndef SBR_DEFAULT_DATA_DIR
#define SBR_DEFAULT_DATA_DIR "./data"
#endif

#define SBR_DEFAULT_CONFIG_NAME "ground.ini"

namespace fs = std::filesystem;

static void PrintUsage(const char* program)
{
	std::cout
		<< "Usage: " << program << " [config_file] [output_dir] [options]\n"
		<< "  config_file  simulation settings (INI format, see SBRSolver_rev2/data/ground.ini)\n"
		<< "               (default: " << SBR_DEFAULT_DATA_DIR "/" SBR_DEFAULT_CONFIG_NAME << ")\n"
		<< "  output_dir   directory for the generated VTK/CSV files\n"
		<< "               (default: <config_file directory>/output, created if missing)\n"
		<< "Options:\n"
		<< "  --view                 open the built-in 3D viewer on the results\n"
		<< "                         (same as viewer = true in [outputs])\n"
		<< "  --screenshot FILE      render the viewer once into FILE (.png, .jpg or .jpeg) and exit\n"
		<< "  --viewer-frames N      close the viewer after N frames, N >= 1 (tests)\n"
		<< "  --viewer-backend NAME  openGL3_glfw (window), openGL3_egl (headless, Linux)\n"
		<< "                         or openGL_mock (no rendering); default: automatic\n"
		<< "File names inside the configuration are relative to the configuration file.\n";
}

// Command line: positional config_file and output_dir plus the options above.
struct CommandLine
{
	std::vector<std::string> positional;
	bool view = false;
	std::string screenshotFile;
	long viewerFrames = -1;
	std::string viewerBackend;
};

static bool ParseCommandLine(int argc, char* argv[], CommandLine& cl, std::string& error)
{
	for (int i = 1; i < argc; i++)
	{
		const std::string arg = argv[i];
		auto value = [&](std::string& out)
		{
			if (i + 1 >= argc) { error = arg + " needs a value"; return false; }
			out = argv[++i];
			return true;
		};
		if (arg == "--view") cl.view = true;
		else if (arg == "--screenshot")
		{
			if (!value(cl.screenshotFile)) return false;
			// The viewer writes PNG or JPEG; anything else would only fail
			// once the window is open.
			std::string extension = fs::path(cl.screenshotFile).extension().string();
			std::transform(extension.begin(), extension.end(), extension.begin(),
			               [](unsigned char c) { return char(std::tolower(c)); });
			if (extension != ".png" && extension != ".jpg" && extension != ".jpeg")
			{
				error = "--screenshot needs a .png, .jpg or .jpeg file name, got '" + cl.screenshotFile + "'";
				return false;
			}
			cl.view = true;
		}
		else if (arg == "--viewer-frames")
		{
			std::string text;
			if (!value(text)) return false;
			// The whole value must be a positive integer.
			size_t parsed = 0;
			long frames = -1;
			try { frames = std::stol(text, &parsed); }
			catch (const std::exception&) { parsed = 0; }
			if (parsed == 0 || parsed != text.size() || frames < 1)
			{
				error = "--viewer-frames needs a positive integer, got '" + text + "'";
				return false;
			}
			cl.viewerFrames = frames;
			cl.view = true;
		}
		else if (arg == "--viewer-backend")
		{
			if (!value(cl.viewerBackend)) return false;
			// Polyscope's backend names; a typo must not cost a whole run.
			if (cl.viewerBackend != "auto" && cl.viewerBackend != "openGL3_glfw" && cl.viewerBackend != "openGL3_egl"
				&& cl.viewerBackend != "openGL_mock")
			{
				error = "--viewer-backend must be openGL3_glfw, openGL3_egl, openGL_mock or auto, got '" + cl.viewerBackend + "'";
				return false;
			}
			if (cl.viewerBackend == "auto") cl.viewerBackend.clear();
			cl.view = true;
		}
		else if (arg.size() > 1 && arg[0] == '-' && arg != "-")
		{
			error = "unknown option " + arg;
			return false;
		}
		else cl.positional.push_back(arg);
	}
	if (cl.positional.size() > 2)
	{
		error = "too many arguments";
		return false;
	}
	return true;
}

static std::string PathString(const fs::path& p)
{
	// Use forward slashes everywhere; every supported OS accepts them.
	return p.generic_string();
}

// Builds the receiver positions described by the configuration.
static bool GenerateReceivers(const ReceiverConfig& rx, VecVec3& output)
{
	switch (rx.layout)
	{
	case ReceiverLayout::Line:
		Preprocessor::GenerateRXLine(rx.lineStart, rx.lineEnd, rx.resolution, output);
		break;
	case ReceiverLayout::Plane:
		Preprocessor::GenerateRxPlane(rx.xMin, rx.yMin, rx.xMax, rx.yMax, rx.height, rx.resolution, output);
		break;
	case ReceiverLayout::File:
		if (!Preprocessor::ReadLocationFile(rx.file, output))
		{
			std::cerr << "Error: cannot read receiver location file " << rx.file << std::endl;
			return false;
		}
		break;
	}
	if (output.empty())
	{
		std::cerr << "Error: the [receivers] section produced no receiver positions" << std::endl;
		return false;
	}
	return true;
}

int main(int argc, char* argv[])
{
	for (int i = 1; i < argc; i++)
	{
		const std::string arg = argv[i];
		if (arg == "-h" || arg == "--help")
		{
			PrintUsage(argv[0]);
			return 0;
		}
	}

	CommandLine cl;
	std::string clError;
	if (!ParseCommandLine(argc, argv, cl, clError))
	{
		std::cerr << "Error: " << clError << std::endl;
		PrintUsage(argv[0]);
		return 1;
	}
	const fs::path configPath = (cl.positional.size() > 0) ? fs::path(cl.positional[0]) : fs::path(SBR_DEFAULT_DATA_DIR) / SBR_DEFAULT_CONFIG_NAME;
	if (configPath.empty() || (cl.positional.size() > 1 && cl.positional[1].empty()))
	{
		std::cerr << "Error: empty path argument" << std::endl;
		PrintUsage(argv[0]);
		return 1;
	}

	std::cout << "[Entering] main ..." << std::endl;
	std::cout << "Loading simulation parameters from " << PathString(configPath) << std::endl;

	SimulationConfig config;
	std::string configError;
	if (!LoadSimulationConfig(PathString(configPath), config, configError))
	{
		std::cerr << "Error: invalid configuration:\n" << configError;
		PrintUsage(argv[0]);
		return 1;
	}

	// config.configFile is already absolute and normalised.
	const fs::path outputDir = (cl.positional.size() > 1) ? fs::path(cl.positional[1]) : fs::path(config.configFile).parent_path() / "output";
#if defined(SBR_HAS_VIEWER)
	const bool openViewer = cl.view || config.viewer;
#else
	if (cl.view)
	{
		std::cerr << "Error: this build has no viewer (configure with -DSBR_BUILD_VIEWER=ON)" << std::endl;
		return 1;
	}
	if (config.viewer)
	{
		std::cerr << "Warning: [outputs] viewer = true, but this build has no viewer (configure with -DSBR_BUILD_VIEWER=ON)" << std::endl;
	}
#endif

	std::error_code ec;
	fs::create_directories(outputDir, ec);
	if (ec)
	{
		std::cerr << "Error: cannot create output directory " << PathString(outputDir) << ": " << ec.message() << std::endl;
		return 1;
	}

	// Output filenames (fixed, so post-processing scripts and the CI reference
	// comparison keep working)
	const std::string rxLocationOutputFileName = PathString(outputDir / "RxLocations.vtk");
	const std::string rayPathFileName = PathString(outputDir / "RayPath.vtk");
	const std::string icosahedronFileName = PathString(outputDir / "Icosahedron.vtk");
	const std::string csvFieldFileName = PathString(outputDir / "ElectricField.csv");

	std::cout << "\tOutput Directory         -> " << PathString(outputDir) << std::endl;
	PrintSimulationConfig(config);

	Timer timer;

	std::cout << "[Entering] Preprocessor" << std::endl;
	// Preprocessor Begin
	VecVec3 rxLocation;
	std::vector<Triangle*> triangle_mesh;
	BVH<Triangle> bvh;

	timer.start();
	{
	SBR_ZONE_NAMED("Preprocessor");
	if (!config.patternFile.empty())
	{
		// The measured pattern is validated here but the field computation still
		// uses the analytic antenna_type: applying a tabulated gain is not
		// implemented yet (see FieldCompute::FieldForPath).
		GainMap txPattern;
		if (!Preprocessor::ReadPatternFile(config.patternFile, txPattern) || txPattern.empty())
		{
			std::cerr << "Error: cannot read transmitter pattern file " << config.patternFile << std::endl;
			return 1;
		}
		std::cout << "\tWarning: pattern_file was read (" << txPattern.size()
			<< " samples) but is not applied; the analytic antenna_type is used instead." << std::endl;
	}
	if (!GenerateReceivers(config.receivers, rxLocation))
	{
		return 1;
	}
	std::vector<Preprocessor::StlSource> sources;
	sources.push_back({ config.stlFile, config.sceneMaterialId });
	for (const GeometryConfig& g : config.geometries)
	{
		sources.push_back({ g.file, g.materialId });
	}
	if (!Preprocessor::StlToGeometry(sources, config.surfaceMaterials, triangle_mesh, config.saveEdges, config.saveFaces, PathString(outputDir) + "/") || triangle_mesh.empty())
	{
		std::cerr << "Error: failed to load the scene geometry" << std::endl;
		return 1;
	}
	bvh.ConstructBVH(triangle_mesh);
	// The solver looks surfaces up as triangle_mesh[ray.hitSurfaceID], so the
	// list must still be in creation order after the BVH is built.
	for (size_t i = 0; i < triangle_mesh.size(); i++)
	{
		if (triangle_mesh[i]->triangleId != int(i))
		{
			std::cerr << "Error: internal error, triangle list is not in id order after BVH construction" << std::endl;
			return 1;
		}
	}
	}
	// Preprocessor End

	std::cout << "\tReceivers                -> " << rxLocation.size() << std::endl;
	std::cout << "\tTriangles                -> " << triangle_mesh.size() << std::endl;
	{
		int surfaces = 0;
		std::vector<size_t> perMaterial(config.materials.size(), 0);
		for (const Triangle* t : triangle_mesh)
		{
			surfaces = std::max(surfaces, t->coplanarId + 1);
			perMaterial[t->materialId]++;
		}
		std::cout << "\tCoplanar surfaces        -> " << surfaces << std::endl;
		for (size_t i = 0; i < perMaterial.size(); i++)
		{
			std::cout << "\tMaterial " << i << " (" << config.materials[i].name << ") -> " << perMaterial[i] << " triangles" << std::endl;
		}
	}
	std::cout << "\tTotal Time in Preprocessor -> " << timer.getTime() << std::endl;
	std::cout << "[Leaving] Preprocessor" << std::endl;

	// Materials in file order; triangles carry the index they were assigned.
	std::vector<MaterialProperties> materials;
	for (const MaterialConfig& m : config.materials)
	{
		materials.push_back(m.properties);
	}

	RTSolver* rayTracer = new RTSolver();
	rayTracer->Init(materials.data(), int(materials.size()), bvh);
	rayTracer->SetThreadCount(config.threads);

	int totalPaths = 0;
	try
	{
		totalPaths = rayTracer->ExecuteRayTracing
		(
			config.txPosition,
			config.txCoordSys,
			config.maxReflections,
			config.maxTransmissions,
			config.frequencyHz,
			config.txPowerWatts,
			config.useFresnelCoefficients,
			int(config.antennaType),
			rxLocation,
			triangle_mesh,
			config.tessellation
		);
	}
	catch (const std::bad_alloc&)
	{
		std::cerr << "Error: out of memory while tracing " << rxLocation.size() << " receivers x "
			<< (10 * (1 << (2 * config.tessellation)) + 2) << " launch rays; reduce the receiver count, "
			<< "the tessellation or the number of threads" << std::endl;
		return 1;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}

	std::cout << "Total ray paths traced -> " << totalPaths << std::endl;

	{
		SBR_ZONE_NAMED("Outputs");
		// rayTracer->CmdLineDebug();
		if (config.writeRayPaths) rayTracer->SavePathsAsVtk(rayPathFileName);
		if (config.writeIcosahedron) rayTracer->SaveIcosahedronAsVtk(icosahedronFileName, config.txPosition, config.tessellation);
		if (config.writeReceivers) rayTracer->SaveReceiversAsVtk(rxLocationOutputFileName);
		if (config.writeFieldCsv) rayTracer->SaveFieldAsCsv(csvFieldFileName);
	}

#if defined(SBR_HAS_VIEWER)
	// Everything the viewer needs is copied out before the solver's memory is
	// released, so a large run does not hold both at once.
	sbr::viewer::SceneData scene;
	if (openViewer)
	{
		SBR_ZONE_NAMED("CollectSceneData");
		std::vector<std::string> materialNames;
		for (const MaterialConfig& m : config.materials) materialNames.push_back(m.name);
		scene.title = PathString(configPath.filename());
		scene.outputDir = PathString(outputDir);
		scene.mesh = sbr::viewer::BuildSceneMesh(triangle_mesh, materialNames);
		scene.edges = sbr::viewer::BuildSceneEdges(scene.mesh);
		scene.txPosition = config.txPosition;
		scene.txCoordSys = config.txCoordSys;
		scene.antennaType = AntennaTypeName(config.antennaType);
		scene.frequencyHz = config.frequencyHz;
		scene.txPowerWatts = config.txPowerWatts;
		scene.receivers.assign(rxLocation.begin(), rxLocation.end());
		scene.fields = rayTracer->Fields();
		scene.paths = rayTracer->CollectPaths();
	}
#endif

	delete rayTracer;

	for (Triangle* triangle : triangle_mesh)
	{
		delete triangle;
	}
	triangle_mesh.clear();

#if defined(SBR_HAS_VIEWER)
	if (openViewer)
	{
		SBR_ZONE_NAMED("Viewer");
		sbr::viewer::ViewerOptions options;
		options.backend = cl.viewerBackend;
		options.screenshotFile = cl.screenshotFile;
		options.frames = cl.viewerFrames;
		std::cout << "[Entering] Viewer" << std::endl;
		const int status = sbr::viewer::RunViewer(scene, options);
		std::cout << "[Leaving] Viewer" << std::endl;
		return status;
	}
#endif
	return 0;
}
