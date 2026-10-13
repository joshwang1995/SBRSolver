// SBRSolver_rev2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <filesystem>
#include <iostream>
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
		<< "Usage: " << program << " [config_file] [output_dir]\n"
		<< "  config_file  simulation settings (INI format, see SBRSolver_rev2/data/ground.ini)\n"
		<< "               (default: " << SBR_DEFAULT_DATA_DIR "/" SBR_DEFAULT_CONFIG_NAME << ")\n"
		<< "  output_dir   directory for the generated VTK/CSV files\n"
		<< "               (default: <config_file directory>/output, created if missing)\n"
		<< "File names inside the configuration are relative to the configuration file.\n";
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

	const fs::path configPath = (argc > 1) ? fs::path(argv[1]) : fs::path(SBR_DEFAULT_DATA_DIR) / SBR_DEFAULT_CONFIG_NAME;
	if (configPath.empty() || (argc > 2 && std::string(argv[2]).empty()))
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
	const fs::path outputDir = (argc > 2) ? fs::path(argv[2]) : fs::path(config.configFile).parent_path() / "output";

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
	if (!Preprocessor::StlToGeometry(config.stlFile, triangle_mesh, config.saveEdges, config.saveFaces, PathString(outputDir) + "/") || triangle_mesh.empty())
	{
		std::cerr << "Error: failed to load any triangles from " << config.stlFile << std::endl;
		return 1;
	}
	bvh.ConstructBVH(triangle_mesh);
	// Preprocessor End

	std::cout << "\tReceivers                -> " << rxLocation.size() << std::endl;
	std::cout << "\tTriangles                -> " << triangle_mesh.size() << std::endl;
	std::cout << "\tTotal Time in Preprocessor -> " << timer.getTime() << std::endl;
	std::cout << "[Leaving] Preprocessor" << std::endl;

	// Material 0 is the default material of every surface; the others are
	// available for surfaces that carry a different material id.
	std::vector<MaterialProperties> materials;
	for (const MaterialConfig& m : config.materials)
	{
		materials.push_back(m.properties);
	}

	RTSolver* rayTracer = new RTSolver();
	rayTracer->Init(materials.data(), int(materials.size()), bvh);

	int totalPaths = rayTracer->ExecuteRayTracing
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

	std::cout << "Total ray paths traced -> " << totalPaths << std::endl;

	// rayTracer->CmdLineDebug();
	if (config.writeRayPaths) rayTracer->SavePathsAsVtk(rayPathFileName);
	if (config.writeIcosahedron) rayTracer->SaveIcosahedronAsVtk(icosahedronFileName, config.txPosition, config.tessellation);
	if (config.writeReceivers) rayTracer->SaveReceiversAsVtk(rxLocationOutputFileName);
	if (config.writeFieldCsv) rayTracer->SaveFieldAsCsv(csvFieldFileName);

	delete rayTracer;

	for (Triangle* triangle : triangle_mesh)
	{
		delete triangle;
	}
	triangle_mesh.clear();

	return 0;
}
