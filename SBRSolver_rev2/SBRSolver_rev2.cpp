// SBRSolver_rev2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <filesystem>
#include <iostream>
#include <string>
#include "BVH/BVH.h"
#include "common/VecMatDef.h"
#include "Icosahedron.h"
#include "StlReader.h"
#include "Triangle.h"
#include "Preprocessor.h"
#include "timer/Timer.h"
#include "Ray.h"
#include "RTSolver.h"

#define DEBUG 1

// Directory that holds stl_files/, the RX_*.dat receiver files and the TX pattern.
// The build system points this at <repo>/SBRSolver_rev2/data so the executable can
// be started from any working directory; it can be overridden on the command line.
#ifndef SBR_DEFAULT_DATA_DIR
#define SBR_DEFAULT_DATA_DIR "./data"
#endif

namespace fs = std::filesystem;

static void PrintUsage(const char* program)
{
	std::cout
		<< "Usage: " << program << " [data_dir] [output_dir]\n"
		<< "  data_dir    directory containing stl_files/ and the RX_*.dat files\n"
		<< "              (default: " << SBR_DEFAULT_DATA_DIR << ")\n"
		<< "  output_dir  directory for the generated VTK/CSV files\n"
		<< "              (default: <data_dir>/output, created if missing)\n";
}

static std::string PathString(const fs::path& p)
{
	// Use forward slashes everywhere; every supported OS accepts them.
	return p.generic_string();
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

	const fs::path dataDir = (argc > 1) ? fs::path(argv[1]) : fs::path(SBR_DEFAULT_DATA_DIR);
	const fs::path outputDir = (argc > 2) ? fs::path(argv[2]) : dataDir / "output";

	std::cout << "[Entering] main ..." << std::endl;
	std::cout << "Initializing Simulation Parameters" << std::endl;

	// Input filenames
	std::string stlFileName = PathString(dataDir / "stl_files" / "ground.stl");
	std::string rxLocationFileName = PathString(dataDir / "RX_Ground.dat");
	std::string txPatternFileName = PathString(dataDir / "TxPatternTest.dat");

	// Output filenames
	std::string rxLocationOutputFileName = PathString(outputDir / "RxLocations.vtk");
	std::string bvhFileName = PathString(outputDir / "BVH.vtk");
	std::string rayPathFileName = PathString(outputDir / "RayPath.vtk");
	std::string icosahedronFileName = PathString(outputDir / "Icosahedron.vtk");
	std::string edgesFileName = PathString(outputDir / "Edges.vtk");
	std::string csvFieldFileName = PathString(outputDir / "ElectricField.csv");

	if (!fs::exists(stlFileName))
	{
		std::cerr << "Error: STL file not found: " << stlFileName << std::endl;
		PrintUsage(argv[0]);
		return 1;
	}

	std::error_code ec;
	fs::create_directories(outputDir, ec);
	if (ec)
	{
		std::cerr << "Error: cannot create output directory " << PathString(outputDir) << ": " << ec.message() << std::endl;
		return 1;
	}


	Timer timer;
	double freq = 1.8e9; // frequency
	double Pt = 1; // transmit power in Watt
	int maxReflection = 1; 
	int maxTransmission = 0;

#if DEBUG
	std::cout << "\tData Directory        -> " << PathString(dataDir) << std::endl;
	std::cout << "\tOutput Directory      -> " << PathString(outputDir) << std::endl;
	std::cout << "\tTX Pattern File Name  -> " << txPatternFileName << std::endl;
	std::cout << "\tRX Location File Name -> " << rxLocationFileName << std::endl;
	std::cout << "\tSTL File Name         -> " << stlFileName << std::endl;
	std::cout << "\tFrequency             -> " << freq << std::endl;
	std::cout << "\tTransmit Power        -> " << Pt << std::endl;
	std::cout << "\tMax Reflection        -> " << maxReflection << std::endl;
	std::cout << "\tMax Transmission      -> " << maxTransmission << std::endl;
#endif

	std::cout << "[Entering] Preprocessor" << std::endl;
	// Preprocessor Begin
	// GainMap txPattern;
	VecVec3 rxLocation;
	std::vector<Triangle*> triangle_mesh;
	BVH<Triangle> bvh;

	timer.start();
	// Preprocessor::ReadPatternFile(txPatternFileName, txPattern);
	// Preprocessor::GenerateRxPlane(-18, -40, -17, 0, 1, 1, rxLocation);
	// Preprocessor::GenerateRxPlane(-1 , 0 , 1, 10, 3, 1, rxLocation); 
	//Preprocessor::GenerateRxPlane(0, 0, 1, 10, 3, 1, rxLocation);
	// Preprocessor::GenerateRxPlane(-5, -5, 5, 5, 7, 1, rxLocation);
	//Preprocessor::ReadLocationFile(rxLocationFileName, rxLocation);
	Preprocessor::GenerateRXLine(Vec3(0, 0, 3), Vec3(0, 3000, 3), 1, rxLocation);
	if (!Preprocessor::StlToGeometry(stlFileName, triangle_mesh, true, true, PathString(outputDir) + "/") || triangle_mesh.empty())
	{
		std::cerr << "Error: failed to load any triangles from " << stlFileName << std::endl;
		return 1;
	}
	bvh.ConstructBVH(triangle_mesh);
	// Preprocessor End

#if DEBUG
	std::cout << "\tTotal Time in Preprocessor -> " << timer.getTime() << std::endl;
	// bvh.SaveAsVtk(bvhFileName);
#endif

	std::cout << "[Leaving] Preprocessor" << std::endl;

	
	// Vec3 rayOrig{ -10, 0,1 }; // for bahen stl file
	// Vec3 rayOrig{ 0.835938, 4.53906, 2.5 }; // for ibwave office
	Vec3 rayOrig{ 0, 0, 5 }; // for ground
	// Vec3 rayOrig{ 0,0,2.5 }; // for corner

	MaterialProperties materials[4];
	// Material 0 [Metal 5mm] -> Default material
	materials[0].frequency = 1.8e9;
	materials[0].reflectionLoss = 0.05;
	materials[0].transmissionLoss = 221.13;
	materials[0].relConductivity = 1e6;
	materials[0].relPermittivityRe = 1;
	materials[0].relPermittivityIm = 0;

	// Material 2 [plaster board]
	materials[1].frequency = 1.8e9;
	materials[1].reflectionLoss = -1;
	materials[1].transmissionLoss = -1;
	materials[1].relConductivity = 0.03;
	materials[1].relPermittivityRe = 2.5;
	materials[1].relPermittivityIm = 0.3;

	// Material 3 [wood]
	materials[2].frequency = 1.8e9;
	materials[2].reflectionLoss = -1;
	materials[2].transmissionLoss = -1;
	materials[2].relConductivity = 0.03;
	materials[2].relPermittivityRe = 1.8;
	materials[2].relPermittivityIm = 0.3;

	// Material 4 [Glass]
	materials[3].frequency = 1.8e9;
	materials[3].reflectionLoss = -1;
	materials[3].transmissionLoss = -1;
	materials[3].relConductivity = 0.005;
	materials[3].relPermittivityRe = 6;
	materials[3].relPermittivityIm = 0.05;
	
	int tessllation = 2;
	// Mat3 txCoordSys = Mat3::Identity();
	Mat3 txCoordSys{ {0,1,0},{0,0,1},{1,0,0} };

	RTSolver* rayTracer = new RTSolver();
	rayTracer->Init(materials, 1, bvh);

	int totalPaths = rayTracer->ExecuteRayTracing
	(
		rayOrig,
		txCoordSys,
		maxReflection, 
		maxTransmission,
		freq,
		Pt,
		true,
		rxLocation, 
		triangle_mesh, 
		tessllation
	);

	std::cout << "Total ray paths traced -> " << totalPaths << std::endl;

	// rayTracer->CmdLineDebug();
	rayTracer->SavePathsAsVtk(rayPathFileName);
	// rayTracer->SaveIcosahedronAsVtk(icosahedronFileName,rayOrig, tessllation);
	rayTracer->SaveReceiversAsVtk(rxLocationOutputFileName);
	rayTracer->SaveFieldAsCsv(csvFieldFileName);

	delete rayTracer;

	for (Triangle* triangle : triangle_mesh)
	{
		delete triangle;
	}
	triangle_mesh.clear();

	return 0;
}
