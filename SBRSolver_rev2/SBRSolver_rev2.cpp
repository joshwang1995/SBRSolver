// SBRSolver_rev2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

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

int main()
{
	std::cout << "[Entering] main ..." << std::endl;
	std::cout << "Initializing Simulation Parameters" << std::endl;

	std::string txPatternFileName = "./data/TxPatternTest.dat";
	std::string rxLocationFileName = "./data/RxLocations.dat";
	std::string rxLocationOutputFileName = "./data/output/RxLocations.vtk";
	std::string stlFileName = "./data/stl_files/ground.stl";
	std::string bvhFileName = "./data/output/BVH.vtk";
	std::string rayPathFileName = "./data/output/RayPath.vtk";
	std::string icosahedronFileName = "./data/output/Icosahedron.vtk";
	Timer timer;

	double freq = 1.8e9; // frequency
	double Pt = 1; // transmit power in Watt
	int maxReflection = 2;
	int maxTransmission = 2;

#if DEBUG
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
	GainMap txPattern;
	VecVec3 rxLocation;
	std::vector<Triangle*> triangle_mesh;
	BVH<Triangle> bvh;

	timer.start();
	Preprocessor::ReadPatternFile(txPatternFileName, txPattern);
	// Preprocessor::GenerateRxPlane(-10, -10, 10, 10, 3, 1, rxLocation);
	Preprocessor::ReadLocationFile(rxLocationFileName, rxLocation);
	Preprocessor::StlToGeometry(stlFileName, triangle_mesh);
	bvh.ConstructBVH(triangle_mesh);
	// Preprocessor End

#if DEBUG
	std::cout << "\tTotal Time in Preprocessor -> " << timer.getTime() << std::endl;
	bvh.SaveAsVtk(bvhFileName);
#endif

	std::cout << "[Leaving] Preprocessor" << std::endl;

	//Vec3 rayOrig{ -0.5, -12.5, 1 }; // for bahen stl file
	Vec3 rayOrig{ 0, 0, 5 };

	MaterialProperties materials[4];
	// Material 0 [Concrete] -> Default material
	materials[0].frequency = 1.8e9;
	materials[0].reflectionLoss = -1;
	materials[0].transmissionLoss = -1;
	materials[0].relConductivity = 0.09;
	materials[0].relPermittivityRe = 9;
	materials[0].relPermittivityIm = 0.9;

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

	int tessllation = 0;
	RTSolver* rayTracer = new RTSolver();
	rayTracer->Init(materials, 4, bvh);
	int total_paths = rayTracer->ExecuteRayTracing(rayOrig, maxReflection, maxTransmission, rxLocation, tessllation);

	rayTracer->CmdLineDebug();
	rayTracer->SavePathsAsVtk(rayPathFileName);
	rayTracer->SaveIcosahedronAsVtk(icosahedronFileName,rayOrig, tessllation);
	Preprocessor::SaveLocationAsVtk(rxLocationOutputFileName, rxLocation);

	delete rayTracer;
}