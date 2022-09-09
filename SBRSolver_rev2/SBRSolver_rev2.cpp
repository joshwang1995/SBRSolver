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
Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ">>");

int main()
{
	std::cout << "[Entering] main ..." << std::endl;
	std::cout << "Initializing Simulation Parameters" << std::endl;
	
	std::string txPatternFileName = "./data/TxPatternTest.dat";
	std::string rxLocationFileName = "./data/RxLocations.dat";
	std::string stlFileName = "./data/stl_files/ground.stl";
	std::string rayPathFileName = "./data/RayPath.vtk";
	Timer timer;

	double freq = 2.45e9; // frequency
	double Pt = 1; // transmit power in Watt
	int maxReflection = 5;

#ifdef DEBUG
	std::cout << "\tTX Pattern File Name  -> " << txPatternFileName << std::endl;
	std::cout << "\tRX Location File Name -> " << rxLocationFileName << std::endl;
	std::cout << "\tSTL File Name         -> " << stlFileName << std::endl;
	std::cout << "\tFrequency             -> " << freq << std::endl;
	std::cout << "\tTransmit Power        -> " << Pt << std::endl;
	std::cout << "\tMax Reflection        -> " << maxReflection << std::endl;
#endif

	std::cout << "[Entering] Preprocessor" << std::endl;
	// Preprocessor Begin
	GainMap txPattern;
	VecVec3 rxLocation;
	std::vector<Triangle*> triangle_mesh;
	BVH<Triangle> bvh;
	
	timer.start();
	Preprocessor::ReadPatternFile(txPatternFileName, txPattern);
	Preprocessor::ReadLocationFile(rxLocationFileName, rxLocation);
	Preprocessor::StlToGeometry(stlFileName, triangle_mesh);
	bvh.ConstructBVH(triangle_mesh);
	// Preprocessor End
#ifdef DEBUG
	std::cout << "\tTotal Time in Preprocessor -> " << timer.getTime() << std::endl;
#endif

	std::cout << "[Leaving] Preprocessor" << std::endl;

	// Vec3 rayOrig{ -0.5, -12.5, 1 }; // for bahen stl file
	Vec3 rayOrig{ 0, 0, 5 };

	MaterialProperties materials[2];
	// Material 1
	materials[0].frequency = 2.4e9;
	materials[0].reflectionLoss = 2;
	materials[0].relConductivity = 58.7e6;
	materials[0].relPermittivityRe = INFINITE;
	materials[0].relPermittivityIm = 0;
	materials[0].transmissionLoss = 0;
	// Material 2
	materials[1].frequency = 2.4e9;
	materials[1].reflectionLoss = 2;
	materials[1].relConductivity = 58.7e6;
	materials[1].relPermittivityRe = INFINITE;
	materials[1].relPermittivityIm = 0;
	materials[1].transmissionLoss = 10;
	

	RTSolver* rayTracer = new RTSolver();
	rayTracer->Init(materials, 2, bvh);
	int total_paths = rayTracer->ExecuteRayTracing(rayOrig, 3, 3, 0);

	rayTracer->DebugFunc();
	rayTracer->SavePathsAsVtk(rayPathFileName);

	delete rayTracer;
}