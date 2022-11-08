#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include "common/VecMatDef.h"
#include "common/DataStructures.h"
#include "Triangle.h"
#include "StlReader.h"


namespace Preprocessor
{
	// Functions to read files
	bool ReadPatternFile(std::string fileName, GainMap& output);
	bool ReadLocationFile(std::string fileName, VecVec3& output);
	bool ReadMaterialsFile(std::string fileName, Materials& output);
	bool StlToGeometry(std::string fileName, std::vector<Triangle*>& output, std::string outputFileName = "");

	// Functions for processing coplanar IDs and edge connectivity
	void InsertEdgeIntoMap
	(
		const size_t& v1, 
		const size_t& v2, 
		const size_t& v3, 
		const size_t& triId, 
		std::map<std::pair<int, int>, 
		std::vector<int>>& edgeMap,
		Eigen::MatrixXi& adjMatrix
	);

	void BfsCoplanarSurface(int root, const Eigen::MatrixXi& adjMatrix, const std::vector<float>& normals);

	// Functions to generate data [to be implemented]
	void GenerateRxPlane(double xMin, double yMin, double xMax, double yMax, double height, double resolution, VecVec3& output);
	void GenerateAntennaPattern(std::string type, double resolution);
	
	// Functions to save data to VTK files
	bool SaveEdgesAsVtk
	(
		std::string fileName,
		const std::map<std::pair<int, int>,std::vector<int>>& edgeMap,
		const std::vector<float>& coords
	);
	// bool SaveFacesAsVtk(std::string fileName, const VecVec3& location);
}

int unvisitedNeighbor(int index, const Eigen::ArrayX<bool>& visited, const Eigen::MatrixXi& adjMatrix);