#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include "common/VecMatDef.h"
#include "common/Global.h"
#include "common/DataStructures.h"
#include "Triangle.h"
#include "StlReader.h"
#include "timer/Timer.h"

// C4267 Warning is suppressed: size_t to int possible loss of data
#pragma warning(push)
#pragma warning(disable: 4267)

namespace Preprocessor
{
	// Functions to read files
	bool ReadPatternFile(std::string fileName, GainMap& output);
	bool ReadLocationFile(std::string fileName, VecVec3& output);
	bool ReadMaterialsFile(std::string fileName, Materials& output);
	bool StlToGeometry
	(
		std::string fileName,
		std::vector<Triangle*>& output,
		bool saveEdges, 
		bool saveFaces,
		std::string dataDir
	);

	// Functions for processing coplanar IDs and edge connectivity
	void InsertEdgeIntoMap
	(
		const size_t& v1, 
		const size_t& v2, 
		const size_t& v3, 
		const size_t& triId,
		const std::vector<Triangle*>& triangles,
		std::map<std::pair<int, int>, 
		std::vector<int>>& edgeMap,
		std::vector<std::vector<int>>& adjacencyList
	);

	void BfsCoplanarSurface(const std::vector<std::vector<int>>& adj, const std::vector<Triangle*>& triangles);
	void BFSUtil
	(
		int u, 
		const std::vector<std::vector<int>>& adj, 
		std::vector<bool>& visited, 
		const std::vector<Triangle*>& triangles,
		int coplanarId
	);

	// Functions to generate data [to be implemented]
	void GenerateRxPlane(double xMin, double yMin, double xMax, double yMax, double height, double resolution, VecVec3& output);
	void GenerateRXLine(Vec3 start, Vec3 end, double resolution, VecVec3& output);
	void GenerateAntennaPattern(std::string type, double resolution);
	
	// Functions to save data to VTK files
	bool SaveEdgesAsVtk
	(
		std::string fileName,
		const std::map<std::pair<int, int>,std::vector<int>>& edgeMap,
		const std::vector<float>& coords
	);
	
	bool SaveFacesAsVtk
	(
		std::string fileName,
		const std::vector<Triangle*>& triangles,
		const std::vector<float>& coords,
		const std::vector<unsigned int> tris
	);
}