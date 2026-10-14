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
#include "timer/Timer.h"

namespace Preprocessor
{
	// Functions to read files
	bool ReadPatternFile(std::string fileName, GainMap& output);
	bool ReadLocationFile(std::string fileName, VecVec3& output);
	bool ReadMaterialsFile(std::string fileName, Materials& output);
	// One STL file and the material id its triangles receive.
	struct StlSource
	{
		std::string fileName;
		int materialId = 0;
	};

	// Loads every source into one triangle mesh (triangle ids are global across
	// the files), builds the coplanar surface ids, applies the per-surface
	// material overrides and optionally writes Edges.vtk / Faces.vtk into
	// dataDir. Returns false (after printing the reason) if a file fails to
	// load or an override names a surface that does not exist.
	bool StlToGeometry
	(
		const std::vector<StlSource>& sources,
		const std::map<int, int>& surfaceMaterials,
		std::vector<Triangle*>& output,
		bool saveEdges, 
		bool saveFaces,
		std::string dataDir
	);

	// Overrides the material of whole coplanar surfaces (key: the CoplanarID
	// field written to Faces.vtk, value: material id). Returns false and fills
	// `error` if an id does not exist in the mesh.
	bool ApplySurfaceMaterials
	(
		const std::map<int, int>& surfaceMaterials,
		std::vector<Triangle*>& triangles,
		std::string& error
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