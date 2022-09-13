#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
	bool StlToGeometry(std::string fileName, std::vector<Triangle*>& output);

	// Functions to generate data [to be implemented]
	void GenerateRxPlane(double xMin, double yMin, double xMax, double yMax, double height, double resolution);
	void GenerateAntennaPattern(std::string type, double resolution);
	
	// Functions to save data to VTK files
	bool SaveLocationAsVtk(std::string fileName, const VecVec3& location);
}