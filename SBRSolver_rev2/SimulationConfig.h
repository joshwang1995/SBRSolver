#pragma once

#include <map>
#include <string>
#include <vector>

#include "common/DataStructures.h"
#include "common/VecMatDef.h"

// Analytic transmitter patterns understood by FieldCompute::GetAnalyticEfieldPattern.
enum class AntennaType
{
	Isotropic = 0,
	HertzianDipole = 1,
	HalfWaveDipole = 2
	// 3 ("array") exists in FieldCompute but has no implementation, so it is
	// not selectable from a configuration file.
};

enum class ReceiverLayout
{
	Line,  // points from line_start to line_end every `resolution` metres
	Plane, // grid over [x_min, x_max] x [y_min, y_max] at `height`, spaced `resolution`
	File   // positions read from an RX_*.dat location file
};

struct ReceiverConfig
{
	ReceiverLayout layout = ReceiverLayout::Line;
	Vec3 lineStart = Vec3::Zero();
	Vec3 lineEnd = Vec3::Zero();
	double xMin = 0.0;
	double yMin = 0.0;
	double xMax = 0.0;
	double yMax = 0.0;
	double height = 0.0;
	double resolution = 1.0;
	std::string file; // absolute path once loaded
};

struct MaterialConfig
{
	std::string name;
	MaterialProperties properties{};
};

// An extra STL file whose triangles all get one material ([geometry] section).
struct GeometryConfig
{
	std::string file;     // absolute path once loaded
	std::string material; // name of a [material]
	int materialId = 0;   // resolved index into SimulationConfig::materials
};

// Every simulation setting that used to be hard-coded in main(), as read from
// an INI-style file (see SBRSolver_rev2/data/ground.ini for the reference).
struct SimulationConfig
{
	std::string configFile; // absolute path of the file this was loaded from

	// [scene]
	std::string stlFile; // absolute path once loaded
	std::string sceneMaterial; // name of the material of stl_file (default: the first [material])
	int sceneMaterialId = 0;
	bool saveEdges = true;
	bool saveFaces = true;

	// [geometry] (zero or more): further STL files, each with its own material
	std::vector<GeometryConfig> geometries;

	// [surface_materials]: coplanar surface id (see Faces.vtk) -> material id
	std::map<int, int> surfaceMaterials;

	// [transmitter]
	Vec3 txPosition = Vec3::Zero();
	Mat3 txCoordSys = Mat3::Identity();
	double txPowerWatts = 1.0;
	double frequencyHz = 0.0;
	int tessellation = 2;
	AntennaType antennaType = AntennaType::HertzianDipole;
	std::string patternFile; // optional, absolute path once loaded

	// [receivers]
	ReceiverConfig receivers;

	// [ray_tracing]
	int maxReflections = 1;
	int maxTransmissions = 0;
	bool useFresnelCoefficients = true;

	// [material] (one section per material, in order of material id; names are unique)
	std::vector<MaterialConfig> materials;

	// [outputs]
	bool writeRayPaths = true;
	bool writeReceivers = true;
	bool writeFieldCsv = true;
	bool writeIcosahedron = false;
};

// Loads `path` into `config`. Relative file names inside the configuration are
// resolved against the directory containing the configuration file. Returns
// false and fills `error` (one "file:line: message" per line) on any problem:
// syntax errors, unknown sections or keys, missing required keys, values of the
// wrong type or out of range, and referenced files that do not exist.
bool LoadSimulationConfig(const std::string& path, SimulationConfig& config, std::string& error);

const char* AntennaTypeName(AntennaType type);
const char* ReceiverLayoutName(ReceiverLayout layout);

// Prints the loaded settings in the same key = value form as the file.
void PrintSimulationConfig(const SimulationConfig& config);
