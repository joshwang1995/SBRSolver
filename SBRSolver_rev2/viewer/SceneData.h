#pragma once

// Everything the built-in viewer shows, gathered from a finished run and kept
// independent of the solver's live objects and of the rendering library, so
// the solver can release its memory before the window opens and the data can
// be unit-tested without a display.
#include <array>
#include <complex>
#include <string>
#include <vector>

#include "common/VecMatDef.h"

class Triangle;

namespace sbr
{
namespace viewer
{

// The scene geometry with welded vertices. Face i is triangle i of the solver
// (the FaceID of Faces.vtk), so coplanarId / materialId line up with it.
struct SceneMesh
{
	std::vector<Vec3> vertices;
	std::vector<std::array<int, 3>> faces;
	std::vector<int> coplanarId;             // per face
	std::vector<int> materialId;             // per face
	std::vector<std::string> materialNames;  // index = material id
};

// An edge of the welded mesh with the faces that share it.
struct SceneEdge
{
	int v0 = 0;
	int v1 = 0;
	int face0 = -1;               // first two faces sharing the edge (-1: none)
	int face1 = -1;
	int faceCount = 0;            // 1: boundary edge, 2: regular edge, more: non-manifold
	double dihedralAngleDeg = 0;  // angle between the two face normals (0 or 180: flat)
	bool wedge = false;           // two faces that are not coplanar: a diffraction candidate
	bool diffracting = false;     // set by the diffraction model (future MR); false for now
};

enum class InteractionType
{
	Reflection = 0,
	Transmission = 1
};

// One image-corrected path from the transmitter to a receiver.
struct ScenePath
{
	int receiverId = -1;
	int launchId = -1;
	std::vector<Vec3> points;                   // transmitter, interaction points ..., receiver
	std::vector<InteractionType> interactions;  // one per interior point
	Vec3c field = Vec3c::Zero();                // this path's contribution at the receiver [V/m]
	double length = 0.0;                        // unfolded length [m]
};

struct SceneData
{
	std::string title;       // configuration file name
	std::string outputDir;
	SceneMesh mesh;
	std::vector<SceneEdge> edges;

	Vec3 txPosition = Vec3::Zero();
	Mat3 txCoordSys = Mat3::Identity();  // rows: antenna x, y, z axes in global coordinates
	std::string antennaType;
	double frequencyHz = 0.0;
	double txPowerWatts = 0.0;

	std::vector<Vec3> receivers;
	std::vector<Vec3c> fields;   // total field per receiver [V/m], same order as receivers
	std::vector<ScenePath> paths;
};

// Scalar quantities derived from a complex field vector.
enum class FieldQuantity
{
	Magnitude = 0,   // |E| [V/m]
	MagnitudeDb,     // 20 log10(|E| / 1 uV/m) [dBuV/m]
	PowerDbm,        // received power of a polarisation-matched isotropic antenna [dBm]
	PhaseX,          // arg(Ex) [deg]
	PhaseY,
	PhaseZ,
	MagnitudeX,      // |Ex| [V/m]
	MagnitudeY,
	MagnitudeZ,
	Count
};

const char* FieldQuantityName(FieldQuantity quantity);   // with unit, e.g. "|E| [V/m]"
bool FieldQuantityIsPhase(FieldQuantity quantity);       // periodic quantities (-180..180 deg)

// Received power [W] of an isotropic, polarisation-matched antenna:
// |E|^2 lambda^2 / (4 pi eta0), the convention of the RayTubeTracing reference.
double ReceivedPowerWatts(const Vec3c& field, double frequencyHz);

// Value of a quantity; a zero field gives -300 for the dB quantities and 0 for
// the phases so that nothing is non-finite.
double EvaluateQuantity(const Vec3c& field, FieldQuantity quantity, double frequencyHz);
std::vector<double> EvaluateQuantity(const std::vector<Vec3c>& fields, FieldQuantity quantity, double frequencyHz);

// Builders used by the solver after a run.
SceneMesh BuildSceneMesh(const std::vector<Triangle*>& triangles, const std::vector<std::string>& materialNames);
// Edges of the welded mesh; two faces whose normals differ by more than
// `coplanarToleranceDeg` (and less than 180 - tolerance) form a wedge.
std::vector<SceneEdge> BuildSceneEdges(const SceneMesh& mesh, double coplanarToleranceDeg = 0.01);

} // namespace viewer
} // namespace sbr
