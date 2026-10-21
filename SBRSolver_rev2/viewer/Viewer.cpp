#include "viewer/Viewer.h"
#include "viewer/FieldGrid.h"
#include "common/Profiling.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// The rendering library stays inside this file: nothing else in the solver
// includes Polyscope, Dear ImGui or glm.
#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/screenshot.h"
#include "polyscope/slice_plane.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include "polyscope/volume_mesh.h"
#include "imgui.h"
#include <glm/glm.hpp>

namespace sbr
{
namespace viewer
{

namespace
{
	// Polyscope keys structures by type and name; the names also appear in
	// its structure panel.
	const char* const kSceneMeshName = "Scene";
	const char* const kEdgesName = "Edges";
	const char* const kWedgeEdgesName = "Wedge edges";
	const char* const kDiffractingEdgesName = "Diffracting edges";
	const char* const kTransmitterName = "Transmitter";
	const char* const kTransmitterAxesName = "Transmitter axes";
	const char* const kReceiversName = "Receivers";
	const char* const kPathsName = "Ray paths";
	const char* const kPlaneViewName = "Field plane (2D)";
	const char* const kVolumeViewName = "Field volume (3D)";
	const char* const kMockBackend = "openGL_mock";

	const int kQuantityCount = static_cast<int>(FieldQuantity::Count);
	// Which side of the scene faces is hidden. Front faces are culled by
	// registering the mesh with reversed winding, which turns them into back
	// faces for the renderer; with the outward-facing normals of a modelled
	// room that removes the walls between the camera and the interior.
	const int kCullNone = 0;
	const int kCullBack = 1;
	const int kCullFront = 2;
	const char* const kFaceCullingNames[] = { "None (both sides drawn)", "Back faces", "Front faces (look into a room)" };
	// Upper limits of the grid resolutions: about a million nodes each, which
	// interpolates in seconds and fits in memory comfortably.
	// (The slice of a volume mesh needs its tetrahedra, about five per cell,
	// so the volume stays smaller.)
	const int kMaxPlaneResolution = 1000;
	const int kMaxVolumeResolution = 64;
	// EvaluateQuantity's value of the dB quantities for a zero field
	// (SceneData.h); kept out of the colour range so it does not flatten it.
	const double kZeroFieldDb = -300.0;

	glm::vec3 ToGlm(const Vec3& v)
	{
		return glm::vec3(static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z()));
	}

	std::vector<glm::vec3> ToGlm(const std::vector<Vec3>& points)
	{
		std::vector<glm::vec3> result;
		result.reserve(points.size());
		for (const Vec3& p : points)
		{
			result.push_back(ToGlm(p));
		}
		return result;
	}

	const char* MethodName(InterpolationMethod method)
	{
		return method == InterpolationMethod::Nearest ? "Nearest sample" : "Inverse distance";
	}

	// Colour map and range of a scalar quantity derived from the field: phases
	// are periodic, the dB quantities ignore the zero-field sentinel.
	template <class ScalarQuantityT>
	void ConfigureFieldScalar(ScalarQuantityT* scalar, FieldQuantity quantity, const std::vector<double>& values)
	{
		if (FieldQuantityIsPhase(quantity))
		{
			scalar->setColorMap("phase");
			scalar->setMapRange({ -180.0, 180.0 });
			return;
		}
		if (quantity != FieldQuantity::MagnitudeDb && quantity != FieldQuantity::PowerDbm)
		{
			return;
		}
		double low = std::numeric_limits<double>::infinity();
		double high = -std::numeric_limits<double>::infinity();
		for (double v : values)
		{
			if (v > kZeroFieldDb + 1.0 && std::isfinite(v))
			{
				low = std::min(low, v);
				high = std::max(high, v);
			}
		}
		if (low < high)
		{
			scalar->setMapRange({ low, high });
		}
	}

	// Registers one scalar quantity per FieldQuantity through `add(name, values)`
	// (a lambda calling the structure's addScalarQuantity flavour) and enables
	// the selected one.
	template <class AddFunction>
	void AddFieldQuantities(const std::vector<Vec3c>& fields, double frequencyHz, FieldQuantity selected, AddFunction add)
	{
		for (int q = 0; q < kQuantityCount; ++q)
		{
			const FieldQuantity quantity = static_cast<FieldQuantity>(q);
			const std::vector<double> values = EvaluateQuantity(fields, quantity, frequencyHz);
			auto* scalar = add(std::string(FieldQuantityName(quantity)), values);
			ConfigureFieldScalar(scalar, quantity, values);
			scalar->setEnabled(quantity == selected);
		}
	}

	// Checkbox bound to a structure's visibility. The state is read back from
	// Polyscope every frame so toggles made in its own structure panel show.
	void OverlayCheckbox(const char* label, polyscope::Structure* structure, bool& shown)
	{
		if (structure == nullptr)
		{
			bool off = false;
			ImGui::BeginDisabled();
			ImGui::Checkbox(label, &off);
			ImGui::EndDisabled();
			return;
		}
		shown = structure->isEnabled();
		if (ImGui::Checkbox(label, &shown))
		{
			structure->setEnabled(shown);
		}
	}

	std::string FormatFrequency(double hz)
	{
		char buffer[64];
		if (hz >= 1e9)
		{
			std::snprintf(buffer, sizeof(buffer), "%.6g GHz", hz / 1e9);
		}
		else if (hz >= 1e6)
		{
			std::snprintf(buffer, sizeof(buffer), "%.6g MHz", hz / 1e6);
		}
		else
		{
			std::snprintf(buffer, sizeof(buffer), "%.6g Hz", hz);
		}
		return buffer;
	}

	// Everything registered with Polyscope for one run, plus the panel state.
	class ViewerSession
	{
	public:
		ViewerSession(const SceneData& scene, const ViewerOptions& options);

		void RegisterAll();
		void FrameCamera();
		void BuildGui();

	private:
		void ComputeBounds();
		void RegisterSceneMesh();
		void RegisterEdges();
		void RegisterTransmitter();
		void RegisterReceivers();
		void RegisterPaths();
		void RegisterPlaneView();
		void RegisterVolumeView();
		void RegisterSlicePlane();
		void SyncSliceState();
		void SyncSliceInspection();
		void ApplySliceIgnores();
		void SelectQuantity(int quantityIndex);
		void SetSlicePose();

		void GuiRunInfo();
		void GuiOverlays();
		void GuiScene();
		void GuiPaths();
		void GuiInterpolation();
		void GuiPlaneView();
		void GuiVolumeView();

		const SceneData& _scene;
		ViewerOptions _options;
		std::vector<Vec3c> _receiverFields;   // scene.fields padded to the receiver count
		std::vector<int> _pathsPerReceiver;

		Vec3 _boundMin = Vec3::Zero();         // everything shown: mesh, receivers, transmitter
		Vec3 _boundMax = Vec3::Zero();
		double _boundDiagonal = 1.0;

		// Structures (null when not registered)
		polyscope::SurfaceMesh* _sceneMesh = nullptr;
		polyscope::CurveNetwork* _edges = nullptr;
		polyscope::CurveNetwork* _wedgeEdges = nullptr;
		polyscope::CurveNetwork* _diffractingEdges = nullptr;
		polyscope::PointCloud* _transmitter = nullptr;
		polyscope::CurveNetwork* _transmitterAxes = nullptr;
		polyscope::PointCloud* _receivers = nullptr;
		polyscope::CurveNetwork* _paths = nullptr;
		polyscope::SurfaceMesh* _planeView = nullptr;
		polyscope::VolumeMesh* _volumeView = nullptr;
		polyscope::SlicePlane* _slicePlane = nullptr;

		// Overlay visibility
		bool _showSceneMesh = true;
		bool _showEdges = false;
		bool _showWedgeEdges = true;
		bool _showDiffractingEdges = true;
		bool _showTransmitter = true;
		bool _showTransmitterAxes = true;
		bool _showReceivers = true;
		bool _showPaths = true;
		bool _showPlaneView = true;
		bool _showVolumeView = false;
		bool _showSlicePlane = false;

		// Panel state
		int _selectedQuantity = static_cast<int>(FieldQuantity::MagnitudeDb);
		// Half transparent by default so that the inside of a closed room (the
		// receivers, the paths and the field views) shows from outside.
		float _sceneTransparency = 0.5f;
		int _faceCulling = kCullNone;
		int _pathReceiverFilter = -1;
		size_t _pathsShown = 0;
		std::string _pathQuantity = "Path |E| [dBuV/m]";   // quantity shown on the paths, kept across re-registrations
		InterpolationSettings _interpolation;
		int _planeResolution = 120;
		int _volumeResolution = 32;
		PlaneFit _planeFit;
		int _planeNu = 0;
		int _planeNv = 0;
		std::string _planeStatus;
		VolumeGrid _volumeGrid;
		std::string _volumeStatus;
		bool _volumeInitialised = false;   // first registration chooses the default visibility
		float _volumeTransparency = 1.0f;
		bool _sliceCutsScene = false;
		bool _sliceDrawPlane = false;
		bool _sliceDrawWidget = true;
		glm::vec3 _sliceCentre = glm::vec3(0.0f);
		glm::vec3 _sliceNormal = glm::vec3(0.0f, 0.0f, 1.0f);
		int _sliceAxis = 2;
	};

	ViewerSession::ViewerSession(const SceneData& scene, const ViewerOptions& options)
		: _scene(scene), _options(options)
	{
		_planeResolution = std::max(2, std::min(kMaxPlaneResolution, options.planeResolution));
		_volumeResolution = std::max(2, std::min(kMaxVolumeResolution, options.volumeResolution));

		_receiverFields = scene.fields;
		_receiverFields.resize(scene.receivers.size(), Vec3c::Zero());

		_pathsPerReceiver.assign(scene.receivers.size(), 0);
		for (const ScenePath& path : scene.paths)
		{
			if (path.receiverId >= 0 && path.receiverId < static_cast<int>(_pathsPerReceiver.size()))
			{
				++_pathsPerReceiver[path.receiverId];
			}
		}
		ComputeBounds();
	}

	void ViewerSession::ComputeBounds()
	{
		bool first = true;
		auto extend = [&](const Vec3& p)
		{
			if (first)
			{
				_boundMin = p;
				_boundMax = p;
				first = false;
			}
			else
			{
				_boundMin = _boundMin.cwiseMin(p);
				_boundMax = _boundMax.cwiseMax(p);
			}
		};
		for (const Vec3& v : _scene.mesh.vertices)
		{
			extend(v);
		}
		for (const Vec3& r : _scene.receivers)
		{
			extend(r);
		}
		extend(_scene.txPosition);
		_boundDiagonal = std::max((_boundMax - _boundMin).norm(), 1e-3);
	}

	void ViewerSession::RegisterAll()
	{
		RegisterSceneMesh();
		RegisterEdges();
		RegisterTransmitter();
		RegisterReceivers();
		RegisterPaths();
		RegisterPlaneView();
		RegisterVolumeView();
		RegisterSlicePlane();
	}

	// 1. Scene mesh with face ids, material ids and surface (coplanar) ids
	void ViewerSession::RegisterSceneMesh()
	{
		const SceneMesh& mesh = _scene.mesh;
		if (mesh.vertices.empty() || mesh.faces.empty())
		{
			return;
		}
		// Re-registration (a culling change) keeps the live state: the
		// visibility and transparency may have been changed in Polyscope's own
		// panel, and the face quantity in use stays selected.
		const char* const faceQuantities[] = { "Face ID", "Material ID", "Surface ID" };
		const char* enabledQuantity = "Material ID";
		if (_sceneMesh != nullptr)
		{
			_showSceneMesh = _sceneMesh->isEnabled();
			_sceneTransparency = _sceneMesh->getTransparency();
			for (const char* name : faceQuantities)
			{
				polyscope::Quantity* quantity = _sceneMesh->getQuantity(name);
				if (quantity != nullptr && quantity->isEnabled())
				{
					enabledQuantity = name;
				}
			}
		}
		std::vector<std::array<size_t, 3>> faces;
		faces.reserve(mesh.faces.size());
		for (const std::array<int, 3>& f : mesh.faces)
		{
			if (_faceCulling == kCullFront)
			{
				faces.push_back({ size_t(f[0]), size_t(f[2]), size_t(f[1]) });
			}
			else
			{
				faces.push_back({ size_t(f[0]), size_t(f[1]), size_t(f[2]) });
			}
		}
		_sceneMesh = polyscope::registerSurfaceMesh(kSceneMeshName, ToGlm(mesh.vertices), faces);

		std::vector<int> faceId(mesh.faces.size());
		for (size_t i = 0; i < faceId.size(); ++i)
		{
			faceId[i] = static_cast<int>(i);
		}
		std::vector<int> materialId = mesh.materialId;
		materialId.resize(mesh.faces.size(), 0);
		std::vector<int> surfaceId = mesh.coplanarId;
		surfaceId.resize(mesh.faces.size(), 0);

		_sceneMesh->addFaceScalarQuantity("Face ID", faceId, polyscope::DataType::CATEGORICAL);
		_sceneMesh->addFaceScalarQuantity("Material ID", materialId, polyscope::DataType::CATEGORICAL);
		_sceneMesh->addFaceScalarQuantity("Surface ID", surfaceId, polyscope::DataType::CATEGORICAL);
		_sceneMesh->getQuantity(enabledQuantity)->setEnabled(true);
		_sceneMesh->setBackFacePolicy(_faceCulling == kCullNone ? polyscope::BackFacePolicy::Different : polyscope::BackFacePolicy::Cull);
		_sceneMesh->setTransparency(_sceneTransparency);
		_sceneMesh->setEnabled(_showSceneMesh);
		ApplySliceIgnores();
	}

	// 2. Edges, wedge edges (diffraction candidates) and diffracting edges
	void ViewerSession::RegisterEdges()
	{
		if (_scene.mesh.vertices.empty() || _scene.edges.empty())
		{
			return;
		}
		auto registerSubset = [&](const char* name, bool (*select)(const SceneEdge&), float radius, glm::vec3 colour,
		                          bool shown) -> polyscope::CurveNetwork*
		{
			// Only the vertices of the selected edges become nodes: a curve
			// network draws a sphere at every node, used by an edge or not.
			std::vector<glm::vec3> nodes;
			std::vector<int> localIndex(_scene.mesh.vertices.size(), -1);
			auto local = [&](int v)
			{
				if (localIndex[v] < 0)
				{
					localIndex[v] = static_cast<int>(nodes.size());
					nodes.push_back(ToGlm(_scene.mesh.vertices[v]));
				}
				return size_t(localIndex[v]);
			};
			std::vector<std::array<size_t, 2>> edges;
			std::vector<double> dihedral;
			std::vector<int> faceCount;
			for (const SceneEdge& e : _scene.edges)
			{
				if (select(e))
				{
					edges.push_back({ local(e.v0), local(e.v1) });
					dihedral.push_back(e.dihedralAngleDeg);
					faceCount.push_back(e.faceCount);
				}
			}
			if (edges.empty())
			{
				return nullptr;
			}
			polyscope::CurveNetwork* network = polyscope::registerCurveNetwork(name, nodes, edges);
			network->addEdgeScalarQuantity("Dihedral angle [deg]", dihedral);
			network->addEdgeScalarQuantity("Face count", faceCount, polyscope::DataType::CATEGORICAL);
			network->setRadius(radius);
			network->setColor(colour);
			network->setEnabled(shown);
			return network;
		};

		_edges = registerSubset(kEdgesName, [](const SceneEdge&) { return true; }, 0.0012f,
		                        glm::vec3(0.35f, 0.35f, 0.35f), _showEdges);
		_wedgeEdges = registerSubset(kWedgeEdgesName, [](const SceneEdge& e) { return e.wedge; }, 0.003f,
		                             glm::vec3(1.0f, 0.55f, 0.0f), _showWedgeEdges);
		// Only present once the diffraction model (a later MR) marks edges.
		_diffractingEdges = registerSubset(kDiffractingEdgesName, [](const SceneEdge& e) { return e.diffracting; },
		                                   0.004f, glm::vec3(0.9f, 0.1f, 0.9f), _showDiffractingEdges);
	}

	// 3a. Transmitter position and antenna axes
	void ViewerSession::RegisterTransmitter()
	{
		const std::vector<glm::vec3> position{ ToGlm(_scene.txPosition) };
		_transmitter = polyscope::registerPointCloud(kTransmitterName, position);
		_transmitter->setPointColor(glm::vec3(0.9f, 0.1f, 0.1f));
		_transmitter->setPointRadius(0.012);
		_transmitter->setEnabled(_showTransmitter);

		// The antenna frame: rows of txCoordSys are its x, y, z axes.
		const double length = std::max(0.08 * _boundDiagonal, 0.05);
		std::vector<glm::vec3> nodes{ ToGlm(_scene.txPosition) };
		std::vector<std::array<size_t, 2>> edges;
		std::vector<glm::vec3> colours;
		const glm::vec3 axisColours[3] = { glm::vec3(0.9f, 0.2f, 0.2f), glm::vec3(0.2f, 0.8f, 0.2f), glm::vec3(0.2f, 0.3f, 0.95f) };
		std::vector<int> axisId;
		for (int axis = 0; axis < 3; ++axis)
		{
			Vec3 direction = _scene.txCoordSys.row(axis).transpose();
			if (direction.norm() > 0.0)
			{
				direction.normalize();
			}
			nodes.push_back(ToGlm(_scene.txPosition + length * direction));
			edges.push_back({ size_t(0), size_t(axis + 1) });
			colours.push_back(axisColours[axis]);
			axisId.push_back(axis);
		}
		_transmitterAxes = polyscope::registerCurveNetwork(kTransmitterAxesName, nodes, edges);
		_transmitterAxes->addEdgeColorQuantity("Axis colour (x red, y green, z blue)", colours)->setEnabled(true);
		_transmitterAxes->addEdgeScalarQuantity("Axis", axisId, polyscope::DataType::CATEGORICAL);
		_transmitterAxes->setRadius(0.004f);
		_transmitterAxes->setEnabled(_showTransmitterAxes);
	}

	// 3b. Receivers coloured by any field quantity
	void ViewerSession::RegisterReceivers()
	{
		if (_scene.receivers.empty())
		{
			return;
		}
		_receivers = polyscope::registerPointCloud(kReceiversName, ToGlm(_scene.receivers));
		_receivers->setPointRadius(0.004);
		polyscope::PointCloud* cloud = _receivers;
		AddFieldQuantities(_receiverFields, _scene.frequencyHz, static_cast<FieldQuantity>(_selectedQuantity),
		                   [cloud](const std::string& name, const std::vector<double>& values)
		                   {
			                   return cloud->addScalarQuantity(name, values);
		                   });
		_receivers->addScalarQuantity("Path count", _pathsPerReceiver);
		_receivers->setEnabled(_showReceivers);
	}

	// 4. Ray paths: one polyline per path, filtered by receiver
	void ViewerSession::RegisterPaths()
	{
		std::vector<glm::vec3> nodes;
		std::vector<std::array<size_t, 2>> edges;
		std::vector<int> receiverId;
		std::vector<int> interactions;
		std::vector<int> launchId;
		std::vector<double> edgeLength;
		std::vector<double> magnitude;
		std::vector<double> magnitudeDb;
		std::vector<int> nodeType;   // 0 transmitter, 1 reflection, 2 transmission, 3 receiver
		_pathsShown = 0;

		// Re-registration (a filter change) keeps the visibility and the
		// quantity in use, also across a filter value without any path.
		const char* const pathQuantities[] = { "Receiver", "Interactions", "Launch ID", "Segment length [m]", "Path |E| [V/m]",
		                                       "Path |E| [dBuV/m]", "Node type (0 Tx, 1 refl, 2 trans, 3 Rx)" };
		if (_paths != nullptr)
		{
			_showPaths = _paths->isEnabled();
			for (const char* name : pathQuantities)
			{
				polyscope::Quantity* quantity = _paths->getQuantity(name);
				if (quantity != nullptr && quantity->isEnabled())
				{
					_pathQuantity = name;
				}
			}
		}

		for (const ScenePath& path : _scene.paths)
		{
			if (path.points.size() < 2)
			{
				continue;
			}
			if (_pathReceiverFilter >= 0 && path.receiverId != _pathReceiverFilter)
			{
				continue;
			}
			++_pathsShown;
			const size_t first = nodes.size();
			const int interior = static_cast<int>(path.points.size()) - 2;
			const double pathMagnitude = path.field.norm();
			const double pathMagnitudeDb = EvaluateQuantity(path.field, FieldQuantity::MagnitudeDb, _scene.frequencyHz);
			for (size_t i = 0; i < path.points.size(); ++i)
			{
				nodes.push_back(ToGlm(path.points[i]));
				magnitude.push_back(pathMagnitude);
				magnitudeDb.push_back(pathMagnitudeDb);
				int type = 1;
				if (i == 0)
				{
					type = 0;
				}
				else if (i + 1 == path.points.size())
				{
					type = 3;
				}
				else if (i - 1 < path.interactions.size())
				{
					type = path.interactions[i - 1] == InteractionType::Transmission ? 2 : 1;
				}
				nodeType.push_back(type);
			}
			for (size_t i = 0; i + 1 < path.points.size(); ++i)
			{
				edges.push_back({ first + i, first + i + 1 });
				receiverId.push_back(path.receiverId);
				interactions.push_back(interior);
				launchId.push_back(path.launchId);
				edgeLength.push_back((path.points[i + 1] - path.points[i]).norm());
			}
		}

		if (edges.empty())
		{
			if (_paths != nullptr)
			{
				polyscope::removeCurveNetwork(kPathsName);
				_paths = nullptr;
			}
			return;
		}

		// Re-registering under the same name replaces the network.
		_paths = polyscope::registerCurveNetwork(kPathsName, nodes, edges);
		_paths->addEdgeScalarQuantity("Receiver", receiverId, polyscope::DataType::CATEGORICAL);
		_paths->addEdgeScalarQuantity("Interactions", interactions, polyscope::DataType::CATEGORICAL);
		_paths->addEdgeScalarQuantity("Launch ID", launchId, polyscope::DataType::CATEGORICAL);
		_paths->addEdgeScalarQuantity("Segment length [m]", edgeLength);
		_paths->addNodeScalarQuantity("Path |E| [V/m]", magnitude);
		_paths->addNodeScalarQuantity("Path |E| [dBuV/m]", magnitudeDb);
		_paths->addNodeScalarQuantity("Node type (0 Tx, 1 refl, 2 trans, 3 Rx)", nodeType, polyscope::DataType::CATEGORICAL);
		if (polyscope::Quantity* quantity = _paths->getQuantity(_pathQuantity))
		{
			quantity->setEnabled(true);
		}
		_paths->setRadius(0.0025f);
		_paths->setEnabled(_showPaths);
		ApplySliceIgnores();
	}

	// 5. Interpolated field on the receivers' plane
	void ViewerSession::RegisterPlaneView()
	{
		if (_planeView != nullptr)
		{
			_showPlaneView = _planeView->isEnabled();
		}
		_planeView = nullptr;
		_planeNu = 0;
		_planeNv = 0;
		if (_scene.receivers.empty())
		{
			_planeStatus = "no receivers";
			return;
		}
		_planeFit = FitPlane(_scene.receivers);
		if (!_planeFit.ok)
		{
			_planeStatus = "plane fit failed";
			return;
		}
		const PlaneGrid grid = InterpolateOnPlane(_scene.receivers, _receiverFields, _planeFit, _planeResolution, _interpolation);
		_planeNu = grid.nu;
		_planeNv = grid.nv;
		if (grid.nu < 2 || grid.nv < 2 || grid.nodes.size() != size_t(grid.nu) * size_t(grid.nv))
		{
			_planeStatus = "grid too small";
			polyscope::removeSurfaceMesh(kPlaneViewName);
			return;
		}

		std::vector<std::array<size_t, 4>> quads;
		quads.reserve(size_t(grid.nu - 1) * size_t(grid.nv - 1));
		for (int j = 0; j + 1 < grid.nv; ++j)
		{
			for (int i = 0; i + 1 < grid.nu; ++i)
			{
				quads.push_back({ grid.Index(i, j), grid.Index(i + 1, j), grid.Index(i + 1, j + 1), grid.Index(i, j + 1) });
			}
		}
		_planeView = polyscope::registerSurfaceMesh(kPlaneViewName, ToGlm(grid.nodes), quads);
		polyscope::SurfaceMesh* mesh = _planeView;
		AddFieldQuantities(grid.values, _scene.frequencyHz, static_cast<FieldQuantity>(_selectedQuantity),
		                   [mesh](const std::string& name, const std::vector<double>& values)
		                   {
			                   return mesh->addVertexScalarQuantity(name, values);
		                   });
		const std::vector<int> valid(grid.valid.begin(), grid.valid.end());
		_planeView->addVertexScalarQuantity("Within radius", valid, polyscope::DataType::CATEGORICAL);
		_planeView->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
		_planeView->setEnabled(_showPlaneView);
		char buffer[160];
		std::snprintf(buffer, sizeof(buffer), "%d x %d nodes, extents %.2f x %.2f m, thickness %.3f m",
		              grid.nu, grid.nv, _planeFit.extentU, _planeFit.extentV, _planeFit.thickness);
		_planeStatus = buffer;
		ApplySliceIgnores();
	}

	// 6. Interpolated field in the receivers' bounding box, cut by a slice plane
	void ViewerSession::RegisterVolumeView()
	{
		if (_volumeView != nullptr)
		{
			_showVolumeView = _volumeView->isEnabled();
			_volumeTransparency = _volumeView->getTransparency();
		}
		_volumeView = nullptr;
		if (_scene.receivers.empty())
		{
			_volumeStatus = "no receivers";
			return;
		}
		_volumeGrid = InterpolateInVolume(_scene.receivers, _receiverFields, _volumeResolution, _interpolation);
		const VolumeGrid& grid = _volumeGrid;
		if (grid.nx < 2 || grid.ny < 2 || grid.nz < 2 || grid.nodes.size() != size_t(grid.nx) * size_t(grid.ny) * size_t(grid.nz))
		{
			char buffer[160];
			std::snprintf(buffer, sizeof(buffer), "flat receiver layout (%d x %d x %d nodes): no volume, see the 2D view",
			              grid.nx, grid.ny, grid.nz);
			_volumeStatus = buffer;
			if (_slicePlane != nullptr)
			{
				_slicePlane->setVolumeMeshToInspect("");   // before the mesh it may refer to goes
			}
			polyscope::removeVolumeMesh(kVolumeViewName);
			SyncSliceInspection();
			return;
		}

		// Hexahedra in Polyscope's (VTK) order: the four corners of the k face
		// counter-clockwise, then the same four corners of the k + 1 face.
		std::vector<std::array<size_t, 8>> hexes;
		hexes.reserve(size_t(grid.nx - 1) * size_t(grid.ny - 1) * size_t(grid.nz - 1));
		for (int k = 0; k + 1 < grid.nz; ++k)
		{
			for (int j = 0; j + 1 < grid.ny; ++j)
			{
				for (int i = 0; i + 1 < grid.nx; ++i)
				{
					hexes.push_back({ grid.Index(i, j, k), grid.Index(i + 1, j, k), grid.Index(i + 1, j + 1, k), grid.Index(i, j + 1, k),
					                  grid.Index(i, j, k + 1), grid.Index(i + 1, j, k + 1), grid.Index(i + 1, j + 1, k + 1), grid.Index(i, j + 1, k + 1) });
				}
			}
		}
		if (_slicePlane != nullptr)
		{
			_slicePlane->setVolumeMeshToInspect("");   // the mesh it refers to is replaced below
		}
		_volumeView = polyscope::registerHexMesh(kVolumeViewName, ToGlm(grid.nodes), hexes);
		polyscope::VolumeMesh* mesh = _volumeView;
		AddFieldQuantities(grid.values, _scene.frequencyHz, static_cast<FieldQuantity>(_selectedQuantity),
		                   [mesh](const std::string& name, const std::vector<double>& values)
		                   {
			                   return mesh->addVertexScalarQuantity(name, values);
		                   });
		const std::vector<int> valid(grid.valid.begin(), grid.valid.end());
		_volumeView->addVertexScalarQuantity("Within radius", valid, polyscope::DataType::CATEGORICAL);
		_volumeView->setTransparency(_volumeTransparency);
		if (!_volumeInitialised)
		{
			// A 3D receiver layout is best explored with the slice; a flat one
			// (the common case) keeps the 2D view only.
			_volumeInitialised = true;
			_showVolumeView = true;
			_showSlicePlane = true;
		}
		_volumeView->setEnabled(_showVolumeView);
		SyncSliceInspection();
		char buffer[160];
		std::snprintf(buffer, sizeof(buffer), "%d x %d x %d nodes, %zu cells", grid.nx, grid.ny, grid.nz, hexes.size());
		_volumeStatus = buffer;
	}

	void ViewerSession::RegisterSlicePlane()
	{
		_slicePlane = polyscope::addSceneSlicePlane(false);
		_slicePlane->setDrawPlane(_sliceDrawPlane);
		_slicePlane->setDrawWidget(_sliceDrawWidget);
		Vec3 centre = 0.5 * (_boundMin + _boundMax);
		if (!_scene.receivers.empty())
		{
			centre = Vec3::Zero();
			for (const Vec3& r : _scene.receivers)
			{
				centre += r;
			}
			centre /= double(_scene.receivers.size());
		}
		_sliceCentre = ToGlm(centre);
		_sliceAxis = 2;
		_sliceNormal = glm::vec3(0.0f, 0.0f, 1.0f);
		SetSlicePose();
		_slicePlane->setEnabled(_showSlicePlane);
		SyncSliceInspection();
		ApplySliceIgnores();
	}

	// Mirrors of the slice plane's live state, which the user may also change
	// in Polyscope's own panel or with the gizmo. A disabled plane reports a
	// pose at infinity, so the pose is only followed while it is enabled.
	void ViewerSession::SyncSliceState()
	{
		if (_slicePlane == nullptr)
		{
			return;
		}
		_sliceDrawPlane = _slicePlane->getDrawPlane();
		_sliceDrawWidget = _slicePlane->getDrawWidget();
		if (_slicePlane->getEnabled())
		{
			_sliceCentre = _slicePlane->getCenter();
			_sliceNormal = _slicePlane->getNormal();
		}
	}

	// The slice shows the field volume's cross-section only while the volume
	// is shown (Polyscope draws the section whether the mesh is enabled or
	// not), and inspecting a mesh resets the plane drawing, so that is
	// re-applied.
	void ViewerSession::SyncSliceInspection()
	{
		if (_slicePlane == nullptr)
		{
			return;
		}
		const bool wanted = _volumeView != nullptr && _volumeView->isEnabled();
		const bool inspecting = !_slicePlane->getVolumeMeshToInspect().empty();
		if (wanted != inspecting)
		{
			const bool drawPlane = _slicePlane->getDrawPlane();
			const bool drawWidget = _slicePlane->getDrawWidget();
			_slicePlane->setVolumeMeshToInspect(wanted ? kVolumeViewName : "");
			_slicePlane->setDrawPlane(drawPlane);
			_slicePlane->setDrawWidget(drawWidget);
		}
	}

	// Only the field volume is cut unless the user asks for the whole scene.
	void ViewerSession::ApplySliceIgnores()
	{
		if (_slicePlane == nullptr)
		{
			return;
		}
		polyscope::Structure* structures[] = { _sceneMesh, _edges, _wedgeEdges, _diffractingEdges, _transmitter,
		                                       _transmitterAxes, _receivers, _paths, _planeView };
		for (polyscope::Structure* s : structures)
		{
			if (s != nullptr)
			{
				s->setIgnoreSlicePlane(_slicePlane->name, !_sliceCutsScene);
			}
		}
	}

	void ViewerSession::SetSlicePose()
	{
		if (_slicePlane == nullptr)
		{
			return;
		}
		// Never hand Polyscope a non-finite pose (a disabled plane reports an
		// infinite centre).
		if (!std::isfinite(_sliceCentre.x) || !std::isfinite(_sliceCentre.y) || !std::isfinite(_sliceCentre.z))
		{
			_sliceCentre = ToGlm(0.5 * (_boundMin + _boundMax));
		}
		const float length = glm::length(_sliceNormal);
		if (!std::isfinite(length) || length < 1e-6f)
		{
			_sliceNormal = glm::vec3(0.0f, 0.0f, 1.0f);
		}
		_slicePlane->setPose(_sliceCentre, _sliceNormal);
	}

	// A three-quarter view from the front-right and above (z up) that frames
	// everything shown; Polyscope's own home view is a flat front elevation.
	void ViewerSession::FrameCamera()
	{
		const Vec3 centre = 0.5 * (_boundMin + _boundMax);
		const Vec3 eye = centre + _boundDiagonal * Vec3(0.9, -1.1, 0.7);
		polyscope::view::lookAt(ToGlm(eye), ToGlm(centre));
	}

	// The quantity selector is shared by the receivers, the 2D and the 3D view.
	void ViewerSession::SelectQuantity(int quantityIndex)
	{
		_selectedQuantity = std::max(0, std::min(kQuantityCount - 1, quantityIndex));
		const std::string name = FieldQuantityName(static_cast<FieldQuantity>(_selectedQuantity));
		polyscope::Structure* structures[] = { _receivers, _planeView, _volumeView };
		for (polyscope::Structure* s : structures)
		{
			if (s == nullptr)
			{
				continue;
			}
			polyscope::Quantity* quantity = s->getQuantity(name);
			if (quantity != nullptr)
			{
				quantity->setEnabled(true);   // scalar quantities are exclusive: the previous one is disabled
			}
		}
	}

	// 7. Control panel (inside Polyscope's user window)
	void ViewerSession::BuildGui()
	{
		// Polyscope's own panels may have changed what the controls mirror.
		SyncSliceState();
		SyncSliceInspection();
		ImGui::PushItemWidth(150.0f);
		ImGui::TextUnformatted("SBRSolver viewer");
		ImGui::Separator();
		if (ImGui::CollapsingHeader("Run", ImGuiTreeNodeFlags_DefaultOpen))
		{
			GuiRunInfo();
		}
		if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
		{
			GuiOverlays();
		}
		if (ImGui::CollapsingHeader("Field quantity", ImGuiTreeNodeFlags_DefaultOpen))
		{
			std::array<const char*, kQuantityCount> names;
			for (int q = 0; q < kQuantityCount; ++q)
			{
				names[q] = FieldQuantityName(static_cast<FieldQuantity>(q));
			}
			int selected = _selectedQuantity;
			if (ImGui::Combo("Quantity", &selected, names.data(), kQuantityCount))
			{
				SelectQuantity(selected);
			}
			ImGui::TextDisabled("shown on the receivers, the 2D and the 3D view");
		}
		if (ImGui::CollapsingHeader("Scene"))
		{
			GuiScene();
		}
		if (ImGui::CollapsingHeader("Ray paths"))
		{
			GuiPaths();
		}
		if (ImGui::CollapsingHeader("Interpolation"))
		{
			GuiInterpolation();
		}
		if (ImGui::CollapsingHeader("2D field view"))
		{
			GuiPlaneView();
		}
		if (ImGui::CollapsingHeader("3D slice view"))
		{
			GuiVolumeView();
		}
		ImGui::PopItemWidth();
	}

	void ViewerSession::GuiRunInfo()
	{
		ImGui::TextWrapped("%s", _scene.title.empty() ? "(untitled run)" : _scene.title.c_str());
		if (!_scene.outputDir.empty())
		{
			ImGui::TextWrapped("Output: %s", _scene.outputDir.c_str());
		}
		ImGui::Text("Frequency: %s", FormatFrequency(_scene.frequencyHz).c_str());
		if (_scene.txPowerWatts > 0.0)
		{
			ImGui::Text("Tx power: %.4g W (%.1f dBm)", _scene.txPowerWatts, 10.0 * std::log10(_scene.txPowerWatts * 1e3));
		}
		else
		{
			ImGui::Text("Tx power: %.4g W", _scene.txPowerWatts);
		}
		ImGui::Text("Antenna: %s", _scene.antennaType.empty() ? "-" : _scene.antennaType.c_str());
		ImGui::Text("Tx position: (%.3f, %.3f, %.3f) m", _scene.txPosition.x(), _scene.txPosition.y(), _scene.txPosition.z());
		size_t wedges = 0;
		size_t diffracting = 0;
		for (const SceneEdge& e : _scene.edges)
		{
			wedges += e.wedge ? 1 : 0;
			diffracting += e.diffracting ? 1 : 0;
		}
		ImGui::Text("Mesh: %zu vertices, %zu faces, %zu surfaces", _scene.mesh.vertices.size(), _scene.mesh.faces.size(),
		            size_t(_scene.mesh.coplanarId.empty() ? 0 : *std::max_element(_scene.mesh.coplanarId.begin(), _scene.mesh.coplanarId.end()) + 1));
		ImGui::Text("Edges: %zu (%zu wedges, %zu diffracting)", _scene.edges.size(), wedges, diffracting);
		ImGui::Text("Receivers: %zu", _scene.receivers.size());
		ImGui::Text("Paths: %zu", _scene.paths.size());
	}

	void ViewerSession::GuiOverlays()
	{
		// Some labels repeat the section headers' names; the id scope keeps
		// the checkboxes apart from them.
		ImGui::PushID("overlays");
		OverlayCheckbox("Scene mesh", _sceneMesh, _showSceneMesh);
		OverlayCheckbox("Edges", _edges, _showEdges);
		OverlayCheckbox("Wedge edges", _wedgeEdges, _showWedgeEdges);
		OverlayCheckbox("Diffracting edges", _diffractingEdges, _showDiffractingEdges);
		OverlayCheckbox("Transmitter", _transmitter, _showTransmitter);
		OverlayCheckbox("Transmitter axes", _transmitterAxes, _showTransmitterAxes);
		OverlayCheckbox("Receivers", _receivers, _showReceivers);
		OverlayCheckbox("Ray paths", _paths, _showPaths);
		OverlayCheckbox("2D field view", _planeView, _showPlaneView);
		OverlayCheckbox("3D field volume", _volumeView, _showVolumeView);
		if (_slicePlane != nullptr)
		{
			_showSlicePlane = _slicePlane->getEnabled();
			if (ImGui::Checkbox("Slice plane", &_showSlicePlane))
			{
				_slicePlane->setEnabled(_showSlicePlane);
				if (_showSlicePlane)
				{
					SetSlicePose();   // back where it was, not where a disabled plane sits
				}
			}
		}
		ImGui::PopID();
	}

	void ViewerSession::GuiScene()
	{
		if (_sceneMesh != nullptr)
		{
			_sceneTransparency = _sceneMesh->getTransparency();   // may be changed in Polyscope's own panel
			if (ImGui::SliderFloat("Transparency", &_sceneTransparency, 0.05f, 1.0f))
			{
				_sceneMesh->setTransparency(_sceneTransparency);
			}
			if (ImGui::Combo("Hidden faces", &_faceCulling, kFaceCullingNames, 3))
			{
				RegisterSceneMesh();
			}
		}
		ImGui::TextUnformatted("Materials (id: name, faces)");
		std::vector<int> facesPerMaterial(_scene.mesh.materialNames.size(), 0);
		for (int id : _scene.mesh.materialId)
		{
			if (id >= 0 && id < static_cast<int>(facesPerMaterial.size()))
			{
				++facesPerMaterial[id];
			}
		}
		if (_scene.mesh.materialNames.empty())
		{
			ImGui::BulletText("(none)");
		}
		for (size_t id = 0; id < _scene.mesh.materialNames.size(); ++id)
		{
			ImGui::BulletText("%zu: %s (%d)", id, _scene.mesh.materialNames[id].c_str(), facesPerMaterial[id]);
		}
	}

	void ViewerSession::GuiPaths()
	{
		const int maxReceiver = static_cast<int>(_scene.receivers.size()) - 1;
		int filter = _pathReceiverFilter;
		if (ImGui::InputInt("Receiver filter (-1: all)", &filter))
		{
			_pathReceiverFilter = std::max(-1, std::min(maxReceiver, filter));
			RegisterPaths();
		}
		if (_pathReceiverFilter >= 0 && _pathReceiverFilter < static_cast<int>(_scene.receivers.size()))
		{
			const Vec3& r = _scene.receivers[_pathReceiverFilter];
			ImGui::Text("Receiver %d at (%.3f, %.3f, %.3f) m", _pathReceiverFilter, r.x(), r.y(), r.z());
			ImGui::Text("|E| = %.4g V/m (%.1f dBuV/m)", _receiverFields[_pathReceiverFilter].norm(),
			            EvaluateQuantity(_receiverFields[_pathReceiverFilter], FieldQuantity::MagnitudeDb, _scene.frequencyHz));
		}
		ImGui::Text("Shown: %zu of %zu paths", _pathsShown, _scene.paths.size());
	}

	void ViewerSession::GuiInterpolation()
	{
		int method = static_cast<int>(_interpolation.method);
		const char* methods[] = { MethodName(InterpolationMethod::Nearest), MethodName(InterpolationMethod::InverseDistance) };
		if (ImGui::Combo("Method", &method, methods, 2))
		{
			_interpolation.method = static_cast<InterpolationMethod>(method);
		}
		if (ImGui::InputInt("Neighbours", &_interpolation.neighbours))
		{
			_interpolation.neighbours = std::max(1, _interpolation.neighbours);
		}
		float power = static_cast<float>(_interpolation.power);
		if (ImGui::SliderFloat("Power", &power, 0.5f, 6.0f))
		{
			_interpolation.power = power;
		}
		float radius = static_cast<float>(_interpolation.radiusFactor);
		if (ImGui::SliderFloat("Radius (x sample spacing, 0: none)", &radius, 0.0f, 10.0f))
		{
			_interpolation.radiusFactor = radius;
		}
		ImGui::Text("Sample spacing: %.3f m", SampleSpacing(_scene.receivers));
		if (ImGui::Button("Recompute both views"))
		{
			RegisterPlaneView();
			RegisterVolumeView();
		}
	}

	void ViewerSession::GuiPlaneView()
	{
		if (ImGui::InputInt("Resolution (2D)", &_planeResolution))
		{
			_planeResolution = std::max(2, std::min(kMaxPlaneResolution, _planeResolution));
		}
		if (ImGui::Button("Recompute 2D view"))
		{
			RegisterPlaneView();
		}
		ImGui::TextWrapped("%s", _planeStatus.c_str());
		if (_planeFit.ok)
		{
			ImGui::Text("Plane normal: (%.3f, %.3f, %.3f)", _planeFit.normal.x(), _planeFit.normal.y(), _planeFit.normal.z());
		}
	}

	void ViewerSession::GuiVolumeView()
	{
		if (ImGui::InputInt("Resolution (3D)", &_volumeResolution))
		{
			_volumeResolution = std::max(2, std::min(kMaxVolumeResolution, _volumeResolution));
		}
		if (ImGui::Button("Recompute 3D view"))
		{
			RegisterVolumeView();
		}
		ImGui::TextWrapped("%s", _volumeStatus.c_str());
		if (_volumeView != nullptr)
		{
			_volumeTransparency = _volumeView->getTransparency();   // may be changed in Polyscope's own panel
			if (ImGui::SliderFloat("Volume transparency", &_volumeTransparency, 0.05f, 1.0f))
			{
				_volumeView->setTransparency(_volumeTransparency);
			}
		}
		if (_slicePlane == nullptr)
		{
			return;
		}
		ImGui::SeparatorText("Slice plane");
		if (ImGui::Checkbox("Draw plane", &_sliceDrawPlane))
		{
			_slicePlane->setDrawPlane(_sliceDrawPlane);
		}
		ImGui::SameLine();
		if (ImGui::Checkbox("Gizmo", &_sliceDrawWidget))
		{
			_slicePlane->setDrawWidget(_sliceDrawWidget);
		}
		if (ImGui::Checkbox("Cut the scene too (not only the volume)", &_sliceCutsScene))
		{
			ApplySliceIgnores();
		}
		// The pose mirrors follow the gizmo while the plane is enabled
		// (SyncSliceState); a disabled plane cannot be moved.
		const bool sliceEnabled = _slicePlane->getEnabled();
		if (!sliceEnabled)
		{
			ImGui::TextDisabled("Enable the slice plane (Overlays) to move it");
		}
		ImGui::BeginDisabled(!sliceEnabled);
		ImGui::TextUnformatted("Normal:");
		const char* axisLabels[3] = { "X", "Y", "Z" };
		for (int axis = 0; axis < 3; ++axis)
		{
			ImGui::SameLine();
			if (ImGui::Button(axisLabels[axis]))
			{
				_sliceAxis = axis;
				const float sign = _sliceNormal[axis] < 0.0f ? -1.0f : 1.0f;
				_sliceNormal = glm::vec3(0.0f);
				_sliceNormal[axis] = sign;
				SetSlicePose();
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("Flip"))
		{
			_sliceNormal = -_sliceNormal;
			SetSlicePose();
		}
		bool moved = false;
		const Vec3 margin = 0.05 * (_boundMax - _boundMin) + Vec3::Constant(0.01);
		const Vec3 low = _boundMin - margin;
		const Vec3 high = _boundMax + margin;
		const char* positionLabels[3] = { "Centre x [m]", "Centre y [m]", "Centre z [m]" };
		for (int axis = 0; axis < 3; ++axis)
		{
			moved |= ImGui::SliderFloat(positionLabels[axis], &_sliceCentre[axis], static_cast<float>(low[axis]), static_cast<float>(high[axis]));
		}
		if (moved)
		{
			SetSlicePose();
		}
		ImGui::EndDisabled();
		ImGui::Text("Normal (%.2f, %.2f, %.2f); the side it points away from is cut", _sliceNormal.x, _sliceNormal.y, _sliceNormal.z);
		ImGui::TextDisabled("Volume bounds: (%.2f, %.2f, %.2f) - (%.2f, %.2f, %.2f)", _volumeGrid.boundMin.x(), _volumeGrid.boundMin.y(),
		                    _volumeGrid.boundMin.z(), _volumeGrid.boundMax.x(), _volumeGrid.boundMax.y(), _volumeGrid.boundMax.z());
	}

	void ConfigurePolyscope(bool interactive)
	{
		polyscope::options::programName = "SBRSolver";
		polyscope::options::usePrefsFile = false;
		polyscope::options::verbosity = 0;
		polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;   // the scenes have their own floors
		polyscope::options::automaticallyComputeSceneExtents = true;
		// Without a display the automatic backend choice may fall back to
		// EGL (headless): fine for a screenshot or a fixed number of frames,
		// but an interactive viewer would then run forever with no window, so
		// it fails instead. Errors in the non-interactive modes throw (and end
		// the run) rather than wait for a click in a popup nobody can see.
		polyscope::options::allowHeadlessBackends = !interactive;
		polyscope::options::errorsThrowExceptions = !interactive;
		polyscope::options::displayMessagePopups = interactive;
	}
} // namespace

int RunViewer(const SceneData& scene, const ViewerOptions& options)
{
	const bool screenshot = !options.screenshotFile.empty();
	const bool interactive = !screenshot && options.frames < 0;
	bool initialised = false;
	try
	{
		ConfigurePolyscope(interactive);
		polyscope::init(options.backend);
		initialised = true;
		if (interactive && polyscope::isHeadless())
		{
			// An explicitly requested headless backend (EGL, mock) has no
			// window to close.
			std::cerr << "Viewer: no display is available (backend " << polyscope::state::backend
			          << " opens no window); use --screenshot FILE or --viewer-frames N to run without one" << std::endl;
			polyscope::shutdown();
			return 1;
		}
		polyscope::view::setUpDir(polyscope::UpDir::ZUp);   // z is up in the scenes

		ViewerSession session(scene, options);
		session.RegisterAll();
		polyscope::state::userCallback = [&session]()
		{
			SBR_FRAME_MARK();   // one profiler frame per rendered frame
			session.BuildGui();
		};
		polyscope::view::resetCameraToHomeView();
		session.FrameCamera();

		if (screenshot)
		{
			polyscope::show(1);   // one frame, so the window and the UI exist
			if (polyscope::state::backend == kMockBackend)
			{
				polyscope::screenshotToBuffer();   // the render path, without a file
				std::cout << "Viewer: mock backend, screenshot " << options.screenshotFile << " not written" << std::endl;
			}
			else
			{
				// Success means a fresh, non-empty file: a stale file is removed
				// first (Polyscope does not report a failed write), and a
				// directory is never touched.
				auto fail = [](const std::string& message)
				{
					std::cerr << "Viewer: " << message << std::endl;
					polyscope::state::userCallback = nullptr;
					polyscope::shutdown();
					return 1;
				};
				const std::filesystem::path file(options.screenshotFile);
				std::error_code ec;
				if (std::filesystem::is_directory(file, ec))
				{
					return fail(options.screenshotFile + " is a directory");
				}
				std::filesystem::remove(file, ec);
				if (std::filesystem::exists(file, ec))
				{
					return fail("cannot replace the existing file " + options.screenshotFile);
				}
				polyscope::ScreenshotOptions screenshotOptions;
				screenshotOptions.transparentBackground = false;
				polyscope::screenshot(options.screenshotFile, screenshotOptions);
				ec.clear();
				const std::uintmax_t size = std::filesystem::file_size(file, ec);
				if (ec || size == 0 || !std::filesystem::is_regular_file(file, ec))
				{
					return fail("could not write the screenshot " + options.screenshotFile);
				}
				std::cout << "Viewer: screenshot written to " << options.screenshotFile << std::endl;
			}
		}
		else if (options.frames >= 0)
		{
			polyscope::show(static_cast<size_t>(options.frames));
		}
		else
		{
			polyscope::show();
		}

		polyscope::state::userCallback = nullptr;
		polyscope::shutdown();   // so a second RunViewer call starts afresh
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Viewer: " << e.what() << std::endl;
		if (interactive && !initialised)
		{
			std::cerr << "Viewer: a window could not be opened; without a display use --screenshot FILE or --viewer-frames N"
			          << std::endl;
		}
		try
		{
			polyscope::state::userCallback = nullptr;
			if (polyscope::isInitialized())
			{
				polyscope::shutdown(true);
			}
		}
		catch (const std::exception&)
		{
		}
		return 1;
	}
}

} // namespace viewer
} // namespace sbr
