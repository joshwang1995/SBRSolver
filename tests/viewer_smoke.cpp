// Smoke test of the built-in viewer: a synthetic run (a box room with a
// partition, a transmitter, a plane of receivers with a synthetic field and a
// few paths) goes through every overlay of RunViewer. With Polyscope's mock
// backend nothing is rendered, so it runs without a display (ctest
// viewer_smoke_mock); with a display or xvfb-run, --screenshot renders it.
//
//   sbr_viewer_smoke [--backend NAME] [--frames N] [--runs N] [--screenshot FILE]
//                    [--receiver-layers N] [--plane-resolution N] [--volume-resolution N]
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "viewer/SceneData.h"
#include "viewer/Viewer.h"

namespace
{
	using namespace sbr::viewer;

	const double kPi = 3.14159265358979323846;

	void AddQuad(SceneMesh& mesh, int a, int b, int c, int d, int surfaceId, int materialId)
	{
		mesh.faces.push_back({ a, b, c });
		mesh.faces.push_back({ a, c, d });
		for (int i = 0; i < 2; ++i)
		{
			mesh.coplanarId.push_back(surfaceId);
			mesh.materialId.push_back(materialId);
		}
	}

	// A 6 x 5 x 3 m room (normals pointing inward, as the solver's indoor
	// scenes do) with a 1 x 1 m partition standing on the floor.
	SceneMesh BuildRoom()
	{
		SceneMesh mesh;
		mesh.materialNames = { "Concrete", "Plasterboard" };
		const double lx = 6.0, ly = 5.0, lz = 3.0;
		mesh.vertices = {
			Vec3(0, 0, 0), Vec3(lx, 0, 0), Vec3(lx, ly, 0), Vec3(0, ly, 0),
			Vec3(0, 0, lz), Vec3(lx, 0, lz), Vec3(lx, ly, lz), Vec3(0, ly, lz)
		};
		AddQuad(mesh, 0, 1, 2, 3, 0, 0);   // floor (normal +z)
		AddQuad(mesh, 4, 7, 6, 5, 1, 0);   // ceiling (normal -z)
		AddQuad(mesh, 0, 4, 5, 1, 2, 0);   // wall y = 0 (normal +y)
		AddQuad(mesh, 1, 5, 6, 2, 3, 0);   // wall x = lx (normal -x)
		AddQuad(mesh, 2, 6, 7, 3, 4, 0);   // wall y = ly (normal -y)
		AddQuad(mesh, 3, 7, 4, 0, 5, 0);   // wall x = 0 (normal +x)

		// Partition in the plane x = 3, y in [2, 3], z in [0, 1]
		const int base = static_cast<int>(mesh.vertices.size());
		mesh.vertices.push_back(Vec3(3, 2, 0));
		mesh.vertices.push_back(Vec3(3, 3, 0));
		mesh.vertices.push_back(Vec3(3, 3, 1));
		mesh.vertices.push_back(Vec3(3, 2, 1));
		AddQuad(mesh, base, base + 1, base + 2, base + 3, 6, 1);
		return mesh;
	}

	// A spherical wave from the transmitter, mostly z-polarised.
	Vec3c SyntheticField(const Vec3& tx, const Vec3& rx, double k)
	{
		const double d = std::max((rx - tx).norm(), 0.1);
		const std::complex<double> phasor = std::polar(1.0 / d, -k * d);
		return Vec3c(0.2 * phasor, 0.3 * phasor * std::polar(1.0, 0.5), 1.0 * phasor);
	}

	// `layers` planes of 10 x 8 receivers, 0.5 m apart around z = 1.5 m; one
	// layer is the flat layout (no field volume), three give a 3D one.
	SceneData BuildScene(int layers)
	{
		SceneData scene;
		scene.title = "viewer_smoke (synthetic)";
		scene.outputDir = "(none)";
		scene.mesh = BuildRoom();
		scene.edges = BuildSceneEdges(scene.mesh);
		for (SceneEdge& e : scene.edges)
		{
			if (e.wedge)   // exercise the diffracting-edge overlay on one wedge
			{
				e.diffracting = true;
				break;
			}
		}

		scene.txPosition = Vec3(1.0, 1.0, 2.5);
		scene.txCoordSys = Mat3::Identity();
		scene.antennaType = "Hertzian dipole (synthetic)";
		scene.frequencyHz = 2.4e9;
		scene.txPowerWatts = 0.1;
		const double k = 2.0 * kPi * scene.frequencyHz / 299792458.0;

		// 10 x 8 receivers per layer, the single layer at z = 1.5 m
		for (int layer = 0; layer < std::max(1, layers); ++layer)
		{
			const double z = 1.5 + 0.5 * (layer - 0.5 * (std::max(1, layers) - 1));
			for (int j = 0; j < 8; ++j)
			{
				for (int i = 0; i < 10; ++i)
				{
					const Vec3 rx(0.5 + 5.0 * i / 9.0, 0.5 + 4.0 * j / 7.0, z);
					scene.receivers.push_back(rx);
					// a few receivers without any field, as behind an opaque wall
					const bool shadowed = (i == 9 && j >= 6);
					scene.fields.push_back(shadowed ? Vec3c(Vec3c::Zero()) : SyntheticField(scene.txPosition, rx, k));
				}
			}
		}

		// A few paths: direct, one floor reflection and one transmission
		// through the partition, to the first receivers of each row.
		for (int receiver = 0; receiver < static_cast<int>(scene.receivers.size()); receiver += 10)
		{
			const Vec3& rx = scene.receivers[receiver];
			ScenePath direct;
			direct.receiverId = receiver;
			direct.launchId = 0;
			direct.points = { scene.txPosition, rx };
			direct.field = SyntheticField(scene.txPosition, rx, k);
			direct.length = (rx - scene.txPosition).norm();
			scene.paths.push_back(direct);

			ScenePath reflected;
			reflected.receiverId = receiver;
			reflected.launchId = 1;
			// image of the transmitter in the floor; the image ray crosses z = 0 at the bounce
			const Vec3 image(scene.txPosition.x(), scene.txPosition.y(), -scene.txPosition.z());
			const double t = scene.txPosition.z() / (scene.txPosition.z() + rx.z());
			const Vec3 bounce = image + t * (rx - image);
			reflected.points = { scene.txPosition, Vec3(bounce.x(), bounce.y(), 0.0), rx };
			reflected.interactions = { InteractionType::Reflection };
			reflected.field = 0.5 * SyntheticField(image, rx, k);
			reflected.length = (rx - image).norm();
			scene.paths.push_back(reflected);

			ScenePath transmitted;
			transmitted.receiverId = receiver;
			transmitted.launchId = 2;
			transmitted.points = { scene.txPosition, Vec3(3.0, 2.5, 0.5), rx };
			transmitted.interactions = { InteractionType::Transmission };
			transmitted.field = 0.1 * SyntheticField(scene.txPosition, rx, k);
			transmitted.length = (Vec3(3.0, 2.5, 0.5) - scene.txPosition).norm() + (rx - Vec3(3.0, 2.5, 0.5)).norm();
			scene.paths.push_back(transmitted);
		}
		return scene;
	}

	int Usage(const char* program)
	{
		std::cerr << "usage: " << program
		          << " [--backend NAME] [--frames N] [--runs N] [--screenshot FILE] [--receiver-layers N]"
		             " [--plane-resolution N] [--volume-resolution N]\n";
		return 2;
	}
} // namespace

int main(int argc, char** argv)
{
	ViewerOptions options;
	options.backend = "openGL_mock";
	options.frames = 2;
	int runs = 1;
	int layers = 1;
	bool backendGiven = false;

	for (int i = 1; i < argc; ++i)
	{
		const std::string arg = argv[i];
		const bool hasValue = i + 1 < argc;
		if (arg == "--backend" && hasValue)
		{
			options.backend = argv[++i];
			backendGiven = true;
		}
		else if (arg == "--frames" && hasValue)
		{
			options.frames = std::atol(argv[++i]);
		}
		else if (arg == "--runs" && hasValue)
		{
			runs = std::atoi(argv[++i]);
		}
		else if (arg == "--receiver-layers" && hasValue)
		{
			layers = std::atoi(argv[++i]);
		}
		else if (arg == "--screenshot" && hasValue)
		{
			options.screenshotFile = argv[++i];
		}
		else if (arg == "--plane-resolution" && hasValue)
		{
			options.planeResolution = std::atoi(argv[++i]);
		}
		else if (arg == "--volume-resolution" && hasValue)
		{
			options.volumeResolution = std::atoi(argv[++i]);
		}
		else
		{
			return Usage(argv[0]);
		}
	}
	if (!options.screenshotFile.empty() && !backendGiven)
	{
		options.backend = "";   // a real backend: the window, or EGL without a display
	}

	const SceneData scene = BuildScene(layers);
	std::cout << "viewer_smoke: " << scene.mesh.faces.size() << " faces, " << scene.edges.size() << " edges, "
	          << scene.receivers.size() << " receivers in " << std::max(1, layers) << " layer(s), " << scene.paths.size()
	          << " paths; backend '"
	          << options.backend << "', frames " << options.frames << ", runs " << runs << std::endl;

	for (int run = 0; run < std::max(1, runs); ++run)
	{
		const int result = RunViewer(scene, options);
		if (result != 0)
		{
			std::cerr << "viewer_smoke: RunViewer returned " << result << " in run " << run + 1 << std::endl;
			return result;
		}
	}
	std::cout << "viewer_smoke: ok" << std::endl;
	return 0;
}
