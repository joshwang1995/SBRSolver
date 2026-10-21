// Unit tests of the viewer's scene data (viewer/SceneData.h) and field
// interpolation (viewer/FieldGrid.h), with a self-contained harness: CHECK counts failures and main returns non-zero when
// any test failed.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "Triangle.h"
#include "common/Constants.h"
#include "viewer/FieldGrid.h"
#include "viewer/SceneData.h"

using namespace sbr::viewer;

namespace
{
	int g_checks = 0;
	int g_failures = 0;
	const char* g_test = "";

	void Fail(const char* file, int line, const char* expression)
	{
		g_failures++;
		std::printf("  FAILED %s:%d (%s): %s\n", file, line, g_test, expression);
	}

#define CHECK(condition)                                               \
	do                                                                 \
	{                                                                  \
		g_checks++;                                                    \
		if (!(condition)) Fail(__FILE__, __LINE__, #condition);        \
	} while (0)

#define CHECK_NEAR(a, b, tolerance)                                    \
	do                                                                 \
	{                                                                  \
		g_checks++;                                                    \
		const double _a = (a);                                         \
		const double _b = (b);                                         \
		if (!(std::abs(_a - _b) <= (tolerance)))                       \
		{                                                              \
			Fail(__FILE__, __LINE__, #a " ~= " #b);                    \
			std::printf("    %.17g vs %.17g\n", _a, _b);               \
		}                                                              \
	} while (0)

	void Run(const char* name, void (*test)())
	{
		g_test = name;
		const int before = g_failures;
		std::printf("[ RUN  ] %s\n", name);
		test();
		std::printf("[ %s ] %s\n", g_failures == before ? " OK " : "FAIL", name);
	}

	typedef std::complex<double> Complex;
	typedef std::chrono::steady_clock Clock;

	double Seconds(Clock::time_point start)
	{
		return std::chrono::duration<double>(Clock::now() - start).count();
	}

	// Regular grid x = 0..nx-1, y = 0..ny-1 at height z (1 m spacing).
	std::vector<Vec3> PlaneSamples(int nx, int ny, double z)
	{
		std::vector<Vec3> samples;
		for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++) samples.push_back(Vec3(double(i), double(j), z));
		}
		return samples;
	}

	// Regular box x = 0..nx-1, y = 0..ny-1, z = 0..nz-1 (1 m spacing).
	std::vector<Vec3> BoxSamples(int nx, int ny, int nz)
	{
		std::vector<Vec3> samples;
		for (int k = 0; k < nz; k++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int i = 0; i < nx; i++) samples.push_back(Vec3(double(i), double(j), double(k)));
			}
		}
		return samples;
	}

	// Receiver line x = 1, z = 2, y = 0..n-1.
	std::vector<Vec3> LineSamples(int n)
	{
		std::vector<Vec3> samples;
		for (int j = 0; j < n; j++) samples.push_back(Vec3(1.0, double(j), 2.0));
		return samples;
	}

	// A complex field linear in x and y, different in every component.
	Vec3c LinearField(const Vec3& p)
	{
		return Vec3c(Complex(1.0 + 0.5 * p.x(), -2.0 * p.y()), Complex(p.x() - p.y(), 3.0),
		             Complex(0.25, p.x() + 2.0 * p.y()));
	}

	std::vector<Vec3c> LinearFieldAt(const std::vector<Vec3>& samples)
	{
		std::vector<Vec3c> values;
		for (const Vec3& p : samples) values.push_back(LinearField(p));
		return values;
	}

	std::vector<Vec3c> RandomField(size_t n, unsigned seed)
	{
		std::mt19937 rng(seed);
		std::uniform_real_distribution<double> uniform(-1.0, 1.0);
		std::vector<Vec3c> values(n);
		for (Vec3c& v : values)
		{
			for (int c = 0; c < 3; c++) v[c] = Complex(uniform(rng), uniform(rng));
		}
		return values;
	}

	std::vector<Vec3> RandomPoints(size_t n, double size, bool flat, unsigned seed)
	{
		std::mt19937 rng(seed);
		std::uniform_real_distribution<double> uniform(0.0, size);
		std::vector<Vec3> points(n);
		for (Vec3& p : points) p = Vec3(uniform(rng), uniform(rng), flat ? 0.0 : uniform(rng));
		return points;
	}

	double RelativeError(const Vec3c& actual, const Vec3c& expected)
	{
		return (actual - expected).norm() / std::max(1e-300, expected.norm());
	}

	size_t NearestIndex(const std::vector<Vec3>& samples, const Vec3& q)
	{
		size_t best = 0;
		for (size_t i = 1; i < samples.size(); i++)
		{
			if ((samples[i] - q).squaredNorm() < (samples[best] - q).squaredNorm()) best = i;
		}
		return best;
	}


	// ---------------------------------------------------------------------
	// SceneData

	// A floor square (z = 0) and a wall square (y = 0) sharing the edge from
	// (0,0,0) to (1,0,0): 4 triangles, 6 distinct corners, one 90 degree wedge.
	struct FloorAndWall
	{
		std::vector<Triangle> triangles;
		std::vector<Triangle*> pointers;
		FloorAndWall()
		{
			triangles.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), 0));
			triangles.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 1, 0), Vec3(0, 1, 0), 1));
			triangles.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 1), Vec3(1, 0, 0), 2));
			triangles.push_back(Triangle(Vec3(0, 0, 0), Vec3(0, 0, 1), Vec3(1, 0, 1), 3));
			for (size_t i = 0; i < triangles.size(); i++)
			{
				triangles[i].coplanarId = i < 2 ? 0 : 1;
				triangles[i].materialId = i < 2 ? 0 : 1;
				pointers.push_back(&triangles[i]);
			}
		}
	};

	void TestBuildSceneMesh()
	{
		FloorAndWall scene;
		const std::vector<std::string> names = { "concrete", "glass" };
		const SceneMesh mesh = BuildSceneMesh(scene.pointers, names);
		CHECK(mesh.vertices.size() == 6);   // welded by exact coordinates
		CHECK(mesh.faces.size() == 4);
		CHECK(mesh.coplanarId.size() == 4);
		CHECK(mesh.materialId.size() == 4);
		CHECK(mesh.materialNames == names);
		for (size_t f = 0; f < mesh.faces.size(); f++)
		{
			const Triangle& t = scene.triangles[f];
			// Face f is triangle f with its corners in the triangle's order.
			CHECK(mesh.vertices[mesh.faces[f][0]] == t.v1);
			CHECK(mesh.vertices[mesh.faces[f][1]] == t.v2);
			CHECK(mesh.vertices[mesh.faces[f][2]] == t.v3);
			CHECK(mesh.coplanarId[f] == t.coplanarId);
			CHECK(mesh.materialId[f] == t.materialId);
		}
		// The shared corners are the same vertex in both squares.
		CHECK(mesh.faces[0][0] == mesh.faces[2][0]);   // (0,0,0)
		CHECK(mesh.faces[0][1] == mesh.faces[2][2]);   // (1,0,0)
		// The first corner seen gets index 0, later corners follow in order.
		CHECK(mesh.faces[0][0] == 0);
		CHECK(mesh.faces[0][1] == 1);
		CHECK(mesh.faces[0][2] == 2);
		CHECK(mesh.faces[1][2] == 3);

		const SceneMesh empty = BuildSceneMesh(std::vector<Triangle*>(), names);
		CHECK(empty.vertices.empty());
		CHECK(empty.faces.empty());
		CHECK(empty.materialNames == names);
	}

	void TestBuildSceneEdges()
	{
		FloorAndWall scene;
		const SceneMesh mesh = BuildSceneMesh(scene.pointers, std::vector<std::string>());
		const std::vector<SceneEdge> edges = BuildSceneEdges(mesh);
		CHECK(edges.size() == 9);   // 5 per square, the bottom edge of the wall shared
		int boundary = 0;
		int flat = 0;
		int wedges = 0;
		for (const SceneEdge& e : edges)
		{
			CHECK(e.v0 < e.v1);
			CHECK(e.v0 >= 0 && e.v1 < int(mesh.vertices.size()));
			CHECK(!e.diffracting);   // nothing marks edges yet
			if (e.faceCount == 1)
			{
				boundary++;
				CHECK(e.face0 >= 0 && e.face1 == -1);
				CHECK(!e.wedge);
			}
			else
			{
				CHECK(e.faceCount == 2);
				CHECK(e.face0 >= 0 && e.face1 >= 0 && e.face0 != e.face1);
				if (e.wedge)
				{
					wedges++;
					CHECK_NEAR(e.dihedralAngleDeg, 90.0, 1e-9);
					CHECK(mesh.vertices[e.v0] == Vec3(0, 0, 0));
					CHECK(mesh.vertices[e.v1] == Vec3(1, 0, 0));
					CHECK((e.face0 < 2) != (e.face1 < 2));   // one floor face, one wall face
				}
				else
				{
					flat++;
					// A square's diagonal: both faces lie in one plane.
					CHECK(std::abs(e.dihedralAngleDeg) < 1e-9 || std::abs(e.dihedralAngleDeg - 180.0) < 1e-9);
					CHECK((e.face0 < 2) == (e.face1 < 2));
				}
			}
		}
		CHECK(boundary == 6);
		CHECK(flat == 2);
		CHECK(wedges == 1);

		// A fold of 0.005 degrees is coplanar at the default tolerance of
		// 0.01 degrees and a wedge below it.
		const double angle = 0.005 * PI / 180.0;
		std::vector<Triangle> fold;
		fold.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), 0));
		fold.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, -std::cos(angle), -std::sin(angle)), Vec3(1, 0, 0), 1));   // same winding sense as the first
		std::vector<Triangle*> foldPointers = { &fold[0], &fold[1] };
		const SceneMesh foldMesh = BuildSceneMesh(foldPointers, std::vector<std::string>());
		int foldWedges = 0;
		for (const SceneEdge& e : BuildSceneEdges(foldMesh))
		{
			if (e.faceCount == 2) CHECK_NEAR(e.dihedralAngleDeg, 0.005, 1e-6);
			foldWedges += e.wedge ? 1 : 0;
		}
		CHECK(foldWedges == 0);
		foldWedges = 0;
		for (const SceneEdge& e : BuildSceneEdges(foldMesh, 0.001)) foldWedges += e.wedge ? 1 : 0;
		CHECK(foldWedges == 1);

		// Three faces on one edge (a T junction) are non-manifold, not a wedge.
		std::vector<Triangle> junction = fold;
		junction.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 0, 1), 2));
		std::vector<Triangle*> junctionPointers = { &junction[0], &junction[1], &junction[2] };
		const SceneMesh junctionMesh = BuildSceneMesh(junctionPointers, std::vector<std::string>());
		int nonManifold = 0;
		for (const SceneEdge& e : BuildSceneEdges(junctionMesh, 0.001))
		{
			if (e.faceCount == 3)
			{
				nonManifold++;
				CHECK(!e.wedge);
			}
		}
		CHECK(nonManifold == 1);

		CHECK(BuildSceneEdges(SceneMesh()).empty());

		// A zero-area face has no normal, so its edges are never wedges.
		std::vector<Triangle> sliver;
		sliver.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), 0));
		sliver.push_back(Triangle(Vec3(0, 0, 0), Vec3(0.5, 0, 0), Vec3(1, 0, 0), 1));   // collinear corners
		std::vector<Triangle*> sliverPointers = { &sliver[0], &sliver[1] };
		const SceneMesh sliverMesh = BuildSceneMesh(sliverPointers, std::vector<std::string>());
		for (const SceneEdge& e : BuildSceneEdges(sliverMesh))
		{
			CHECK(!e.wedge);
			CHECK(std::isfinite(e.dihedralAngleDeg));
		}
	}

	void TestFieldQuantities()
	{
		const double frequency = 1.0e9;
		const double lambda = SPEED_OF_LIGHT / frequency;
		const Vec3c field(Complex(3.0, 4.0), Complex(0.0, 0.0), Complex(-1.0, 0.0));
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::Magnitude, frequency), std::sqrt(26.0), 1e-12);
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::MagnitudeDb, frequency), 20.0 * std::log10(std::sqrt(26.0) * 1e6), 1e-9);
		CHECK_NEAR(ReceivedPowerWatts(field, frequency), 26.0 * lambda * lambda / (4.0 * PI * ETA0), 1e-20);
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::PowerDbm, frequency),
		           10.0 * std::log10(ReceivedPowerWatts(field, frequency) * 1e3), 1e-9);
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::PhaseX, frequency), std::atan2(4.0, 3.0) * 180.0 / PI, 1e-9);
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::PhaseY, frequency), 0.0, 0.0);   // zero component
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::PhaseZ, frequency), 180.0, 1e-9);
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::MagnitudeX, frequency), 5.0, 1e-12);
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::MagnitudeY, frequency), 0.0, 0.0);
		CHECK_NEAR(EvaluateQuantity(field, FieldQuantity::MagnitudeZ, frequency), 1.0, 1e-12);
		CHECK_NEAR(ReceivedPowerWatts(field, 0.0), 0.0, 0.0);

		// A zero field stays finite: -300 dB and zero phase.
		const Vec3c zero = Vec3c::Zero();
		for (int q = 0; q < int(FieldQuantity::Count); q++)
		{
			const FieldQuantity quantity = FieldQuantity(q);
			const double value = EvaluateQuantity(zero, quantity, frequency);
			CHECK(std::isfinite(value));
			if (quantity == FieldQuantity::MagnitudeDb || quantity == FieldQuantity::PowerDbm)
			{
				CHECK_NEAR(value, -300.0, 0.0);
			}
			else
			{
				CHECK_NEAR(value, 0.0, 0.0);
			}
			CHECK(FieldQuantityName(quantity) != nullptr && std::string(FieldQuantityName(quantity)) != "?");
			CHECK(FieldQuantityIsPhase(quantity)
			      == (quantity == FieldQuantity::PhaseX || quantity == FieldQuantity::PhaseY || quantity == FieldQuantity::PhaseZ));
		}

		const std::vector<Vec3c> fields = { field, zero, 2.0 * field };
		const std::vector<double> magnitudes = EvaluateQuantity(fields, FieldQuantity::Magnitude, frequency);
		CHECK(magnitudes.size() == 3);
		CHECK_NEAR(magnitudes[0], std::sqrt(26.0), 1e-12);
		CHECK_NEAR(magnitudes[1], 0.0, 0.0);
		CHECK_NEAR(magnitudes[2], 2.0 * std::sqrt(26.0), 1e-12);
	}

	// ---------------------------------------------------------------------
	// FieldGrid

	void TestSampleSpacing()
	{
		CHECK_NEAR(SampleSpacing(PlaneSamples(11, 6, 1.5)), 1.0, 1e-12);
		CHECK_NEAR(SampleSpacing(BoxSamples(5, 4, 3)), 1.0, 1e-12);

		std::vector<Vec3> scaled = PlaneSamples(11, 6, 0.0);
		for (Vec3& p : scaled) p *= 0.25;
		CHECK_NEAR(SampleSpacing(scaled), 0.25, 1e-12);

		CHECK_NEAR(SampleSpacing({Vec3(0, 0, 0), Vec3(3, 4, 0)}), 5.0, 1e-12);
		CHECK_NEAR(SampleSpacing({Vec3(1, 2, 3)}), 0.0, 0.0);
		CHECK_NEAR(SampleSpacing({}), 0.0, 0.0);

		// Irregular: the median is robust to one outlier and one close pair.
		std::vector<Vec3> line = LineSamples(9);
		line.push_back(Vec3(1.0, 100.0, 2.0));
		line.push_back(Vec3(1.0, 0.001, 2.0));
		CHECK_NEAR(SampleSpacing(line), 1.0, 1e-12);
	}

	void TestFitPlane()
	{
		// Horizontal plane, longest extent along x.
		PlaneFit fit = FitPlane(PlaneSamples(11, 6, 1.5));
		CHECK(fit.ok);
		CHECK_NEAR(std::abs(fit.normal.z()), 1.0, 1e-12);
		CHECK_NEAR(std::abs(fit.u.x()), 1.0, 1e-12);
		CHECK_NEAR(std::abs(fit.v.y()), 1.0, 1e-12);
		CHECK_NEAR(fit.extentU, 10.0, 1e-12);
		CHECK_NEAR(fit.extentV, 5.0, 1e-12);
		CHECK_NEAR(fit.thickness, 0.0, 1e-12);
		CHECK_NEAR((fit.origin - Vec3(5.0, 2.5, 1.5)).norm(), 0.0, 1e-12);
		CHECK_NEAR(fit.u.cross(fit.v).dot(fit.normal), 1.0, 1e-12);   // right-handed
		CHECK_NEAR(fit.u.dot(fit.v), 0.0, 1e-12);
		CHECK_NEAR(fit.u.norm(), 1.0, 1e-12);
		CHECK_NEAR(fit.v.norm(), 1.0, 1e-12);
		CHECK_NEAR(fit.normal.norm(), 1.0, 1e-12);

		// A receiver line: degenerate but ok, v any perpendicular direction.
		fit = FitPlane(LineSamples(10));
		CHECK(fit.ok);
		CHECK_NEAR(std::abs(fit.u.y()), 1.0, 1e-12);
		CHECK_NEAR(fit.extentU, 9.0, 1e-12);
		CHECK_NEAR(fit.extentV, 0.0, 1e-12);
		CHECK_NEAR(fit.thickness, 0.0, 1e-12);
		CHECK_NEAR(fit.u.dot(fit.v), 0.0, 1e-12);
		CHECK_NEAR(fit.v.norm(), 1.0, 1e-12);
		CHECK_NEAR(fit.u.cross(fit.v).dot(fit.normal), 1.0, 1e-12);

		// Tilted plane spanned by an orthonormal pair.
		const Vec3 e1 = Vec3(1.0, 2.0, 2.0) / 3.0;
		const Vec3 e2 = Vec3(2.0, 1.0, -2.0) / 3.0;
		const Vec3 centre(1.0, -2.0, 3.0);
		std::vector<Vec3> tilted;
		for (int b = 0; b <= 5; b++)
		{
			for (int a = 0; a <= 10; a++) tilted.push_back(centre + double(a) * e1 + double(b) * e2);
		}
		fit = FitPlane(tilted);
		CHECK(fit.ok);
		CHECK_NEAR(std::abs(fit.normal.dot(e1.cross(e2))), 1.0, 1e-9);
		CHECK_NEAR(std::abs(fit.u.dot(e1)), 1.0, 1e-9);
		CHECK_NEAR(fit.extentU, 10.0, 1e-9);
		CHECK_NEAR(fit.extentV, 5.0, 1e-9);
		CHECK_NEAR(fit.thickness, 0.0, 1e-9);

		// Degenerate inputs.
		fit = FitPlane({Vec3(1, 2, 3)});
		CHECK(!fit.ok);
		CHECK_NEAR((fit.origin - Vec3(1, 2, 3)).norm(), 0.0, 0.0);
		CHECK(!FitPlane({}).ok);
		CHECK(!FitPlane({Vec3(1, 2, 3), Vec3(1, 2, 3)}).ok);
		CHECK(FitPlane({Vec3(0, 0, 0), Vec3(1, 1, 1)}).ok);
	}

	void TestInterpolate()
	{
		const std::vector<Vec3> samples = PlaneSamples(11, 6, 0.0);
		const std::vector<Vec3c> random = RandomField(samples.size(), 1);
		std::vector<Vec3c> out;
		std::vector<unsigned char> valid;
		InterpolationSettings idw;
		InterpolationSettings nearest;
		nearest.method = InterpolationMethod::Nearest;

		// Exact recovery at the sample positions, both methods.
		Interpolate(samples, random, samples, idw, out, valid);
		CHECK(out.size() == samples.size() && valid.size() == samples.size());
		bool exact = true, allValid = true;
		for (size_t i = 0; i < samples.size(); i++)
		{
			exact = exact && (out[i] - random[i]).norm() == 0.0;
			allValid = allValid && valid[i] == 1;
		}
		CHECK(exact);
		CHECK(allValid);
		Interpolate(samples, random, samples, nearest, out, valid);
		exact = true;
		for (size_t i = 0; i < samples.size(); i++) exact = exact && (out[i] - random[i]).norm() == 0.0;
		CHECK(exact);

		// A constant field is reproduced everywhere.
		const Vec3c constant(Complex(1.0, 2.0), Complex(3.0, -1.0), Complex(-0.5, 0.25));
		const std::vector<Vec3c> constants(samples.size(), constant);
		const std::vector<Vec3> queries = RandomPoints(200, 10.0, false, 2);
		Interpolate(samples, constants, queries, idw, out, valid);
		double worst = 0.0;
		for (const Vec3c& v : out) worst = std::max(worst, RelativeError(v, constant));
		CHECK(worst < 1e-12);
		Interpolate(samples, constants, queries, nearest, out, valid);
		worst = 0.0;
		for (const Vec3c& v : out) worst = std::max(worst, RelativeError(v, constant));
		CHECK(worst == 0.0);

		// A linear field at edge midpoints: the 8 nearest samples are symmetric
		// about the midpoint on a regular grid, so Shepard reproduces it.
		const std::vector<Vec3c> linear = LinearFieldAt(samples);
		std::vector<Vec3> midpoints;
		for (int j = 1; j <= 4; j++)
		{
			for (int i = 2; i <= 7; i++) midpoints.push_back(Vec3(i + 0.5, double(j), 0.0));
		}
		Interpolate(samples, linear, midpoints, idw, out, valid);
		worst = 0.0;
		for (size_t i = 0; i < midpoints.size(); i++) worst = std::max(worst, RelativeError(out[i], LinearField(midpoints[i])));
		CHECK(worst < 1e-6);

		// Cell centres with the 4 nearest samples (symmetric again).
		InterpolationSettings four;
		four.neighbours = 4;
		std::vector<Vec3> centres;
		for (int j = 1; j <= 3; j++)
		{
			for (int i = 1; i <= 8; i++) centres.push_back(Vec3(i + 0.5, j + 0.5, 0.0));
		}
		Interpolate(samples, linear, centres, four, out, valid);
		worst = 0.0;
		for (size_t i = 0; i < centres.size(); i++) worst = std::max(worst, RelativeError(out[i], LinearField(centres[i])));
		CHECK(worst < 1e-6);

		// Nearest picks the closest sample.
		const std::vector<Vec3> offGrid = {Vec3(3.3, 1.8, 0.1), Vec3(-0.4, 5.4, 0.0), Vec3(9.51, 2.49, 0.0)};
		Interpolate(samples, random, offGrid, nearest, out, valid);
		for (size_t i = 0; i < offGrid.size(); i++)
		{
			CHECK((out[i] - random[NearestIndex(samples, offGrid[i])]).norm() == 0.0);
			CHECK(valid[i] == 1);
		}

		// Beyond the radius: flagged, and the nearest sample's value.
		const std::vector<Vec3> far = {Vec3(30.0, 30.0, 0.0), Vec3(5.0, 2.0, 4.0), Vec3(5.0, 2.0, 2.0)};
		Interpolate(samples, random, far, idw, out, valid);
		CHECK(valid[0] == 0);
		CHECK(valid[1] == 0);
		CHECK(valid[2] == 1);   // 2 m away: within 3 spacings
		CHECK((out[0] - random[NearestIndex(samples, far[0])]).norm() == 0.0);
		CHECK((out[1] - random[NearestIndex(samples, far[1])]).norm() == 0.0);
		Interpolate(samples, random, far, nearest, out, valid);
		CHECK(valid[0] == 0 && valid[1] == 0 && valid[2] == 1);
		CHECK((out[0] - random[NearestIndex(samples, far[0])]).norm() == 0.0);
		InterpolationSettings unlimited;
		unlimited.radiusFactor = 0.0;
		Interpolate(samples, random, far, unlimited, out, valid);
		CHECK(valid[0] == 1 && valid[1] == 1 && valid[2] == 1);

		// One sample: every query gets it.
		const std::vector<Vec3> one = {Vec3(1, 2, 3)};
		const std::vector<Vec3c> oneValue = {constant};
		Interpolate(one, oneValue, queries, idw, out, valid);
		CHECK(out.size() == queries.size());
		exact = true;
		allValid = true;
		for (size_t i = 0; i < out.size(); i++)
		{
			exact = exact && (out[i] - constant).norm() == 0.0;
			allValid = allValid && valid[i] == 1;
		}
		CHECK(exact);
		CHECK(allValid);
		Interpolate(one, oneValue, queries, nearest, out, valid);
		exact = true;
		for (size_t i = 0; i < out.size(); i++) exact = exact && (out[i] - constant).norm() == 0.0;
		CHECK(exact);

		// No samples: empty output.
		Interpolate({}, {}, queries, idw, out, valid);
		CHECK(out.empty());
		CHECK(valid.empty());

		// No queries: empty output of the right sizes.
		Interpolate(samples, random, {}, idw, out, valid);
		CHECK(out.empty());
		CHECK(valid.empty());
	}

	void TestInterpolateOnPlane()
	{
		const std::vector<Vec3> samples = PlaneSamples(11, 6, 1.5);
		const std::vector<Vec3c> linear = LinearFieldAt(samples);
		const PlaneFit plane = FitPlane(samples);
		const InterpolationSettings settings;

		PlaneGrid grid = InterpolateOnPlane(samples, linear, plane, 11, settings);
		CHECK(grid.nu == 11);
		CHECK(grid.nv >= 2);
		CHECK(grid.nv <= 11);
		CHECK(grid.nodes.size() == size_t(grid.nu) * size_t(grid.nv));
		CHECK(grid.values.size() == grid.nodes.size());
		CHECK(grid.valid.size() == grid.nodes.size());
		CHECK(grid.plane.ok);

		// Nodes lie in the plane, cover every sample and are all valid.
		double offPlane = 0.0, nodeSMin = 1e300, nodeSMax = -1e300, nodeTMin = 1e300, nodeTMax = -1e300;
		bool allValid = true;
		for (size_t n = 0; n < grid.nodes.size(); n++)
		{
			const Vec3 d = grid.nodes[n] - plane.origin;
			offPlane = std::max(offPlane, std::abs(d.dot(plane.normal)));
			nodeSMin = std::min(nodeSMin, d.dot(plane.u));
			nodeSMax = std::max(nodeSMax, d.dot(plane.u));
			nodeTMin = std::min(nodeTMin, d.dot(plane.v));
			nodeTMax = std::max(nodeTMax, d.dot(plane.v));
			allValid = allValid && grid.valid[n] == 1;
		}
		CHECK(offPlane < 1e-9);
		CHECK(allValid);
		bool covered = true;
		for (const Vec3& p : samples)
		{
			const Vec3 d = p - plane.origin;
			covered = covered && d.dot(plane.u) >= nodeSMin && d.dot(plane.u) <= nodeSMax
			          && d.dot(plane.v) >= nodeTMin && d.dot(plane.v) <= nodeTMax;
		}
		CHECK(covered);
		// Half a spacing of padding on every side.
		CHECK_NEAR(nodeSMax - nodeSMin, 11.0, 1e-9);
		CHECK_NEAR(nodeTMax - nodeTMin, 6.0, 1e-9);

		// Index layout: i runs along u, j along v.
		CHECK((grid.nodes[grid.Index(1, 0)] - grid.nodes[grid.Index(0, 0)]).dot(plane.u) > 0.0);
		CHECK((grid.nodes[grid.Index(0, 1)] - grid.nodes[grid.Index(0, 0)]).dot(plane.v) > 0.0);
		CHECK_NEAR((grid.nodes[grid.Index(1, 0)] - grid.nodes[grid.Index(0, 0)]).dot(plane.v), 0.0, 1e-9);
		CHECK_NEAR((grid.nodes[grid.Index(0, 1)] - grid.nodes[grid.Index(0, 0)]).dot(plane.u), 0.0, 1e-9);
		CHECK(grid.Index(3, 2) == 2 * 11 + 3);

		// With 0.5 m nodes every sample coincides with a node: exact values.
		grid = InterpolateOnPlane(samples, linear, plane, 23, settings);
		CHECK(grid.nu == 23);
		CHECK(grid.nv == 13);
		int coincident = 0;
		double worst = 0.0;
		for (size_t s = 0; s < samples.size(); s++)
		{
			for (size_t n = 0; n < grid.nodes.size(); n++)
			{
				if ((grid.nodes[n] - samples[s]).norm() < 1e-9)
				{
					coincident++;
					worst = std::max(worst, RelativeError(grid.values[n], linear[s]));
				}
			}
		}
		CHECK(coincident == int(samples.size()));
		CHECK(worst < 1e-9);
		// Interior nodes between samples (at least 1.5 m from the sample
		// boundary): odd i, j coincide with samples; an edge midpoint (one even
		// index) sees a symmetric set of 8 neighbours and reproduces the linear
		// field; a cell centre (both even) takes 4 of 8 equidistant samples,
		// so it is only close.
		double worstEdge = 0.0, worstCentre = 0.0;
		for (int j = 4; j <= 8; j++)
		{
			for (int i = 4; i <= 18; i++)
			{
				const size_t n = grid.Index(i, j);
				const double error = RelativeError(grid.values[n], LinearField(grid.nodes[n]));
				if (i % 2 == 0 && j % 2 == 0)
				{
					worstCentre = std::max(worstCentre, error);
				}
				else
				{
					worstEdge = std::max(worstEdge, error);
				}
			}
		}
		std::printf("  linear field on the 0.5 m plane grid: worst error %.2e at edge midpoints, %.2e at cell centres\n",
		            worstEdge, worstCentre);
		CHECK(worstEdge < 1e-6);
		CHECK(worstCentre < 0.25);   // which 4 of the 8 tied samples win is library-dependent

		// A receiver line gives a strip two nodes wide.
		const std::vector<Vec3> line = LineSamples(10);
		const PlaneFit lineFit = FitPlane(line);
		grid = InterpolateOnPlane(line, LinearFieldAt(line), lineFit, 10, settings);
		CHECK(grid.nu == 10);
		CHECK(grid.nv == 2);
		CHECK(grid.nodes.size() == 20);
		offPlane = 0.0;
		allValid = true;
		for (size_t n = 0; n < grid.nodes.size(); n++)
		{
			offPlane = std::max(offPlane, std::abs((grid.nodes[n] - lineFit.origin).dot(lineFit.normal)));
			allValid = allValid && grid.valid[n] == 1;
		}
		CHECK(offPlane < 1e-9);
		CHECK(allValid);
		CHECK_NEAR((grid.nodes[grid.Index(0, 1)] - grid.nodes[grid.Index(0, 0)]).norm(), 1.0, 1e-9);

		// Degenerate inputs.
		grid = InterpolateOnPlane({}, {}, plane, 11, settings);
		CHECK(grid.nu == 0 && grid.nv == 0 && grid.nodes.empty() && grid.values.empty() && grid.valid.empty());
		const std::vector<Vec3> one = {Vec3(1, 2, 3)};
		grid = InterpolateOnPlane(one, {linear[0]}, FitPlane(one), 5, settings);
		CHECK(grid.nu == 5 && grid.nv == 5);
		CHECK(grid.values.size() == 25 && (grid.values[12] - linear[0]).norm() == 0.0);
		grid = InterpolateOnPlane(samples, linear, plane, 1, settings);   // resolution clamped to 2
		CHECK(grid.nu == 2 && grid.nv == 2);
	}

	void TestInterpolateInVolume()
	{
		const InterpolationSettings settings;

		// A flat layout: a single layer of nodes at the samples' height.
		std::vector<Vec3> samples = PlaneSamples(11, 6, 1.5);
		std::vector<Vec3c> linear = LinearFieldAt(samples);
		VolumeGrid grid = InterpolateInVolume(samples, linear, 11, settings);
		CHECK(grid.nx == 11);
		CHECK(grid.ny >= 2 && grid.ny <= 11);
		CHECK(grid.nz == 1);
		CHECK(grid.nodes.size() == size_t(grid.nx) * size_t(grid.ny) * size_t(grid.nz));
		CHECK(grid.values.size() == grid.nodes.size() && grid.valid.size() == grid.nodes.size());
		CHECK_NEAR((grid.boundMin - Vec3(-0.5, -0.5, 1.0)).norm(), 0.0, 1e-12);
		CHECK_NEAR((grid.boundMax - Vec3(10.5, 5.5, 2.0)).norm(), 0.0, 1e-12);
		double worstZ = 0.0;
		bool allValid = true;
		for (size_t n = 0; n < grid.nodes.size(); n++)
		{
			worstZ = std::max(worstZ, std::abs(grid.nodes[n].z() - 1.5));
			allValid = allValid && grid.valid[n] == 1;
		}
		CHECK(worstZ < 1e-12);
		CHECK(allValid);

		// A box: proportional node counts, index layout matching the axes.
		samples = BoxSamples(11, 5, 3);
		linear = LinearFieldAt(samples);
		grid = InterpolateInVolume(samples, linear, 11, settings);
		CHECK(grid.nx == 11);
		CHECK(grid.ny >= 5 && grid.ny <= 6);
		CHECK(grid.nz >= 3 && grid.nz <= 4);
		const double step = (grid.boundMax.x() - grid.boundMin.x()) / double(grid.nx - 1);
		CHECK(std::abs((grid.ny - 1) * step - (grid.boundMax.y() - grid.boundMin.y())) <= step);
		CHECK(std::abs((grid.nz - 1) * step - (grid.boundMax.z() - grid.boundMin.z())) <= step);
		CHECK(grid.Index(2, 3, 1) == (size_t(1) * size_t(grid.ny) + 3) * size_t(grid.nx) + 2);
		bool layout = true;
		for (int k = 0; k < grid.nz; k++)
		{
			for (int j = 0; j < grid.ny; j++)
			{
				for (int i = 0; i < grid.nx; i++)
				{
					const Vec3& node = grid.nodes[grid.Index(i, j, k)];
					if (i > 0) layout = layout && node.x() > grid.nodes[grid.Index(i - 1, j, k)].x();
					if (j > 0) layout = layout && node.y() > grid.nodes[grid.Index(i, j - 1, k)].y();
					if (k > 0) layout = layout && node.z() > grid.nodes[grid.Index(i, j, k - 1)].z();
					layout = layout && node.x() == grid.nodes[grid.Index(i, 0, 0)].x();
					layout = layout && node.y() == grid.nodes[grid.Index(0, j, 0)].y();
					layout = layout && node.z() == grid.nodes[grid.Index(0, 0, k)].z();
				}
			}
		}
		CHECK(layout);
		CHECK_NEAR(grid.nodes[grid.Index(0, 0, 0)].x(), -0.5, 1e-12);
		CHECK_NEAR(grid.nodes[grid.Index(grid.nx - 1, grid.ny - 1, grid.nz - 1)].z(), 2.5, 1e-12);

		// 0.5 m nodes: every sample coincides with a node, exact values there.
		grid = InterpolateInVolume(samples, linear, 23, settings);
		CHECK(grid.nx == 23 && grid.ny == 11 && grid.nz == 7);
		int coincident = 0;
		double worst = 0.0;
		for (size_t s = 0; s < samples.size(); s++)
		{
			const Vec3 q = samples[s];
			const int i = int(std::lround((q.x() - grid.boundMin.x()) / 0.5));
			const int j = int(std::lround((q.y() - grid.boundMin.y()) / 0.5));
			const int k = int(std::lround((q.z() - grid.boundMin.z()) / 0.5));
			const size_t n = grid.Index(i, j, k);
			if ((grid.nodes[n] - q).norm() < 1e-9)
			{
				coincident++;
				worst = std::max(worst, RelativeError(grid.values[n], linear[s]));
			}
		}
		CHECK(coincident == int(samples.size()));
		CHECK(worst < 1e-9);

		// A receiver line: n x 1 x 1.
		const std::vector<Vec3> line = LineSamples(10);
		grid = InterpolateInVolume(line, LinearFieldAt(line), 10, settings);
		CHECK(grid.nx == 1 && grid.ny == 10 && grid.nz == 1);
		CHECK(grid.nodes.size() == 10);
		CHECK_NEAR(grid.nodes[0].x(), 1.0, 1e-12);
		CHECK_NEAR(grid.nodes[0].z(), 2.0, 1e-12);
		CHECK_NEAR(grid.nodes[0].y(), -0.5, 1e-12);
		CHECK_NEAR(grid.nodes[9].y(), 9.5, 1e-12);

		// Degenerate inputs.
		grid = InterpolateInVolume({}, {}, 11, settings);
		CHECK(grid.nx == 0 && grid.ny == 0 && grid.nz == 0 && grid.nodes.empty());
		const std::vector<Vec3> one = {Vec3(1, 2, 3)};
		grid = InterpolateInVolume(one, {linear[0]}, 5, settings);
		CHECK(grid.nodes.size() == 5 && (grid.values[2] - linear[0]).norm() == 0.0 && grid.valid[2] == 1);
	}

	void TestPerformance()
	{
		const std::vector<Vec3> flat = RandomPoints(3000, 10.0, true, 3);
		const std::vector<Vec3c> flatField = RandomField(flat.size(), 4);
		const InterpolationSettings settings;

		Clock::time_point start = Clock::now();
		const PlaneFit plane = FitPlane(flat);
		const PlaneGrid planeGrid = InterpolateOnPlane(flat, flatField, plane, 100, settings);
		const double planeSeconds = Seconds(start);
		std::printf("  3000 samples, %d x %d plane nodes: %.3f s\n", planeGrid.nu, planeGrid.nv, planeSeconds);
		CHECK(plane.ok);
		CHECK(planeGrid.nu == 100 || planeGrid.nv == 100);
		CHECK(planeGrid.nodes.size() >= 100 * 90);
		CHECK(planeSeconds < 2.0);

		const std::vector<Vec3> box = RandomPoints(3000, 10.0, false, 5);
		const std::vector<Vec3c> boxField = RandomField(box.size(), 6);
		start = Clock::now();
		const VolumeGrid volume = InterpolateInVolume(box, boxField, 40, settings);
		const double volumeSeconds = Seconds(start);
		std::printf("  3000 samples, %d x %d x %d volume nodes: %.3f s\n", volume.nx, volume.ny, volume.nz, volumeSeconds);
		CHECK(volume.nodes.size() >= size_t(40 * 36 * 36));
		// Loose bounds for the larger cases: they only catch an O(n^2) search
		// (minutes), not the slowdown of an unoptimised build.
		CHECK(volumeSeconds < 10.0);

		// Spacing of a large layout is O(n log n): 1e5 samples.
		const std::vector<Vec3> many = RandomPoints(100000, 100.0, false, 7);
		start = Clock::now();
		const double spacing = SampleSpacing(many);
		const double spacingSeconds = Seconds(start);
		std::printf("  spacing of 100000 samples (%.4f m): %.3f s\n", spacing, spacingSeconds);
		CHECK(spacing > 0.0);
		CHECK(spacingSeconds < 10.0);
	}
}

int main()
{
	Run("BuildSceneMesh", TestBuildSceneMesh);
	Run("BuildSceneEdges", TestBuildSceneEdges);
	Run("FieldQuantities", TestFieldQuantities);
	Run("SampleSpacing", TestSampleSpacing);
	Run("FitPlane", TestFitPlane);
	Run("Interpolate", TestInterpolate);
	Run("InterpolateOnPlane", TestInterpolateOnPlane);
	Run("InterpolateInVolume", TestInterpolateInVolume);
	Run("Performance", TestPerformance);

	std::printf("%d checks, %d failures\n", g_checks, g_failures);
	return g_failures == 0 ? 0 : 1;
}
