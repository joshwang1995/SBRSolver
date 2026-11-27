// Unit tests of the capture stage of the CPU tracer (CpuTracer.h): the
// receiver grid must give the same candidates as a test of every receiver
// against every segment. The harness is the one of viewer_tests.cpp.
#include <algorithm>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

#include "CpuTracer.h"
#include "common/Constants.h"

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

	void Run(const char* name, void (*test)())
	{
		g_test = name;
		const int before = g_failures;
		std::printf("[ RUN  ] %s\n", name);
		test();
		std::printf("[ %s ] %s\n", g_failures == before ? " OK " : "FAIL", name);
	}

	const unsigned kSeed = 1234;
	const int kSegments = 500;
	const int kReceivers = 2000;

	enum class Layout
	{
		Box,
		Plane,
		Line,
		Point
	};

	const char* LayoutName(Layout layout)
	{
		switch (layout)
		{
		case Layout::Box: return "box";
		case Layout::Plane: return "plane";
		case Layout::Line: return "line";
		case Layout::Point: return "point";
		}
		return "";
	}

	VecVec3 MakeReceivers(Layout layout, std::mt19937& random)
	{
		std::uniform_real_distribution<double> coordinate(-10.0, 10.0);
		VecVec3 receivers;
		for (int i = 0; i < kReceivers; i++)
		{
			Vec3 r(coordinate(random), coordinate(random), coordinate(random));
			if (layout == Layout::Plane) r.z() = 1.5;
			if (layout == Layout::Line) r = Vec3(3.0, r.y(), 2.0);
			if (layout == Layout::Point) r = Vec3(1.0, 2.0, 3.0);
			receivers.push_back(r);
		}
		return receivers;
	}

	// Random finite and escaping segments, some of them children of an
	// earlier finite segment, so that the unfolded start length varies.
	std::vector<Segment> MakeSegments(std::mt19937& random)
	{
		std::uniform_real_distribution<double> coordinate(-30.0, 30.0);
		std::uniform_real_distribution<double> unit(-1.0, 1.0);
		std::uniform_real_distribution<double> length(0.5, 40.0);
		std::uniform_real_distribution<double> chance(0.0, 1.0);
		std::vector<Segment> segments;
		std::vector<int> finite;
		for (int i = 0; i < kSegments; i++)
		{
			Segment segment;
			segment.launch = i % 3;
			Vec3 direction(unit(random), unit(random), unit(random));
			if (direction.norm() == 0.0) direction = Vec3(1.0, 0.0, 0.0);
			direction.normalize();
			double unfoldedStart = 0.0;
			if (!finite.empty() && chance(random) < 0.5)
			{
				segment.parent = finite[size_t(chance(random) * double(finite.size())) % finite.size()];
				unfoldedStart = segments[segment.parent].ray.pathLength;
				segment.ray.sourcePoint = segments[segment.parent].ray.targetPoint;
			}
			else
			{
				segment.ray.sourcePoint = Vec3(coordinate(random), coordinate(random), coordinate(random));
			}
			if (chance(random) < 0.6)
			{
				const double l = length(random);
				segment.ray.targetPoint = segment.ray.sourcePoint + direction * l;
				segment.ray.pathLength = unfoldedStart + l;
				finite.push_back(i);
			}
			else
			{
				segment.ray.targetPoint = segment.ray.sourcePoint + direction * 1E10;
				segment.ray.pathLength = INF;
			}
			segments.push_back(segment);
		}
		return segments;
	}

	void CheckLayout(Layout layout)
	{
		std::mt19937 random(kSeed);
		const VecVec3 receivers = MakeReceivers(layout, random);
		const std::vector<TransmitterSource> transmitters(2);
		const std::vector<int> launchTx = { 0, 1, 0 };
		const std::vector<Vec3> launches(3, Vec3(0.0, 0.0, 1.0));

		TracerInput input;
		input.transmitters = &transmitters;
		input.launches = &launches;
		input.launchTx = &launchTx;
		input.receivers = &receivers;
		input.tanBeta = { ReceptionSlope(3), ReceptionSlope(5) };
		input.threads = 3;

		TraceResult result;
		result.segments = MakeSegments(random);
		CaptureCandidates(input, result);

		std::vector<std::pair<int, int>> expected;
		for (int s = 0; s < int(result.segments.size()); s++)
		{
			const Segment& segment = result.segments[s];
			const double unfoldedStart = (segment.parent < 0) ? 0.0 : result.segments[segment.parent].ray.pathLength;
			const double tanBeta = input.tanBeta[launchTx[segment.launch]];
			for (int r = 0; r < int(receivers.size()); r++)
			{
				Vec3 capturePoint;
				if (HitReceptionSphere(segment.ray.sourcePoint, segment.ray.targetPoint, receivers[r], unfoldedStart, tanBeta, capturePoint))
				{
					expected.push_back(std::make_pair(r, s));
				}
			}
		}
		std::sort(expected.begin(), expected.end());

		std::printf("  %s: %zu candidates of %zu pairs (seed %u)\n", LayoutName(layout), expected.size(),
			result.segments.size() * receivers.size(), kSeed);
		// A layout that captures nothing, or everything, would not test the cull.
		CHECK(!expected.empty());
		CHECK(expected.size() < result.segments.size() * receivers.size());
		CHECK(result.candidates == expected);
	}

	void TestBox() { CheckLayout(Layout::Box); }
	void TestPlane() { CheckLayout(Layout::Plane); }
	void TestLine() { CheckLayout(Layout::Line); }
	void TestPoint() { CheckLayout(Layout::Point); }
}

int main()
{
	Run("GridCaptureBox", TestBox);
	Run("GridCapturePlane", TestPlane);
	Run("GridCaptureLine", TestLine);
	Run("GridCapturePoint", TestPoint);
	std::printf("%d checks, %d failures\n", g_checks, g_failures);
	return g_failures == 0 ? 0 : 1;
}
