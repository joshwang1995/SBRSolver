// Unit tests of the Vulkan host layer (gpu/GpuContext.h) with the probe
// shader: buffers, uploads, downloads, atomics and the chunked dispatch in
// both precisions. Exit code 77 when no Vulkan device exists, so ctest
// reports the test as skipped. The harness is the one of viewer_tests.cpp.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "gpu/GpuContext.h"
#include "probe_double.spv.h"
#include "probe_single.spv.h"

using sbr::gpu::Buffer;
using sbr::gpu::DeviceUnavailable;
using sbr::gpu::DispatchBase;
using sbr::gpu::GpuContext;
using sbr::gpu::Pipeline;

namespace
{
	int g_checks = 0;
	int g_failures = 0;
	const char* g_test = "";
	std::unique_ptr<GpuContext> g_gpu;

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

	struct ProbePush
	{
		DispatchBase base;
		uint32_t count;
	};

	// Runs the probe over `count` items and checks every value and the
	// appended index set. T is float or double, matching the shader build.
	template <typename T>
	void CheckProbe(const uint32_t* spirv, size_t words, uint32_t count)
	{
		Buffer values = g_gpu->CreateBuffer(uint64_t(count) * sizeof(T));
		Buffer appended = g_gpu->CreateBuffer(sizeof(uint32_t) * (uint64_t(count) + 1));
		g_gpu->Fill(appended, 0);
		Pipeline probe = g_gpu->CreateCompute(spirv, words, 2, sizeof(ProbePush));
		ProbePush push{};
		push.count = count;
		g_gpu->Dispatch(probe, { &values, &appended }, &push, sizeof(push), count);

		std::vector<T> result(count);
		g_gpu->Download(values, result.data(), uint64_t(count) * sizeof(T));
		size_t wrong = 0;
		for (uint32_t i = 0; i < count; i++)
		{
			const T expected = T(i) * T(0.25) + T(1.0);
			if (std::memcmp(&result[i], &expected, sizeof(T)) != 0) wrong++;
		}
		CHECK(wrong == 0);

		std::vector<uint32_t> items(count + 1);
		g_gpu->Download(appended, items.data(), sizeof(uint32_t) * (uint64_t(count) + 1));
		CHECK(items[0] == count);
		std::vector<uint32_t> indices(items.begin() + 1, items.end());
		std::sort(indices.begin(), indices.end());
		std::vector<uint32_t> expectedIndices(count);
		std::iota(expectedIndices.begin(), expectedIndices.end(), 0u);
		CHECK(indices == expectedIndices);
		std::printf("  %u items, %zu wrong values, counter %u\n", count, wrong, items[0]);
	}

	void TestProbeSingle()
	{
		CheckProbe<float>(kSpirv_probe_single, kSpirvWords_probe_single, 100000);
	}

	void TestProbeDouble()
	{
		if (!g_gpu->Info().float64)
		{
			std::printf("  skipped: the device has no 64-bit floats\n");
			return;
		}
		CheckProbe<double>(kSpirv_probe_double, kSpirvWords_probe_double, 100000);
	}

	void TestProbeChunks()
	{
		// More items than one dispatch may hold on this device, so the
		// chunk loop of GpuContext::Dispatch runs at least twice.
		const uint64_t perDispatch = uint64_t(g_gpu->Info().workGroupCount[0]) * GpuContext::kWorkGroupSize;
		const uint32_t count = 4200000;
		CHECK(count > perDispatch);
		CheckProbe<float>(kSpirv_probe_single, kSpirvWords_probe_single, count);
	}

	void TestResolve()
	{
		const bool float64 = g_gpu->Info().float64;
		CHECK(g_gpu->Resolve(GpuPrecision::Auto) == (float64 ? GpuPrecision::Double : GpuPrecision::Single));
		CHECK(g_gpu->Resolve(GpuPrecision::Single) == GpuPrecision::Single);
		bool threw = false;
		try
		{
			CHECK(g_gpu->Resolve(GpuPrecision::Double) == GpuPrecision::Double);
		}
		catch (const std::runtime_error&)
		{
			threw = true;
		}
		CHECK(threw == !float64);
	}

	void TestBufferLimit()
	{
		const uint64_t tooLarge = g_gpu->Info().storageRange + 4;
		std::string message;
		try
		{
			Buffer buffer = g_gpu->CreateBuffer(tooLarge);
		}
		catch (const std::runtime_error& e)
		{
			message = e.what();
		}
		CHECK(message.find(std::to_string(g_gpu->Info().storageRange)) != std::string::npos);
	}
}

int main()
{
	try
	{
		g_gpu = GpuContext::Create(false);
	}
	catch (const DeviceUnavailable& e)
	{
		std::printf("SKIP: %s\n", e.what());
		return 77;
	}
	const sbr::gpu::DeviceInfo& info = g_gpu->Info();
	std::printf("device: %s (%s, Vulkan %u.%u, driver %s, float64 %s, work groups %u x %u, storage range %llu)\n",
		info.name.c_str(), info.type.c_str(), VK_VERSION_MAJOR(info.apiVersion), VK_VERSION_MINOR(info.apiVersion),
		info.driver.c_str(), info.float64 ? "yes" : "no", info.workGroupCount[0], info.workGroupCount[1],
		static_cast<unsigned long long>(info.storageRange));
	if (!info.validation.empty()) std::printf("note: %s\n", info.validation.c_str());

	Run("ProbeSingle", TestProbeSingle);
	Run("ProbeDouble", TestProbeDouble);
	Run("ProbeChunks", TestProbeChunks);
	Run("Resolve", TestResolve);
	Run("BufferLimit", TestBufferLimit);
	std::printf("%d checks, %d failures\n", g_checks, g_failures);
	g_gpu.reset();
	return g_failures == 0 ? 0 : 1;
}
