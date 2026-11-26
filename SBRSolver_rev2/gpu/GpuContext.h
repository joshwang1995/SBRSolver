#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <volk.h>

#include "gpu/GpuTypes.h"

// The Vulkan host layer of the GPU tracer: one device, one compute queue,
// storage buffers, compute pipelines and chunked dispatches. Vulkan 1.1 core
// only, plus the optional 64-bit float feature. Every Vulkan object is
// created here and nowhere else, so a reader can follow GpuContext.cpp from
// top to bottom.
namespace sbr
{
namespace gpu
{

// No usable Vulkan device: the solver exits with 77 instead of failing.
class DeviceUnavailable : public std::runtime_error
{
public:
	using std::runtime_error::runtime_error;
};

// What the chosen device reports.
struct DeviceInfo
{
	std::string name;
	std::string type;                   // discrete, integrated, virtual, cpu or other
	uint32_t apiVersion = 0;
	std::string driver;                 // driver name and version
	bool float64 = false;               // shaders may compute in 64-bit floats
	uint64_t storageRange = 0;          // largest storage buffer binding, bytes
	uint64_t allocationLimit = 0;       // largest single memory allocation, bytes
	uint64_t heapBytes = 0;             // largest device-local memory heap, bytes
	uint32_t workGroupCount[3] = { 0, 0, 0 }; // dispatch limit per dimension, in work groups
	std::string validation;             // what happened to SBR_VULKAN_VALIDATION, empty when unset
};

// A storage buffer in device memory, freed with the object. Destroy every
// buffer before the context that created it.
class Buffer
{
public:
	Buffer() = default;
	~Buffer();
	Buffer(Buffer&& other) noexcept;
	Buffer& operator=(Buffer&& other) noexcept;
	Buffer(const Buffer&) = delete;
	Buffer& operator=(const Buffer&) = delete;

	uint64_t Bytes() const { return _bytes; }

private:
	friend class GpuContext;
	void Release();

	VkDevice _device = VK_NULL_HANDLE;
	VkBuffer _buffer = VK_NULL_HANDLE;
	VkDeviceMemory _memory = VK_NULL_HANDLE;
	uint64_t _bytes = 0;
};

// A compute shader with its storage buffer bindings and its push constants.
class Pipeline
{
public:
	Pipeline() = default;
	~Pipeline();
	Pipeline(Pipeline&& other) noexcept;
	Pipeline& operator=(Pipeline&& other) noexcept;
	Pipeline(const Pipeline&) = delete;
	Pipeline& operator=(const Pipeline&) = delete;

private:
	friend class GpuContext;
	void Release();

	VkDevice _device = VK_NULL_HANDLE;
	VkDescriptorSetLayout _setLayout = VK_NULL_HANDLE;
	VkPipelineLayout _layout = VK_NULL_HANDLE;
	VkPipeline _pipeline = VK_NULL_HANDLE;
	uint32_t _bindings = 0;
	uint32_t _pushBytes = 0;
};

// The first two values of every push constant block (SBR_DISPATCH_BASE in
// common.glsl). Dispatch writes them for every chunk.
struct DispatchBase
{
	uint32_t firstX = 0;
	uint32_t firstY = 0;
};

class GpuContext
{
public:
	// Work group size of every shader in x (local_size_x).
	static const uint32_t kWorkGroupSize = 64;

	/**
	 * @brief Loads the Vulkan loader, creates the instance and opens a device.
	 * @param skipCpuDevices Ignore devices of type CPU, such as lavapipe.
	 * @return The context; its device info is ready.
	 * @note Throws DeviceUnavailable when there is no loader, no device, or
	 * no device that matches SBR_VULKAN_DEVICE.
	 */
	static std::unique_ptr<GpuContext> Create(bool skipCpuDevices);

	GpuContext() = default;
	~GpuContext();
	GpuContext(const GpuContext&) = delete;
	GpuContext& operator=(const GpuContext&) = delete;

	const DeviceInfo& Info() const { return _info; }

	/**
	 * @brief The precision a run computes in: Auto gives Double when the
	 * device has 64-bit floats, else Single.
	 * @note Throws std::runtime_error for Double on a device without it.
	 */
	GpuPrecision Resolve(GpuPrecision requested) const;

	/** @brief One line for the log on why Resolve picked what it picked. */
	std::string PrecisionReason(GpuPrecision requested) const;

	/** @brief A storage buffer of `bytes`; throws when it exceeds a device limit. */
	Buffer CreateBuffer(uint64_t bytes);

	/** @brief Copies host memory into a buffer, in chunks through a staging buffer. */
	void Upload(Buffer& buffer, const void* data, uint64_t bytes, uint64_t offset = 0);

	/** @brief Copies a buffer into host memory, in chunks through a staging buffer. */
	void Download(const Buffer& buffer, void* data, uint64_t bytes, uint64_t offset = 0);

	/** @brief Sets every 32-bit word of a buffer to `value`. */
	void Fill(Buffer& buffer, uint32_t value);

	/**
	 * @brief A compute pipeline from embedded SPIR-V.
	 * @param bindings Storage buffers the shader binds, as set 0, bindings 0 to bindings-1.
	 * @param pushBytes Size of the shader's push constant block, DispatchBase included.
	 */
	Pipeline CreateCompute(const uint32_t* spirv, size_t words, uint32_t bindings, uint32_t pushBytes);

	/**
	 * @brief Runs a pipeline over itemsX by itemsY items and waits for it.
	 * The dispatch runs in chunks the device accepts; each chunk's first
	 * item is written into `push`, which starts with a DispatchBase.
	 * @param buffers One per binding, in binding order.
	 * @param push The push constant block of pipeline's pushBytes.
	 */
	void Dispatch(const Pipeline& pipeline, const std::vector<const Buffer*>& buffers, void* push, uint32_t pushBytes, uint64_t itemsX, uint64_t itemsY = 1);

private:
	void CreateInstance();
	void PickDevice(bool skipCpuDevices);
	void CreateDevice();
	void EnsureStaging();
	uint32_t MemoryType(uint32_t typeBits, VkMemoryPropertyFlags required) const;
	void Begin();
	void Barrier();
	void Submit();

	VkInstance _instance = VK_NULL_HANDLE;
	VkPhysicalDevice _physical = VK_NULL_HANDLE;
	VkDevice _device = VK_NULL_HANDLE;
	VkQueue _queue = VK_NULL_HANDLE;
	uint32_t _queueFamily = 0;
	VkCommandPool _commandPool = VK_NULL_HANDLE;
	VkCommandBuffer _commands = VK_NULL_HANDLE;
	VkFence _fence = VK_NULL_HANDLE;
	VkDescriptorPool _descriptorPool = VK_NULL_HANDLE;
	VkBuffer _staging = VK_NULL_HANDLE;
	VkDeviceMemory _stagingMemory = VK_NULL_HANDLE;
	void* _stagingMapped = nullptr;
	VkPhysicalDeviceMemoryProperties _memoryProperties{};
	DeviceInfo _info;
};

} // namespace gpu
} // namespace sbr
