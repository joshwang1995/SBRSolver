#include "gpu/GpuContext.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>

namespace sbr
{
namespace gpu
{
namespace
{
	// Transfers between the host and a buffer go through a staging buffer of
	// this size, so a large buffer never needs a host copy of its own size.
	const uint64_t kStagingBytes = uint64_t(64) << 20;
	const char* kValidationLayer = "VK_LAYER_KHRONOS_validation";
	// Storage buffers one dispatch may bind.
	const uint32_t kMaxBindings = 16;

	void Check(VkResult result, const char* what)
	{
		if (result != VK_SUCCESS)
		{
			std::ostringstream message;
			message << what << " failed (VkResult " << result << ")";
			throw std::runtime_error(message.str());
		}
	}

	std::string ToLower(std::string text)
	{
		std::transform(text.begin(), text.end(), text.begin(),
			[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		return text;
	}

	const char* TypeName(VkPhysicalDeviceType type)
	{
		switch (type)
		{
		case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU: return "discrete";
		case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU: return "integrated";
		case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU: return "virtual";
		case VK_PHYSICAL_DEVICE_TYPE_CPU: return "cpu";
		default: return "other";
		}
	}

	// Order of preference when no device is named: a discrete GPU first.
	int TypeRank(VkPhysicalDeviceType type)
	{
		switch (type)
		{
		case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU: return 0;
		case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU: return 1;
		case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU: return 2;
		case VK_PHYSICAL_DEVICE_TYPE_CPU: return 3;
		default: return 4;
		}
	}

	std::string VersionText(uint32_t version)
	{
		std::ostringstream text;
		text << VK_VERSION_MAJOR(version) << "." << VK_VERSION_MINOR(version) << "." << VK_VERSION_PATCH(version);
		return text.str();
	}

	const char* kNoDevice = "no Vulkan device (install the graphics driver, or run with --device cpu)";
}

// ---------------------------------------------------------------------------
// Buffer and Pipeline: handles that free themselves
// ---------------------------------------------------------------------------

Buffer::~Buffer()
{
	Release();
}

Buffer::Buffer(Buffer&& other) noexcept
{
	*this = std::move(other);
}

Buffer& Buffer::operator=(Buffer&& other) noexcept
{
	if (this != &other)
	{
		Release();
		_device = other._device;
		_buffer = other._buffer;
		_memory = other._memory;
		_bytes = other._bytes;
		other._buffer = VK_NULL_HANDLE;
		other._memory = VK_NULL_HANDLE;
		other._bytes = 0;
	}
	return *this;
}

void Buffer::Release()
{
	if (_buffer != VK_NULL_HANDLE) vkDestroyBuffer(_device, _buffer, nullptr);
	if (_memory != VK_NULL_HANDLE) vkFreeMemory(_device, _memory, nullptr);
	_buffer = VK_NULL_HANDLE;
	_memory = VK_NULL_HANDLE;
	_bytes = 0;
}

Pipeline::~Pipeline()
{
	Release();
}

Pipeline::Pipeline(Pipeline&& other) noexcept
{
	*this = std::move(other);
}

Pipeline& Pipeline::operator=(Pipeline&& other) noexcept
{
	if (this != &other)
	{
		Release();
		_device = other._device;
		_setLayout = other._setLayout;
		_layout = other._layout;
		_pipeline = other._pipeline;
		_bindings = other._bindings;
		_pushBytes = other._pushBytes;
		other._setLayout = VK_NULL_HANDLE;
		other._layout = VK_NULL_HANDLE;
		other._pipeline = VK_NULL_HANDLE;
	}
	return *this;
}

void Pipeline::Release()
{
	if (_pipeline != VK_NULL_HANDLE) vkDestroyPipeline(_device, _pipeline, nullptr);
	if (_layout != VK_NULL_HANDLE) vkDestroyPipelineLayout(_device, _layout, nullptr);
	if (_setLayout != VK_NULL_HANDLE) vkDestroyDescriptorSetLayout(_device, _setLayout, nullptr);
	_pipeline = VK_NULL_HANDLE;
	_layout = VK_NULL_HANDLE;
	_setLayout = VK_NULL_HANDLE;
}

// ---------------------------------------------------------------------------
// Instance and device
// ---------------------------------------------------------------------------

std::unique_ptr<GpuContext> GpuContext::Create(bool skipCpuDevices)
{
	// volk finds libvulkan.so.1 or vulkan-1.dll; without it there is no
	// Vulkan on this machine.
	if (volkInitialize() != VK_SUCCESS)
	{
		throw DeviceUnavailable("no Vulkan loader (install the graphics driver, or run with --device cpu)");
	}
	std::unique_ptr<GpuContext> context = std::make_unique<GpuContext>();
	context->CreateInstance();
	context->PickDevice(skipCpuDevices);
	context->CreateDevice();
	return context;
}

GpuContext::~GpuContext()
{
	if (_device != VK_NULL_HANDLE)
	{
		vkDeviceWaitIdle(_device);
		if (_staging != VK_NULL_HANDLE) vkDestroyBuffer(_device, _staging, nullptr);
		if (_stagingMemory != VK_NULL_HANDLE) vkFreeMemory(_device, _stagingMemory, nullptr);
		if (_descriptorPool != VK_NULL_HANDLE) vkDestroyDescriptorPool(_device, _descriptorPool, nullptr);
		if (_fence != VK_NULL_HANDLE) vkDestroyFence(_device, _fence, nullptr);
		if (_commandPool != VK_NULL_HANDLE) vkDestroyCommandPool(_device, _commandPool, nullptr);
		vkDestroyDevice(_device, nullptr);
	}
	if (_instance != VK_NULL_HANDLE) vkDestroyInstance(_instance, nullptr);
}

void GpuContext::CreateInstance()
{
	// The instance is the connection to the loader. The validation layer,
	// when asked for and installed, checks every call and prints to stdout.
	std::vector<const char*> layers;
	const char* wantValidation = std::getenv("SBR_VULKAN_VALIDATION");
	if (wantValidation != nullptr && std::string(wantValidation) == "1")
	{
		uint32_t count = 0;
		vkEnumerateInstanceLayerProperties(&count, nullptr);
		std::vector<VkLayerProperties> available(count);
		vkEnumerateInstanceLayerProperties(&count, available.data());
		const bool installed = std::any_of(available.begin(), available.end(),
			[](const VkLayerProperties& layer) { return std::string(layer.layerName) == kValidationLayer; });
		if (installed)
		{
			layers.push_back(kValidationLayer);
			_info.validation = std::string(kValidationLayer) + " loaded";
		}
		else
		{
			_info.validation = std::string("SBR_VULKAN_VALIDATION=1 but ") + kValidationLayer + " is not installed";
		}
	}

	VkApplicationInfo application{};
	application.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
	application.pApplicationName = "sbrsolver";
	application.apiVersion = VK_API_VERSION_1_1;

	VkInstanceCreateInfo info{};
	info.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
	info.pApplicationInfo = &application;
	info.enabledLayerCount = uint32_t(layers.size());
	info.ppEnabledLayerNames = layers.data();
	if (vkCreateInstance(&info, nullptr, &_instance) != VK_SUCCESS)
	{
		throw DeviceUnavailable(kNoDevice);
	}
	volkLoadInstanceOnly(_instance);
}

void GpuContext::PickDevice(bool skipCpuDevices)
{
	uint32_t count = 0;
	Check(vkEnumeratePhysicalDevices(_instance, &count, nullptr), "vkEnumeratePhysicalDevices");
	std::vector<VkPhysicalDevice> devices(count);
	if (count > 0) Check(vkEnumeratePhysicalDevices(_instance, &count, devices.data()), "vkEnumeratePhysicalDevices");
	if (devices.empty()) throw DeviceUnavailable(kNoDevice);

	std::vector<VkPhysicalDeviceProperties> properties(devices.size());
	std::string names;
	for (size_t i = 0; i < devices.size(); i++)
	{
		vkGetPhysicalDeviceProperties(devices[i], &properties[i]);
		names += (i > 0 ? ", " : "") + std::to_string(i) + ": " + properties[i].deviceName;
	}

	// SBR_VULKAN_DEVICE names the device by index or by part of its name;
	// otherwise the best type wins.
	int chosen = -1;
	const char* selector = std::getenv("SBR_VULKAN_DEVICE");
	if (selector != nullptr && *selector != '\0')
	{
		const std::string wanted = selector;
		const bool isIndex = std::all_of(wanted.begin(), wanted.end(), [](unsigned char c) { return std::isdigit(c) != 0; });
		if (isIndex)
		{
			const long index = std::strtol(wanted.c_str(), nullptr, 10);
			if (index >= 0 && index < long(devices.size())) chosen = int(index);
		}
		else
		{
			for (size_t i = 0; i < devices.size() && chosen < 0; i++)
			{
				if (ToLower(properties[i].deviceName).find(ToLower(wanted)) != std::string::npos) chosen = int(i);
			}
		}
		if (chosen < 0)
		{
			throw DeviceUnavailable("no Vulkan device matches SBR_VULKAN_DEVICE=" + wanted + " (devices: " + names
				+ "; unset the variable or pick one of them)");
		}
	}
	else
	{
		bool cpuOnly = false;
		int bestRank = std::numeric_limits<int>::max();
		for (size_t i = 0; i < devices.size(); i++)
		{
			if (properties[i].apiVersion < VK_API_VERSION_1_1) continue;
			if (skipCpuDevices && properties[i].deviceType == VK_PHYSICAL_DEVICE_TYPE_CPU)
			{
				cpuOnly = true;
				continue;
			}
			const int rank = TypeRank(properties[i].deviceType);
			if (rank < bestRank)
			{
				bestRank = rank;
				chosen = int(i);
			}
		}
		if (chosen < 0)
		{
			if (cpuOnly)
			{
				throw DeviceUnavailable("no Vulkan device apart from the CPU device (run with --device gpu to use it)");
			}
			throw DeviceUnavailable(kNoDevice);
		}
	}
	if (properties[chosen].apiVersion < VK_API_VERSION_1_1)
	{
		throw DeviceUnavailable(std::string("Vulkan device ") + properties[chosen].deviceName + " supports Vulkan "
			+ VersionText(properties[chosen].apiVersion) + " only; version 1.1 is needed (update the graphics driver, or run with --device cpu)");
	}

	_physical = devices[chosen];
	const VkPhysicalDeviceProperties& chosenProperties = properties[chosen];
	_info.name = chosenProperties.deviceName;
	_info.type = TypeName(chosenProperties.deviceType);
	_info.apiVersion = chosenProperties.apiVersion;
	_info.storageRange = chosenProperties.limits.maxStorageBufferRange;
	for (int axis = 0; axis < 3; axis++) _info.workGroupCount[axis] = chosenProperties.limits.maxComputeWorkGroupCount[axis];

	// The allocation limit and the driver name come from the extended
	// properties; the driver name exists from Vulkan 1.2 on.
	VkPhysicalDeviceMaintenance3Properties maintenance{};
	maintenance.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_MAINTENANCE_3_PROPERTIES;
	VkPhysicalDeviceDriverProperties driver{};
	driver.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DRIVER_PROPERTIES;
	VkPhysicalDeviceProperties2 properties2{};
	properties2.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2;
	properties2.pNext = &maintenance;
	if (chosenProperties.apiVersion >= VK_API_VERSION_1_2) maintenance.pNext = &driver;
	vkGetPhysicalDeviceProperties2(_physical, &properties2);
	_info.allocationLimit = maintenance.maxMemoryAllocationSize;
	if (chosenProperties.apiVersion >= VK_API_VERSION_1_2)
	{
		_info.driver = std::string(driver.driverName) + " " + driver.driverInfo;
	}
	else
	{
		_info.driver = "version " + VersionText(chosenProperties.driverVersion);
	}

	VkPhysicalDeviceFeatures features{};
	vkGetPhysicalDeviceFeatures(_physical, &features);
	_info.float64 = features.shaderFloat64 == VK_TRUE;

	vkGetPhysicalDeviceMemoryProperties(_physical, &_memoryProperties);
	for (uint32_t i = 0; i < _memoryProperties.memoryHeapCount; i++)
	{
		const VkMemoryHeap& heap = _memoryProperties.memoryHeaps[i];
		if ((heap.flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT) != 0) _info.heapBytes = std::max(_info.heapBytes, uint64_t(heap.size));
	}
}

void GpuContext::CreateDevice()
{
	// One queue of the first family that runs compute shaders; the same
	// queue also copies buffers.
	uint32_t familyCount = 0;
	vkGetPhysicalDeviceQueueFamilyProperties(_physical, &familyCount, nullptr);
	std::vector<VkQueueFamilyProperties> families(familyCount);
	vkGetPhysicalDeviceQueueFamilyProperties(_physical, &familyCount, families.data());
	bool found = false;
	for (uint32_t i = 0; i < familyCount && !found; i++)
	{
		if ((families[i].queueFlags & VK_QUEUE_COMPUTE_BIT) != 0)
		{
			_queueFamily = i;
			found = true;
		}
	}
	if (!found) throw DeviceUnavailable("Vulkan device " + _info.name + " has no compute queue (run with --device cpu)");

	const float priority = 1.0f;
	VkDeviceQueueCreateInfo queue{};
	queue.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
	queue.queueFamilyIndex = _queueFamily;
	queue.queueCount = 1;
	queue.pQueuePriorities = &priority;

	VkPhysicalDeviceFeatures enabled{};
	enabled.shaderFloat64 = _info.float64 ? VK_TRUE : VK_FALSE;

	VkDeviceCreateInfo info{};
	info.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
	info.queueCreateInfoCount = 1;
	info.pQueueCreateInfos = &queue;
	info.pEnabledFeatures = &enabled;
	Check(vkCreateDevice(_physical, &info, nullptr, &_device), "vkCreateDevice");
	volkLoadDevice(_device);
	vkGetDeviceQueue(_device, _queueFamily, 0, &_queue);

	// One command buffer, recorded and submitted for every transfer or
	// dispatch; the fence waits for it, so the host and the device never
	// work at the same time.
	VkCommandPoolCreateInfo pool{};
	pool.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
	pool.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
	pool.queueFamilyIndex = _queueFamily;
	Check(vkCreateCommandPool(_device, &pool, nullptr, &_commandPool), "vkCreateCommandPool");

	VkCommandBufferAllocateInfo commands{};
	commands.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
	commands.commandPool = _commandPool;
	commands.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
	commands.commandBufferCount = 1;
	Check(vkAllocateCommandBuffers(_device, &commands, &_commands), "vkAllocateCommandBuffers");

	VkFenceCreateInfo fence{};
	fence.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
	Check(vkCreateFence(_device, &fence, nullptr, &_fence), "vkCreateFence");

	// Descriptor sets tell a shader which buffers it reads; one set is
	// written per dispatch and the pool is reset after it ran.
	VkDescriptorPoolSize poolSize{};
	poolSize.type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
	poolSize.descriptorCount = kMaxBindings;
	VkDescriptorPoolCreateInfo descriptors{};
	descriptors.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
	descriptors.maxSets = 1;
	descriptors.poolSizeCount = 1;
	descriptors.pPoolSizes = &poolSize;
	Check(vkCreateDescriptorPool(_device, &descriptors, nullptr, &_descriptorPool), "vkCreateDescriptorPool");
}

GpuPrecision GpuContext::Resolve(GpuPrecision requested) const
{
	if (requested == GpuPrecision::Double && !_info.float64)
	{
		throw std::runtime_error("gpu_precision = double, but the Vulkan device " + _info.name
			+ " has no 64-bit floats; use gpu_precision = single or auto");
	}
	if (requested == GpuPrecision::Auto) return _info.float64 ? GpuPrecision::Double : GpuPrecision::Single;
	return requested;
}

std::string GpuContext::PrecisionReason(GpuPrecision requested) const
{
	if (requested != GpuPrecision::Auto) return "requested";
	return _info.float64 ? "float64 supported" : "no float64 on this device";
}

// ---------------------------------------------------------------------------
// Buffers and transfers
// ---------------------------------------------------------------------------

uint32_t GpuContext::MemoryType(uint32_t typeBits, VkMemoryPropertyFlags required) const
{
	for (uint32_t i = 0; i < _memoryProperties.memoryTypeCount; i++)
	{
		const bool allowed = (typeBits & (1u << i)) != 0;
		const bool suitable = (_memoryProperties.memoryTypes[i].propertyFlags & required) == required;
		if (allowed && suitable) return i;
	}
	throw std::runtime_error("the Vulkan device has no memory type for this buffer");
}

Buffer GpuContext::CreateBuffer(uint64_t bytes)
{
	if (bytes > _info.storageRange || bytes > _info.allocationLimit)
	{
		std::ostringstream message;
		message << "a buffer of " << bytes << " bytes exceeds the Vulkan device limit of "
			<< std::min(_info.storageRange, _info.allocationLimit) << " bytes (storage range " << _info.storageRange
			<< ", allocation " << _info.allocationLimit << "); the tracer must use smaller batches";
		throw std::runtime_error(message.str());
	}
	Buffer buffer;
	buffer._device = _device;
	buffer._bytes = bytes;

	VkBufferCreateInfo info{};
	info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
	info.size = bytes;
	info.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
	info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
	Check(vkCreateBuffer(_device, &info, nullptr, &buffer._buffer), "vkCreateBuffer");

	VkMemoryRequirements requirements{};
	vkGetBufferMemoryRequirements(_device, buffer._buffer, &requirements);
	VkMemoryAllocateInfo allocation{};
	allocation.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
	allocation.allocationSize = requirements.size;
	allocation.memoryTypeIndex = MemoryType(requirements.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
	if (vkAllocateMemory(_device, &allocation, nullptr, &buffer._memory) != VK_SUCCESS)
	{
		std::ostringstream message;
		message << "out of Vulkan device memory for a buffer of " << bytes << " bytes (device heap "
			<< _info.heapBytes << " bytes); reduce the tessellation or the receiver count";
		throw std::runtime_error(message.str());
	}
	Check(vkBindBufferMemory(_device, buffer._buffer, buffer._memory, 0), "vkBindBufferMemory");
	return buffer;
}

void GpuContext::EnsureStaging()
{
	if (_staging != VK_NULL_HANDLE) return;
	VkBufferCreateInfo info{};
	info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
	info.size = kStagingBytes;
	info.usage = VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
	info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
	Check(vkCreateBuffer(_device, &info, nullptr, &_staging), "vkCreateBuffer (staging)");

	VkMemoryRequirements requirements{};
	vkGetBufferMemoryRequirements(_device, _staging, &requirements);
	VkMemoryAllocateInfo allocation{};
	allocation.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
	allocation.allocationSize = requirements.size;
	allocation.memoryTypeIndex = MemoryType(requirements.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
	Check(vkAllocateMemory(_device, &allocation, nullptr, &_stagingMemory), "vkAllocateMemory (staging)");
	Check(vkBindBufferMemory(_device, _staging, _stagingMemory, 0), "vkBindBufferMemory (staging)");
	Check(vkMapMemory(_device, _stagingMemory, 0, kStagingBytes, 0, &_stagingMapped), "vkMapMemory (staging)");
}

void GpuContext::Upload(Buffer& buffer, const void* data, uint64_t bytes, uint64_t offset)
{
	if (offset + bytes > buffer._bytes) throw std::logic_error("Upload past the end of a buffer");
	EnsureStaging();
	const char* source = static_cast<const char*>(data);
	for (uint64_t done = 0; done < bytes; )
	{
		const uint64_t chunk = std::min(kStagingBytes, bytes - done);
		std::memcpy(_stagingMapped, source + done, size_t(chunk));
		Begin();
		VkBufferCopy region{};
		region.srcOffset = 0;
		region.dstOffset = offset + done;
		region.size = chunk;
		vkCmdCopyBuffer(_commands, _staging, buffer._buffer, 1, &region);
		Barrier();
		Submit();
		done += chunk;
	}
}

void GpuContext::Download(const Buffer& buffer, void* data, uint64_t bytes, uint64_t offset)
{
	if (offset + bytes > buffer._bytes) throw std::logic_error("Download past the end of a buffer");
	EnsureStaging();
	char* destination = static_cast<char*>(data);
	for (uint64_t done = 0; done < bytes; )
	{
		const uint64_t chunk = std::min(kStagingBytes, bytes - done);
		Begin();
		VkBufferCopy region{};
		region.srcOffset = offset + done;
		region.dstOffset = 0;
		region.size = chunk;
		vkCmdCopyBuffer(_commands, buffer._buffer, _staging, 1, &region);
		Barrier();
		Submit();
		std::memcpy(destination + done, _stagingMapped, size_t(chunk));
		done += chunk;
	}
}

void GpuContext::Fill(Buffer& buffer, uint32_t value)
{
	Begin();
	vkCmdFillBuffer(_commands, buffer._buffer, 0, VK_WHOLE_SIZE, value);
	Barrier();
	Submit();
}

// ---------------------------------------------------------------------------
// Pipelines and dispatches
// ---------------------------------------------------------------------------

Pipeline GpuContext::CreateCompute(const uint32_t* spirv, size_t words, uint32_t bindings, uint32_t pushBytes)
{
	if (bindings > kMaxBindings) throw std::logic_error("a pipeline may bind at most 16 buffers");
	Pipeline pipeline;
	pipeline._device = _device;
	pipeline._bindings = bindings;
	pipeline._pushBytes = pushBytes;

	// The descriptor set layout lists the storage buffers the shader binds,
	// the pipeline layout adds the push constants.
	std::vector<VkDescriptorSetLayoutBinding> layoutBindings(bindings);
	for (uint32_t i = 0; i < bindings; i++)
	{
		layoutBindings[i].binding = i;
		layoutBindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
		layoutBindings[i].descriptorCount = 1;
		layoutBindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
	}
	VkDescriptorSetLayoutCreateInfo setLayout{};
	setLayout.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
	setLayout.bindingCount = bindings;
	setLayout.pBindings = layoutBindings.data();
	Check(vkCreateDescriptorSetLayout(_device, &setLayout, nullptr, &pipeline._setLayout), "vkCreateDescriptorSetLayout");

	VkPushConstantRange pushRange{};
	pushRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
	pushRange.offset = 0;
	pushRange.size = pushBytes;
	VkPipelineLayoutCreateInfo layout{};
	layout.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
	layout.setLayoutCount = 1;
	layout.pSetLayouts = &pipeline._setLayout;
	layout.pushConstantRangeCount = pushBytes > 0 ? 1 : 0;
	layout.pPushConstantRanges = &pushRange;
	Check(vkCreatePipelineLayout(_device, &layout, nullptr, &pipeline._layout), "vkCreatePipelineLayout");

	VkShaderModuleCreateInfo module{};
	module.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
	module.codeSize = words * sizeof(uint32_t);
	module.pCode = spirv;
	VkShaderModule shader = VK_NULL_HANDLE;
	Check(vkCreateShaderModule(_device, &module, nullptr, &shader), "vkCreateShaderModule");

	VkComputePipelineCreateInfo info{};
	info.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
	info.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
	info.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
	info.stage.module = shader;
	info.stage.pName = "main";
	info.layout = pipeline._layout;
	const VkResult result = vkCreateComputePipelines(_device, VK_NULL_HANDLE, 1, &info, nullptr, &pipeline._pipeline);
	vkDestroyShaderModule(_device, shader, nullptr);
	Check(result, "vkCreateComputePipelines");
	return pipeline;
}

void GpuContext::Dispatch(const Pipeline& pipeline, const std::vector<const Buffer*>& buffers, void* push, uint32_t pushBytes, uint64_t itemsX, uint64_t itemsY)
{
	if (buffers.size() != pipeline._bindings) throw std::logic_error("Dispatch got a different number of buffers than the pipeline binds");
	if (pushBytes != pipeline._pushBytes || pushBytes < sizeof(DispatchBase)) throw std::logic_error("Dispatch got push constants of the wrong size");
	if (itemsX > std::numeric_limits<uint32_t>::max() || itemsY > std::numeric_limits<uint32_t>::max())
	{
		throw std::logic_error("Dispatch over more than 2^32 items in one dimension");
	}

	// A device accepts at most maxComputeWorkGroupCount groups per
	// dimension, so a large dispatch runs as a loop of chunks. The shader
	// adds the chunk's first item to its invocation id.
	const uint64_t maxX = uint64_t(_info.workGroupCount[0]) * kWorkGroupSize;
	const uint64_t maxY = _info.workGroupCount[1];
	for (uint64_t firstY = 0; firstY < itemsY; firstY += maxY)
	{
		const uint64_t countY = std::min(maxY, itemsY - firstY);
		for (uint64_t firstX = 0; firstX < itemsX; firstX += maxX)
		{
			const uint64_t countX = std::min(maxX, itemsX - firstX);
			DispatchBase* base = static_cast<DispatchBase*>(push);
			base->firstX = uint32_t(firstX);
			base->firstY = uint32_t(firstY);

			Check(vkResetDescriptorPool(_device, _descriptorPool, 0), "vkResetDescriptorPool");
			VkDescriptorSetAllocateInfo allocation{};
			allocation.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
			allocation.descriptorPool = _descriptorPool;
			allocation.descriptorSetCount = 1;
			allocation.pSetLayouts = &pipeline._setLayout;
			VkDescriptorSet set = VK_NULL_HANDLE;
			Check(vkAllocateDescriptorSets(_device, &allocation, &set), "vkAllocateDescriptorSets");
			std::vector<VkDescriptorBufferInfo> bufferInfos(buffers.size());
			std::vector<VkWriteDescriptorSet> writes(buffers.size());
			for (size_t i = 0; i < buffers.size(); i++)
			{
				bufferInfos[i].buffer = buffers[i]->_buffer;
				bufferInfos[i].offset = 0;
				bufferInfos[i].range = VK_WHOLE_SIZE;
				writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
				writes[i].dstSet = set;
				writes[i].dstBinding = uint32_t(i);
				writes[i].descriptorCount = 1;
				writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
				writes[i].pBufferInfo = &bufferInfos[i];
			}
			if (!writes.empty()) vkUpdateDescriptorSets(_device, uint32_t(writes.size()), writes.data(), 0, nullptr);

			Begin();
			vkCmdBindPipeline(_commands, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline._pipeline);
			vkCmdBindDescriptorSets(_commands, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline._layout, 0, 1, &set, 0, nullptr);
			vkCmdPushConstants(_commands, pipeline._layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, pushBytes, push);
			vkCmdDispatch(_commands, uint32_t((countX + kWorkGroupSize - 1) / kWorkGroupSize), uint32_t(countY), 1);
			Barrier();
			Submit();
		}
	}
}

// ---------------------------------------------------------------------------
// Command buffer helpers
// ---------------------------------------------------------------------------

void GpuContext::Begin()
{
	VkCommandBufferBeginInfo info{};
	info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
	info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
	Check(vkBeginCommandBuffer(_commands, &info), "vkBeginCommandBuffer");
}

void GpuContext::Barrier()
{
	// Every command buffer ends with this barrier, so what one dispatch or
	// copy wrote is visible to the next one and to the host.
	VkMemoryBarrier barrier{};
	barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
	barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_TRANSFER_WRITE_BIT;
	barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_TRANSFER_READ_BIT
		| VK_ACCESS_TRANSFER_WRITE_BIT | VK_ACCESS_HOST_READ_BIT;
	vkCmdPipelineBarrier(_commands, VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, VK_PIPELINE_STAGE_ALL_COMMANDS_BIT | VK_PIPELINE_STAGE_HOST_BIT,
		0, 1, &barrier, 0, nullptr, 0, nullptr);
}

void GpuContext::Submit()
{
	Check(vkEndCommandBuffer(_commands), "vkEndCommandBuffer");
	VkSubmitInfo submit{};
	submit.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
	submit.commandBufferCount = 1;
	submit.pCommandBuffers = &_commands;
	Check(vkQueueSubmit(_queue, 1, &submit, _fence), "vkQueueSubmit");
	Check(vkWaitForFences(_device, 1, &_fence, VK_TRUE, std::numeric_limits<uint64_t>::max()), "vkWaitForFences");
	Check(vkResetFences(_device, 1, &_fence), "vkResetFences");
}

} // namespace gpu
} // namespace sbr
