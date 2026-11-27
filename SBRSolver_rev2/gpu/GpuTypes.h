#pragma once

// Which tracer runs, and in which precision the GPU shaders compute. The
// configuration file and the command line set them; the GPU context
// resolves `Auto`.

enum class TracerDevice
{
	Cpu,   // the reference tracer
	Gpu,   // the Vulkan tracer; the run stops with exit code 77 without a device
	Auto   // the Vulkan tracer when a GPU exists, otherwise the CPU tracer
};

enum class GpuPrecision
{
	Auto,    // double when the device computes in 64-bit floats, else single
	Single,
	Double
};

inline const char* TracerDeviceName(TracerDevice device)
{
	switch (device)
	{
	case TracerDevice::Cpu: return "cpu";
	case TracerDevice::Gpu: return "gpu";
	case TracerDevice::Auto: return "auto";
	}
	return "unknown";
}

inline const char* GpuPrecisionName(GpuPrecision precision)
{
	switch (precision)
	{
	case GpuPrecision::Auto: return "auto";
	case GpuPrecision::Single: return "single";
	case GpuPrecision::Double: return "double";
	}
	return "unknown";
}
