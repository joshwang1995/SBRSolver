#include "SimulationConfig.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

namespace fs = std::filesystem;

namespace
{
	// Sanity limits. The receiver limit keeps the generators' int
	// conversions in range and the per-receiver vectors in memory.
	const int kMaxTessellation = 8;                      // 655,362 launch directions
	const int kMaxThreads = 1024;
	const double kMaxReceivers = 100000000.0;

	int LaunchRayCount(int tessellation)
	{
		return 10 * (1 << (2 * tessellation)) + 2;      // 12, 42, 162, 642, ...
	}

	struct Entry
	{
		std::string value;
		int line = 0;
		bool used = false;
	};

	struct Section
	{
		std::string name;
		int line = 0;
		std::map<std::string, Entry> entries;
	};

	std::string Trim(const std::string& text)
	{
		const char* whitespace = " \t\r\n";
		const std::string::size_type first = text.find_first_not_of(whitespace);
		if (first == std::string::npos) return std::string();
		const std::string::size_type last = text.find_last_not_of(whitespace);
		return text.substr(first, last - first + 1);
	}

	std::string ToLower(std::string text)
	{
		std::transform(text.begin(), text.end(), text.begin(),
			[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		return text;
	}

	// Collects "file:line: message" diagnostics.
	class Diagnostics
	{
	public:
		explicit Diagnostics(std::string file) : _file(std::move(file)) {}

		void Add(int line, const std::string& message)
		{
			std::ostringstream oss;
			oss << _file << ":" << line << ": " << message;
			_messages.push_back(oss.str());
		}

		bool Empty() const { return _messages.empty(); }

		std::string Join() const
		{
			std::string joined;
			for (const std::string& m : _messages) joined += m + "\n";
			return joined;
		}

	private:
		std::string _file;
		std::vector<std::string> _messages;
	};

	bool ParseIniFile(const std::string& path, std::vector<Section>& sections, Diagnostics& diag)
	{
		// Only a directory or an unreadable path is refused here; a pipe or
		// process substitution falls through to std::ifstream as before.
		std::error_code ec;
		const fs::file_status status = fs::status(fs::path(path), ec);
		if (status.type() == fs::file_type::not_found)
		{
			diag.Add(0, "configuration file not found");
			return false;
		}
		if (ec)
		{
			diag.Add(0, "cannot open configuration file: " + ec.message());
			return false;
		}
		if (fs::is_directory(status))
		{
			diag.Add(0, "configuration file is a directory");
			return false;
		}
		std::ifstream in(path);
		if (!in)
		{
			diag.Add(0, "cannot open configuration file");
			return false;
		}

		std::string rawLine;
		int lineNo = 0;
		Section* current = nullptr;
		while (std::getline(in, rawLine))
		{
			lineNo++;
			// Drop a UTF-8 byte order mark; some Windows tools write one (Visual
			// Studio "UTF-8 with signature", Notepad "UTF-8 with BOM",
			// Windows PowerShell 5.1 Out-File -Encoding UTF8).
			if (lineNo == 1 && rawLine.compare(0, 3, "\xEF\xBB\xBF") == 0) rawLine.erase(0, 3);
			std::string line = Trim(rawLine);
			// Whole-line comments start with '#' or ';'. Inline comments start
			// with '#' preceded by whitespace, so ';' stays usable as a row
			// separator inside values.
			if (line.empty() || line.front() == '#' || line.front() == ';') continue;
			for (std::string::size_type i = 1; i < line.size(); i++)
			{
				if (line[i] == '#' && std::isspace(static_cast<unsigned char>(line[i - 1])))
				{
					line.erase(i);
					break;
				}
			}
			line = Trim(line);
			if (line.empty()) continue;

			if (line.front() == '[')
			{
				if (line.back() != ']')
				{
					diag.Add(lineNo, "malformed section header '" + line + "'");
					continue;
				}
				Section section;
				section.name = ToLower(Trim(line.substr(1, line.size() - 2)));
				section.line = lineNo;
				if (section.name.empty())
				{
					diag.Add(lineNo, "empty section name");
					continue;
				}
				sections.push_back(section);
				current = &sections.back();
				continue;
			}

			const std::string::size_type eq = line.find('=');
			if (eq == std::string::npos)
			{
				diag.Add(lineNo, "expected 'key = value', got '" + line + "'");
				continue;
			}
			if (current == nullptr)
			{
				diag.Add(lineNo, "key outside of any [section]");
				continue;
			}
			const std::string key = ToLower(Trim(line.substr(0, eq)));
			const std::string value = Trim(line.substr(eq + 1));
			if (key.empty())
			{
				diag.Add(lineNo, "empty key");
				continue;
			}
			if (current->entries.count(key))
			{
				diag.Add(lineNo, "duplicate key '" + key + "' in section [" + current->name + "]");
				continue;
			}
			current->entries[key] = Entry{ value, lineNo, false };
		}
		return true;
	}

	// Numbers are separated by whitespace (';' additionally separates matrix
	// rows); the decimal mark is always '.'. Every token must be a complete
	// number, so '1.8e9-' or '2,5' are rejected rather than silently truncated.
	bool ParseDoubles(const std::string& text, std::vector<double>& out)
	{
		std::string cleaned = text;
		std::replace(cleaned.begin(), cleaned.end(), ';', ' ');
		std::istringstream iss(cleaned);
		std::string token;
		while (iss >> token)
		{
			std::istringstream t(token);
			double v;
			if (!(t >> v) || !t.eof() || !std::isfinite(v)) return false;
			out.push_back(v);
		}
		return !out.empty();
	}

	bool ParseBoolean(const std::string& text, bool& out)
	{
		const std::string v = ToLower(text);
		if (v == "true" || v == "yes" || v == "on" || v == "1") { out = true; return true; }
		if (v == "false" || v == "no" || v == "off" || v == "0") { out = false; return true; }
		return false;
	}

	// Typed access to one section with "every key must be consumed" checking.
	class SectionReader
	{
	public:
		SectionReader(Section& section, Diagnostics& diag, const fs::path& baseDir)
			: _section(section), _diag(diag), _baseDir(baseDir) {}

		bool Has(const std::string& key) const { return _section.entries.count(key) != 0; }

		int Line(const std::string& key) const
		{
			auto it = _section.entries.find(key);
			return it == _section.entries.end() ? _section.line : it->second.line;
		}

		bool String(const std::string& key, std::string& out, bool required)
		{
			auto it = _section.entries.find(key);
			if (it == _section.entries.end())
			{
				if (required) _diag.Add(_section.line, "missing required key '" + key + "' in section [" + _section.name + "]");
				return false;
			}
			it->second.used = true;
			if (it->second.value.empty())
			{
				_diag.Add(it->second.line, "key '" + key + "' has no value");
				return false;
			}
			out = it->second.value;
			return true;
		}

		bool Double(const std::string& key, double& out, bool required)
		{
			std::string text;
			if (!String(key, text, required)) return false;
			std::vector<double> values;
			if (!ParseDoubles(text, values) || values.size() != 1)
			{
				_diag.Add(Line(key), "key '" + key + "' expects one number (decimal mark '.'), got '" + text + "'");
				return false;
			}
			out = values[0];
			return true;
		}

		bool Integer(const std::string& key, int& out, bool required)
		{
			std::string text;
			if (!String(key, text, required)) return false;
			std::istringstream iss(text);
			int value;
			// operator>>(int&) sets failbit for values outside the int range.
			if (!(iss >> value) || !(iss >> std::ws).eof())
			{
				_diag.Add(Line(key), "key '" + key + "' expects an integer within the int range, got '" + text + "'");
				return false;
			}
			out = value;
			return true;
		}

		bool Boolean(const std::string& key, bool& out, bool required)
		{
			std::string text;
			if (!String(key, text, required)) return false;
			if (!ParseBoolean(text, out))
			{
				_diag.Add(Line(key), "key '" + key + "' expects true or false, got '" + text + "'");
				return false;
			}
			return true;
		}

		bool Vector3(const std::string& key, Vec3& out, bool required)
		{
			std::string text;
			if (!String(key, text, required)) return false;
			std::vector<double> values;
			if (!ParseDoubles(text, values) || values.size() != 3)
			{
				_diag.Add(Line(key), "key '" + key + "' expects three space-separated numbers 'x y z' (decimal mark '.'), got '" + text + "'");
				return false;
			}
			out = Vec3(values[0], values[1], values[2]);
			return true;
		}

		// Nine numbers, row by row; rows may be separated by ';'.
		bool Matrix3(const std::string& key, Mat3& out, bool required)
		{
			std::string text;
			if (!String(key, text, required)) return false;
			std::vector<double> values;
			if (!ParseDoubles(text, values) || values.size() != 9)
			{
				_diag.Add(Line(key), "key '" + key + "' expects nine space-separated numbers (three rows of three, rows may be separated by ';'), got '" + text + "'");
				return false;
			}
			for (int r = 0; r < 3; r++)
				for (int c = 0; c < 3; c++)
					out(r, c) = values[3 * r + c];
			return true;
		}

		// A file name, resolved against the configuration file's directory.
		bool Path(const std::string& key, std::string& out, bool required, bool mustExist)
		{
			std::string text;
			if (!String(key, text, required)) return false;
			fs::path p(text);
			if (p.is_relative()) p = _baseDir / p;
			p = p.lexically_normal();
			if (mustExist)
			{
				// Non-throwing overload: an unreadable path (ENAMETOOLONG, ELOOP,
				// EACCES, ...) must become a diagnostic, not a filesystem_error.
				std::error_code ec;
				const fs::file_status status = fs::status(p, ec);
				if (!fs::is_regular_file(status))
				{
					std::string why;
					if (status.type() == fs::file_type::not_found) why = "file not found";
					else if (ec) why = ec.message();
					else if (fs::is_directory(status)) why = "is a directory, not a file";
					else why = "is not a regular file";
					_diag.Add(Line(key), "key '" + key + "': " + why + ": " + p.generic_string());
					return false;
				}
			}
			out = p.generic_string();
			return true;
		}

		std::vector<std::string> Keys() const
		{
			std::vector<std::string> keys;
			for (const auto& kv : _section.entries) keys.push_back(kv.first);
			return keys;
		}

		void ReportUnusedKeys()
		{
			for (const auto& kv : _section.entries)
			{
				if (!kv.second.used)
					_diag.Add(kv.second.line, "unknown key '" + kv.first + "' in section [" + _section.name + "]");
			}
		}

	private:
		Section& _section;
		Diagnostics& _diag;
		fs::path _baseDir;
	};

	bool ParseAntennaType(const std::string& text, AntennaType& out)
	{
		const std::string v = ToLower(text);
		if (v == "isotropic") { out = AntennaType::Isotropic; return true; }
		if (v == "hertzian_dipole") { out = AntennaType::HertzianDipole; return true; }
		if (v == "half_wave_dipole") { out = AntennaType::HalfWaveDipole; return true; }
		return false;
	}

	bool ParseReceiverLayout(const std::string& text, ReceiverLayout& out)
	{
		const std::string v = ToLower(text);
		if (v == "line") { out = ReceiverLayout::Line; return true; }
		if (v == "plane") { out = ReceiverLayout::Plane; return true; }
		if (v == "file") { out = ReceiverLayout::File; return true; }
		return false;
	}

	int FindMaterial(const SimulationConfig& cfg, const std::string& name)
	{
		for (size_t i = 0; i < cfg.materials.size(); i++)
		{
			if (cfg.materials[i].name == name) return int(i);
		}
		return -1;
	}

	// Resolves a material name from `key`; reports unknown names.
	bool ReadMaterialReference(SectionReader& r, const SimulationConfig& cfg, const std::string& key,
		bool required, std::string& name, int& id, Diagnostics& diag)
	{
		if (!r.String(key, name, required)) return false;
		id = FindMaterial(cfg, name);
		if (id < 0)
		{
			std::string known;
			for (const MaterialConfig& m : cfg.materials) known += (known.empty() ? "" : ", ") + m.name;
			diag.Add(r.Line(key), "key '" + key + "': unknown material '" + name + "' (defined materials: " + known + ")");
			return false;
		}
		return true;
	}

	void ReadScene(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		r.Path("stl_file", cfg.stlFile, true, true);
		if (!cfg.materials.empty())
		{
			cfg.sceneMaterial = cfg.materials[0].name;
			cfg.sceneMaterialId = 0;
		}
		ReadMaterialReference(r, cfg, "material", false, cfg.sceneMaterial, cfg.sceneMaterialId, diag);
		r.Boolean("save_edges", cfg.saveEdges, false);
		r.Boolean("save_faces", cfg.saveFaces, false);
	}

	void ReadGeometry(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		GeometryConfig g;
		r.Path("file", g.file, true, true);
		ReadMaterialReference(r, cfg, "material", true, g.material, g.materialId, diag);
		cfg.geometries.push_back(g);
	}

	// Every key is a coplanar surface id, every value a material name.
	void ReadSurfaceMaterials(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		for (const std::string& key : r.Keys())
		{
			// Plain decimal digits only (no sign, no leading zeros), so every
			// surface id has exactly one spelling and the parser's duplicate-key
			// check catches an id given twice.
			const bool canonical = !key.empty()
				&& std::all_of(key.begin(), key.end(), [](unsigned char c) { return std::isdigit(c) != 0; })
				&& (key.size() == 1 || key[0] != '0');
			std::istringstream iss(key);
			int surfaceId;
			if (!canonical || !(iss >> surfaceId) || !(iss >> std::ws).eof())
			{
				std::string ignored;
				r.String(key, ignored, false); // mark as used so it is not reported twice
				diag.Add(r.Line(key), "[surface_materials] keys must be plain decimal surface ids such as 12 (see the CoplanarID field of Faces.vtk), got '" + key + "'");
				continue;
			}
			std::string name;
			int id;
			if (ReadMaterialReference(r, cfg, key, true, name, id, diag))
			{
				cfg.surfaceMaterials[surfaceId] = id;
			}
		}
	}

	// A transmitter name labels its paths in the viewer and the log, so it is
	// one word: letters, digits, '_', '-' and '.'.
	bool ValidTransmitterName(const std::string& name)
	{
		if (name.empty()) return false;
		return std::all_of(name.begin(), name.end(),
			[](unsigned char c) { return std::isalnum(c) != 0 || c == '_' || c == '-' || c == '.'; });
	}

	// One [transmitter] section. A later section starts from the first
	// one's values, so a second access point needs only its position. The
	// frequency is one for the run, because the fields of all transmitters
	// are summed at every receiver.
	void ReadTransmitter(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		const bool first = cfg.transmitters.empty();
		TransmitterConfig tx = first ? TransmitterConfig() : cfg.transmitters.front();
		tx.position = Vec3::Zero();
		const auto nameUsed = [&cfg](const std::string& name)
		{
			return std::any_of(cfg.transmitters.begin(), cfg.transmitters.end(),
				[&name](const TransmitterConfig& other) { return other.name == name; });
		};
		if (r.String("name", tx.name, false))
		{
			if (!ValidTransmitterName(tx.name))
				diag.Add(r.Line("name"), "name must be one word of letters, digits, '_', '-' and '.'; got '" + tx.name + "'");
			else if (nameUsed(tx.name))
				diag.Add(r.Line("name"), "transmitter name '" + tx.name + "' is used more than once");
		}
		else
		{
			// The first unused default, so it never clashes with a given name.
			for (int n = 0; ; n++)
			{
				tx.name = "tx" + std::to_string(n);
				if (!nameUsed(tx.name)) break;
			}
		}
		r.Vector3("position", tx.position, true);
		if (r.Matrix3("coordinate_system", tx.coordSys, false))
		{
			// The rows are the antenna axes in global coordinates and are used
			// as such (no normalisation), so they must form a right-handed
			// orthonormal frame.
			const Mat3 gram = tx.coordSys * tx.coordSys.transpose();
			if ((gram - Mat3::Identity()).cwiseAbs().maxCoeff() > 1e-9 || tx.coordSys.determinant() < 0.0)
				diag.Add(r.Line("coordinate_system"), "coordinate_system rows must be unit length, mutually perpendicular and right-handed (x, y, z axes of the antenna in global coordinates)");
		}
		if (r.Double("power_watts", tx.powerWatts, false) && tx.powerWatts <= 0.0)
			diag.Add(r.Line("power_watts"), "power_watts must be positive");
		double frequency = cfg.frequencyHz;
		if (r.Double("frequency_hz", frequency, first))
		{
			if (frequency <= 0.0)
				diag.Add(r.Line("frequency_hz"), "frequency_hz must be positive");
			else if (first)
				cfg.frequencyHz = frequency;
			else if (cfg.frequencyHz > 0.0 && frequency != cfg.frequencyHz)
			{
				// A bad first frequency is reported there, so compare only with a valid one.
				std::ostringstream firstValue;
				firstValue << cfg.frequencyHz;
				diag.Add(r.Line("frequency_hz"), "frequency_hz must be the same for every transmitter, because their fields are summed at the receivers; the first transmitter has " + firstValue.str());
			}
		}
		if (r.Integer("tessellation", tx.tessellation, false) && (tx.tessellation < 0 || tx.tessellation > kMaxTessellation))
			diag.Add(r.Line("tessellation"), "tessellation must be between 0 and " + std::to_string(kMaxTessellation)
				+ " (level " + std::to_string(kMaxTessellation) + " already launches " + std::to_string(LaunchRayCount(kMaxTessellation)) + " rays)");
		std::string antenna;
		if (r.String("antenna_type", antenna, false) && !ParseAntennaType(antenna, tx.antennaType))
			diag.Add(r.Line("antenna_type"), "antenna_type must be one of isotropic, hertzian_dipole, half_wave_dipole; got '" + antenna + "'");
		r.Path("pattern_file", tx.patternFile, false, true);
		cfg.transmitters.push_back(tx);
	}

	// Checks that a generated line or plane has a sane, representable number
	// of points before Preprocessor::GenerateRXLine/GenerateRxPlane convert the
	// counts to int. Mirrors the generators: a line yields
	// floor(length / resolution) points, a plane (nx + 1) * (ny + 1) with
	// nx = floor(dx / resolution) and ny = floor(dy / resolution). Only called
	// when every value it uses parsed and passed its own range check.
	void ValidateReceiverCount(SectionReader& r, const ReceiverConfig& rx, Diagnostics& diag)
	{
		double count = 0.0;
		if (rx.layout == ReceiverLayout::Line)
		{
			const double length = (rx.lineEnd - rx.lineStart).norm();
			if (!(length > 0.0))
			{
				diag.Add(r.Line("line_end"), "line_start and line_end must be different points");
				return;
			}
			count = std::floor(length / rx.resolution);
			if (count < 1.0)
			{
				diag.Add(r.Line("resolution"), "the line is shorter than resolution, so no receivers would be generated");
				return;
			}
		}
		else if (rx.layout == ReceiverLayout::Plane)
		{
			const double nx = std::floor((rx.xMax - rx.xMin) / rx.resolution);
			const double ny = std::floor((rx.yMax - rx.yMin) / rx.resolution);
			count = (nx + 1.0) * (ny + 1.0);
		}
		else
		{
			return;
		}
		if (!std::isfinite(count) || count > kMaxReceivers)
		{
			std::ostringstream msg;
			msg << "the receiver layout would generate " << std::fixed << std::setprecision(0) << count
				<< " receivers, more than the limit of " << kMaxReceivers << "; increase resolution or shrink the extent";
			diag.Add(r.Line("resolution"), msg.str());
		}
	}

	// Reads `resolution`; true only when it parsed and is positive.
	bool ReadResolution(SectionReader& r, ReceiverConfig& rx, Diagnostics& diag)
	{
		if (!r.Double("resolution", rx.resolution, true)) return false;
		if (rx.resolution <= 0.0)
		{
			diag.Add(r.Line("resolution"), "resolution must be positive");
			return false;
		}
		return true;
	}

	void ReadReceivers(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		ReceiverConfig& rx = cfg.receivers;
		std::string layout;
		if (!r.String("layout", layout, true)) return;
		if (!ParseReceiverLayout(layout, rx.layout))
		{
			diag.Add(r.Line("layout"), "layout must be one of line, plane, file; got '" + layout + "'");
			return;
		}
		// `ok` tracks only the values the count check derives from, so an
		// error in another section cannot hide a receiver-count problem.
		// (`&=` rather than `&&=` so every key is still read and reported.)
		bool ok = true;
		switch (rx.layout)
		{
		case ReceiverLayout::Line:
			ok &= r.Vector3("line_start", rx.lineStart, true);
			ok &= r.Vector3("line_end", rx.lineEnd, true);
			ok &= ReadResolution(r, rx, diag);
			break;
		case ReceiverLayout::Plane:
			ok &= r.Double("x_min", rx.xMin, true);
			ok &= r.Double("y_min", rx.yMin, true);
			ok &= r.Double("x_max", rx.xMax, true);
			ok &= r.Double("y_max", rx.yMax, true);
			r.Double("height", rx.height, true);
			ok &= ReadResolution(r, rx, diag);
			if (ok && (rx.xMax < rx.xMin || rx.yMax < rx.yMin))
			{
				diag.Add(r.Line("x_max"), "plane bounds must satisfy x_min <= x_max and y_min <= y_max");
				ok = false;
			}
			break;
		case ReceiverLayout::File:
			r.Path("file", rx.file, true, true);
			ok = false;
			break;
		}
		if (ok) ValidateReceiverCount(r, rx, diag);
	}

	void ReadRayTracing(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		if (r.Integer("max_reflections", cfg.maxReflections, false) && cfg.maxReflections < 0)
			diag.Add(r.Line("max_reflections"), "max_reflections must be 0 or greater");
		if (r.Integer("max_transmissions", cfg.maxTransmissions, false) && cfg.maxTransmissions < 0)
			diag.Add(r.Line("max_transmissions"), "max_transmissions must be 0 or greater");
		r.Boolean("use_fresnel_coefficients", cfg.useFresnelCoefficients, false);
		std::string sum;
		if (r.String("transmitter_sum", sum, false))
		{
			const std::string v = ToLower(sum);
			if (v == "coherent") cfg.powerSum = false;
			else if (v == "power") cfg.powerSum = true;
			else diag.Add(r.Line("transmitter_sum"), "transmitter_sum must be coherent or power; got '" + sum + "'");
		}
		if (r.Integer("threads", cfg.threads, false) && (cfg.threads < 0 || cfg.threads > kMaxThreads))
			diag.Add(r.Line("threads"), "threads must be between 0 (all hardware threads) and " + std::to_string(kMaxThreads));
		std::string device;
		if (r.String("device", device, false) && !ParseTracerDevice(device, cfg.device))
			diag.Add(r.Line("device"), "device must be cpu, gpu or auto; got '" + device + "'");
		std::string precision;
		if (r.String("gpu_precision", precision, false) && !ParseGpuPrecision(precision, cfg.gpuPrecision))
			diag.Add(r.Line("gpu_precision"), "gpu_precision must be auto, single or double; got '" + precision + "'");
	}

	void ReadMaterial(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		MaterialConfig m;
		MaterialProperties& p = m.properties;
		if (r.String("name", m.name, true) && FindMaterial(cfg, m.name) >= 0)
			diag.Add(r.Line("name"), "material name '" + m.name + "' is used more than once");
		// Defaults to the transmitter frequency when omitted (the solver only
		// evaluates materials at the simulation frequency).
		p.frequency = cfg.frequencyHz;
		r.Double("frequency_hz", p.frequency, false);
		r.Double("reflection_loss", p.reflectionLoss, true);
		r.Double("transmission_loss", p.transmissionLoss, true);
		if (r.Double("rel_conductivity", p.relConductivity, true) && p.relConductivity < 0.0)
			diag.Add(r.Line("rel_conductivity"), "rel_conductivity must be 0 or greater");
		if (r.Double("rel_permittivity_re", p.relPermittivityRe, true) && p.relPermittivityRe <= 0.0)
			diag.Add(r.Line("rel_permittivity_re"), "rel_permittivity_re must be positive");
		r.Double("rel_permittivity_im", p.relPermittivityIm, true);
		cfg.materials.push_back(m);
	}

	void ReadOutputs(SectionReader& r, SimulationConfig& cfg)
	{
		r.Boolean("ray_paths", cfg.writeRayPaths, false);
		r.Boolean("receivers", cfg.writeReceivers, false);
		r.Boolean("field_csv", cfg.writeFieldCsv, false);
		r.Boolean("icosahedron", cfg.writeIcosahedron, false);
		r.Boolean("paths_csv", cfg.writePathsCsv, false);
		r.Boolean("viewer", cfg.viewer, false);
	}
} // namespace

bool LoadSimulationConfig(const std::string& path, SimulationConfig& config, std::string& error)
{
	error.clear();
	config = SimulationConfig();

	std::error_code ec;
	fs::path absolute = fs::absolute(fs::path(path), ec);
	if (ec) absolute = fs::path(path);
	absolute = absolute.lexically_normal();
	config.configFile = absolute.generic_string();
	const fs::path baseDir = absolute.parent_path();

	Diagnostics diag(config.configFile);
	std::vector<Section> sections;
	if (!ParseIniFile(config.configFile, sections, diag))
	{
		error = diag.Join();
		return false;
	}

	// Sections are read in a fixed order regardless of file order: [material]
	// needs the transmitter frequency for its default, and [scene],
	// [geometry] and [surface_materials] refer to materials by name.
	// Repeated [transmitter] sections keep their file order.
	const char* order[] = { "transmitter", "material", "scene", "geometry", "receivers", "ray_tracing", "surface_materials", "outputs" };
	const int sectionCount = int(sizeof(order) / sizeof(order[0]));
	std::vector<bool> seen(sectionCount, false);
	for (int i = 0; i < sectionCount; i++)
	{
		const std::string name = order[i];
		const bool repeatable = (name == "material" || name == "geometry" || name == "transmitter");
		for (Section& section : sections)
		{
			if (section.name != name) continue;
			if (seen[i] && !repeatable)
			{
				diag.Add(section.line, "section [" + name + "] appears more than once");
				continue;
			}
			seen[i] = true;
			SectionReader reader(section, diag, baseDir);
			if (name == "scene") ReadScene(reader, config, diag);
			else if (name == "transmitter") ReadTransmitter(reader, config, diag);
			else if (name == "receivers") ReadReceivers(reader, config, diag);
			else if (name == "ray_tracing") ReadRayTracing(reader, config, diag);
			else if (name == "material") ReadMaterial(reader, config, diag);
			else if (name == "geometry") ReadGeometry(reader, config, diag);
			else if (name == "surface_materials") ReadSurfaceMaterials(reader, config, diag);
			else ReadOutputs(reader, config);
			reader.ReportUnusedKeys();
		}
	}
	for (const Section& section : sections)
	{
		if (std::find(std::begin(order), std::end(order), section.name) == std::end(order))
			diag.Add(section.line, "unknown section [" + section.name + "]");
	}
	const char* required[] = { "scene", "transmitter", "receivers" };
	for (const char* name : required)
	{
		if (!std::any_of(sections.begin(), sections.end(), [&](const Section& s) { return s.name == name; }))
			diag.Add(0, std::string("missing required section [") + name + "]");
	}
	if (config.materials.empty())
		diag.Add(0, "at least one [material] section is required (the first one is the default for [scene] stl_file)");

	error = diag.Join();
	return diag.Empty();
}

const char* AntennaTypeName(AntennaType type)
{
	switch (type)
	{
	case AntennaType::Isotropic: return "isotropic";
	case AntennaType::HertzianDipole: return "hertzian_dipole";
	case AntennaType::HalfWaveDipole: return "half_wave_dipole";
	}
	return "unknown";
}

const char* ReceiverLayoutName(ReceiverLayout layout)
{
	switch (layout)
	{
	case ReceiverLayout::Line: return "line";
	case ReceiverLayout::Plane: return "plane";
	case ReceiverLayout::File: return "file";
	}
	return "unknown";
}

bool ParseTracerDevice(const std::string& text, TracerDevice& out)
{
	const std::string v = ToLower(text);
	if (v == "cpu") { out = TracerDevice::Cpu; return true; }
	if (v == "gpu") { out = TracerDevice::Gpu; return true; }
	if (v == "auto") { out = TracerDevice::Auto; return true; }
	return false;
}

bool ParseGpuPrecision(const std::string& text, GpuPrecision& out)
{
	const std::string v = ToLower(text);
	if (v == "auto") { out = GpuPrecision::Auto; return true; }
	if (v == "single") { out = GpuPrecision::Single; return true; }
	if (v == "double") { out = GpuPrecision::Double; return true; }
	return false;
}

void PrintSimulationConfig(const SimulationConfig& c)
{
	using std::cout;
	using std::endl;
	const Eigen::IOFormat rowFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
	const Eigen::IOFormat matFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "; ", "", "", "", "");

	cout << "\tconfig_file              -> " << c.configFile << endl;
	cout << "\t[scene] stl_file         -> " << c.stlFile << " (material " << c.sceneMaterialId << ": " << c.sceneMaterial << ")" << endl;
	for (const GeometryConfig& g : c.geometries)
		cout << "\t[geometry] file          -> " << g.file << " (material " << g.materialId << ": " << g.material << ")" << endl;
	for (const auto& kv : c.surfaceMaterials)
		cout << "\t[surface_materials] surface " << kv.first << " -> material " << kv.second << ": " << c.materials[kv.second].name << endl;
	cout << "\t[scene] save_edges       -> " << (c.saveEdges ? "true" : "false") << endl;
	cout << "\t[scene] save_faces       -> " << (c.saveFaces ? "true" : "false") << endl;
	cout << "\t[transmitter] frequency_hz  -> " << c.frequencyHz << endl;
	for (const TransmitterConfig& tx : c.transmitters)
	{
		const std::string section = c.transmitters.size() > 1 ? "[transmitter " + tx.name + "]" : "[transmitter]";
		cout << "\t" << section << " position   -> " << tx.position.transpose().format(rowFmt) << endl;
		cout << "\t" << section << " coordinate_system -> " << tx.coordSys.format(matFmt) << endl;
		cout << "\t" << section << " power_watts   -> " << tx.powerWatts << endl;
		cout << "\t" << section << " tessellation  -> " << tx.tessellation << endl;
		cout << "\t" << section << " antenna_type  -> " << AntennaTypeName(tx.antennaType) << endl;
		if (!tx.patternFile.empty())
			cout << "\t" << section << " pattern_file  -> " << tx.patternFile << endl;
	}
	cout << "\t[receivers] layout       -> " << ReceiverLayoutName(c.receivers.layout) << endl;
	switch (c.receivers.layout)
	{
	case ReceiverLayout::Line:
		cout << "\t[receivers] line_start   -> " << c.receivers.lineStart.transpose().format(rowFmt) << endl;
		cout << "\t[receivers] line_end     -> " << c.receivers.lineEnd.transpose().format(rowFmt) << endl;
		cout << "\t[receivers] resolution   -> " << c.receivers.resolution << endl;
		break;
	case ReceiverLayout::Plane:
		cout << "\t[receivers] x_min y_min  -> " << c.receivers.xMin << " " << c.receivers.yMin << endl;
		cout << "\t[receivers] x_max y_max  -> " << c.receivers.xMax << " " << c.receivers.yMax << endl;
		cout << "\t[receivers] height       -> " << c.receivers.height << endl;
		cout << "\t[receivers] resolution   -> " << c.receivers.resolution << endl;
		break;
	case ReceiverLayout::File:
		cout << "\t[receivers] file         -> " << c.receivers.file << endl;
		break;
	}
	cout << "\t[ray_tracing] max_reflections   -> " << c.maxReflections << endl;
	cout << "\t[ray_tracing] max_transmissions -> " << c.maxTransmissions << endl;
	cout << "\t[ray_tracing] use_fresnel_coefficients -> " << (c.useFresnelCoefficients ? "true" : "false") << endl;
	cout << "\t[ray_tracing] transmitter_sum   -> " << (c.powerSum ? "power" : "coherent") << endl;
	cout << "\t[ray_tracing] threads           -> " << c.threads << (c.threads == 0 ? " (all hardware threads)" : "") << endl;
	cout << "\t[ray_tracing] device            -> " << TracerDeviceName(c.device) << endl;
	cout << "\t[ray_tracing] gpu_precision     -> " << GpuPrecisionName(c.gpuPrecision) << endl;
	for (size_t i = 0; i < c.materials.size(); i++)
	{
		const MaterialProperties& p = c.materials[i].properties;
		cout << "\t[material " << i << "] " << c.materials[i].name
			<< ": frequency_hz=" << p.frequency
			<< " reflection_loss=" << p.reflectionLoss
			<< " transmission_loss=" << p.transmissionLoss
			<< " rel_conductivity=" << p.relConductivity
			<< " rel_permittivity=" << p.relPermittivityRe << (p.relPermittivityIm >= 0 ? "+" : "") << p.relPermittivityIm << "i" << endl;
	}
	cout << "\t[outputs] ray_paths=" << (c.writeRayPaths ? "true" : "false")
		<< " receivers=" << (c.writeReceivers ? "true" : "false")
		<< " field_csv=" << (c.writeFieldCsv ? "true" : "false")
		<< " icosahedron=" << (c.writeIcosahedron ? "true" : "false")
		<< " paths_csv=" << (c.writePathsCsv ? "true" : "false")
		<< " viewer=" << (c.viewer ? "true" : "false") << endl;
}
