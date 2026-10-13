#include "SimulationConfig.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

namespace fs = std::filesystem;

namespace
{
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

	bool ParseDoubles(const std::string& text, std::vector<double>& out)
	{
		std::string cleaned = text;
		std::replace(cleaned.begin(), cleaned.end(), ',', ' ');
		std::replace(cleaned.begin(), cleaned.end(), ';', ' ');
		std::istringstream iss(cleaned);
		double v;
		while (iss >> v) out.push_back(v);
		return iss.eof() && !out.empty();
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
				_diag.Add(Line(key), "key '" + key + "' expects one number, got '" + text + "'");
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
			long long value;
			if (!(iss >> value) || !(iss >> std::ws).eof())
			{
				_diag.Add(Line(key), "key '" + key + "' expects an integer, got '" + text + "'");
				return false;
			}
			out = static_cast<int>(value);
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
				_diag.Add(Line(key), "key '" + key + "' expects three numbers 'x y z', got '" + text + "'");
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
				_diag.Add(Line(key), "key '" + key + "' expects nine numbers (three rows of three), got '" + text + "'");
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
			if (mustExist && !fs::exists(p))
			{
				_diag.Add(Line(key), "key '" + key + "': file not found: " + p.generic_string());
				return false;
			}
			out = p.generic_string();
			return true;
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
		if (v == "array") { out = AntennaType::Array; return true; }
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

	void ReadScene(SectionReader& r, SimulationConfig& cfg)
	{
		r.Path("stl_file", cfg.stlFile, true, true);
		r.Boolean("save_edges", cfg.saveEdges, false);
		r.Boolean("save_faces", cfg.saveFaces, false);
	}

	void ReadTransmitter(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		r.Vector3("position", cfg.txPosition, true);
		r.Matrix3("coordinate_system", cfg.txCoordSys, false);
		if (r.Double("power_watts", cfg.txPowerWatts, false) && cfg.txPowerWatts <= 0.0)
			diag.Add(r.Line("power_watts"), "power_watts must be positive");
		if (r.Double("frequency_hz", cfg.frequencyHz, true) && cfg.frequencyHz <= 0.0)
			diag.Add(r.Line("frequency_hz"), "frequency_hz must be positive");
		if (r.Integer("tessellation", cfg.tessellation, false) && cfg.tessellation < 0)
			diag.Add(r.Line("tessellation"), "tessellation must be 0 or greater");
		std::string antenna;
		if (r.String("antenna_type", antenna, false) && !ParseAntennaType(antenna, cfg.antennaType))
			diag.Add(r.Line("antenna_type"), "antenna_type must be one of isotropic, hertzian_dipole, half_wave_dipole, array; got '" + antenna + "'");
		r.Path("pattern_file", cfg.patternFile, false, true);
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
		switch (rx.layout)
		{
		case ReceiverLayout::Line:
			r.Vector3("line_start", rx.lineStart, true);
			r.Vector3("line_end", rx.lineEnd, true);
			if (r.Double("resolution", rx.resolution, true) && rx.resolution <= 0.0)
				diag.Add(r.Line("resolution"), "resolution must be positive");
			break;
		case ReceiverLayout::Plane:
			r.Double("x_min", rx.xMin, true);
			r.Double("y_min", rx.yMin, true);
			r.Double("x_max", rx.xMax, true);
			r.Double("y_max", rx.yMax, true);
			r.Double("height", rx.height, true);
			if (r.Double("resolution", rx.resolution, true) && rx.resolution <= 0.0)
				diag.Add(r.Line("resolution"), "resolution must be positive");
			if (rx.xMax < rx.xMin || rx.yMax < rx.yMin)
				diag.Add(r.Line("x_max"), "plane bounds must satisfy x_min <= x_max and y_min <= y_max");
			break;
		case ReceiverLayout::File:
			r.Path("file", rx.file, true, true);
			break;
		}
	}

	void ReadRayTracing(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		if (r.Integer("max_reflections", cfg.maxReflections, false) && cfg.maxReflections < 0)
			diag.Add(r.Line("max_reflections"), "max_reflections must be 0 or greater");
		if (r.Integer("max_transmissions", cfg.maxTransmissions, false) && cfg.maxTransmissions < 0)
			diag.Add(r.Line("max_transmissions"), "max_transmissions must be 0 or greater");
		r.Boolean("use_fresnel_coefficients", cfg.useFresnelCoefficients, false);
	}

	void ReadMaterial(SectionReader& r, SimulationConfig& cfg, Diagnostics& diag)
	{
		MaterialConfig m;
		MaterialProperties& p = m.properties;
		r.String("name", m.name, true);
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

	// Sections are read in a fixed order so that [material] can default its
	// frequency to the transmitter frequency regardless of file order.
	const char* order[] = { "scene", "transmitter", "receivers", "ray_tracing", "material", "outputs" };
	bool seen[6] = { false, false, false, false, false, false };
	for (int i = 0; i < 6; i++)
	{
		const std::string name = order[i];
		for (Section& section : sections)
		{
			if (section.name != name) continue;
			if (seen[i] && name != "material")
			{
				diag.Add(section.line, "section [" + name + "] appears more than once");
				continue;
			}
			seen[i] = true;
			SectionReader reader(section, diag, baseDir);
			if (name == "scene") ReadScene(reader, config);
			else if (name == "transmitter") ReadTransmitter(reader, config, diag);
			else if (name == "receivers") ReadReceivers(reader, config, diag);
			else if (name == "ray_tracing") ReadRayTracing(reader, config, diag);
			else if (name == "material") ReadMaterial(reader, config, diag);
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
		diag.Add(0, "at least one [material] section is required (material 0 is the default for every surface)");

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
	case AntennaType::Array: return "array";
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

void PrintSimulationConfig(const SimulationConfig& c)
{
	using std::cout;
	using std::endl;
	const Eigen::IOFormat rowFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
	const Eigen::IOFormat matFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "; ", "", "", "", "");

	cout << "\tconfig_file              -> " << c.configFile << endl;
	cout << "\t[scene] stl_file         -> " << c.stlFile << endl;
	cout << "\t[scene] save_edges       -> " << (c.saveEdges ? "true" : "false") << endl;
	cout << "\t[scene] save_faces       -> " << (c.saveFaces ? "true" : "false") << endl;
	cout << "\t[transmitter] position   -> " << c.txPosition.transpose().format(rowFmt) << endl;
	cout << "\t[transmitter] coordinate_system -> " << c.txCoordSys.format(matFmt) << endl;
	cout << "\t[transmitter] power_watts   -> " << c.txPowerWatts << endl;
	cout << "\t[transmitter] frequency_hz  -> " << c.frequencyHz << endl;
	cout << "\t[transmitter] tessellation  -> " << c.tessellation << endl;
	cout << "\t[transmitter] antenna_type  -> " << AntennaTypeName(c.antennaType) << endl;
	if (!c.patternFile.empty())
		cout << "\t[transmitter] pattern_file  -> " << c.patternFile << endl;
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
		<< " icosahedron=" << (c.writeIcosahedron ? "true" : "false") << endl;
}
