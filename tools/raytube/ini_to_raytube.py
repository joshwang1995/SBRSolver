#!/usr/bin/env python3
"""Convert an SBRSolver INI configuration into a RayTubeTracing run directory.

RayTubeTracing (Takahiro Hashimoto's ray-tube tracer, ``RayTubeTracing.exe``)
is the reference solver that SBRSolver's results are checked against. It
reads a ``Config.xml`` plus a single ASCII STL, a receiver list and antenna
pattern tables from its working directory. This script writes all of them
from an SBRSolver INI file so that both solvers simulate the same scene:

* ``Config.xml``          the reference configuration (same keys as the
                          supplied sample, with the values taken from the INI)
* ``geometry.stl``        the scene STL files merged in SBRSolver's order
* ``receivers.dat``       the receiver positions SBRSolver generates for the
                          INI's ``[receivers]`` layout, in the same order
* ``tx_pattern.dat``      the analytic transmitter pattern of ``antenna_type``
                          tabulated for the reference (``sqrt(gain)``, 0 phase)
* ``rx_pattern.dat``      an isotropic receiver pattern (fields are compared,
                          not received powers)

Mapping notes

* SBRSolver's ``[transmitter] coordinate_system`` rows are the TX frame axes in
  global coordinates; row 2 becomes ``tx_axis.z`` and row 0 ``tx_axis.x``.
* ``power_watts`` becomes ``transmitting_power_dBm``.
* Every ``[material]`` becomes a RayTubeTracing material with Fresnel
  reflection/transmission (type 1) and ``epsilon_setting_mode`` 1 (static
  permittivity ``rel_permittivity_re`` and conductivity ``rel_conductivity``),
  which is exactly how SBRSolver builds its complex permittivity.
* Per-file materials and ``[surface_materials]`` overrides are turned into
  ``face_index`` lists. The coplanar surface ids are recomputed here with the
  same rule SBRSolver uses (triangles sharing an edge, with equal or opposite
  facet normals, in file order); ``--check-faces Faces.vtk`` verifies the
  result against the ``MaterialID`` field SBRSolver wrote.
* ``use_fresnel_coefficients = false`` (every surface a perfect conductor in
  SBRSolver) becomes a material with conductivity 1e30, the Fresnel limit of
  a perfect conductor in the reference.
* ``tessellation`` n (10*4^n+2 launch rays) becomes ``num_tesselate`` 2^n,
  which launches the same number of ray-tube corner rays. Both solvers
  reconstruct the exact image-method path for every captured ray, so the
  result does not depend on this value once every path is found.

Usage: ini_to_raytube.py <config.ini> <output_dir> [--check-faces Faces.vtk]
                         [--tesselate N] [--threads N] [--pattern-step DEG]
"""
import argparse
import math
import os
import re
import struct
import sys
from collections import deque

ETA0 = 4.0 * math.pi * 1e-7 * 299792458.0  # free-space impedance used by both solvers

# Directivities of the analytic antennas (SBRSolver's FieldCompute uses the
# same constants). The half-wave value is 4*pi / integral of |cos(pi/2 cos t)/sin t|^2.
DIRECTIVITY = {"isotropic": 1.0, "hertzian_dipole": 1.5, "half_wave_dipole": 1.6409223769845667}


def fail(message):
    sys.stderr.write("ini_to_raytube: " + message + "\n")
    sys.exit(1)


# ---------------------------------------------------------------------------
# INI reading (same rules as SBRSolver_rev2/SimulationConfig.cpp)
# ---------------------------------------------------------------------------
class Section:
    def __init__(self, name, line):
        self.name = name
        self.line = line
        self.entries = {}

    def get(self, key, default=None):
        return self.entries.get(key, default)

    def require(self, key, ini):
        if key not in self.entries:
            fail("%s:%d: [%s] needs key '%s'" % (ini, self.line, self.name, key))
        return self.entries[key]


REPEATABLE = ("material", "geometry")


def read_ini(path):
    """Same rules as SimulationConfig.cpp: comments, a leading BOM, repeatable
    [material]/[geometry] sections only, no duplicate keys."""
    sections = []
    current = None
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for number, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line[0] in "#;":
                continue
            for i in range(1, len(line)):
                if line[i] == "#" and line[i - 1].isspace():
                    line = line[:i]
                    break
            line = line.strip()
            if not line:
                continue
            if line[0] == "[":
                if line[-1] != "]":
                    fail("%s:%d: malformed section header '%s'" % (path, number, line))
                name = line[1:-1].strip().lower()
                if not name:
                    fail("%s:%d: empty section name" % (path, number))
                if name not in REPEATABLE and any(sec.name == name for sec in sections):
                    fail("%s:%d: section [%s] appears more than once" % (path, number, name))
                current = Section(name, number)
                sections.append(current)
                continue
            if "=" not in line:
                fail("%s:%d: expected 'key = value'" % (path, number))
            if current is None:
                fail("%s:%d: key outside of any [section]" % (path, number))
            key, value = line.split("=", 1)
            key = key.strip().lower()
            if not key:
                fail("%s:%d: empty key" % (path, number))
            if key in current.entries:
                fail("%s:%d: duplicate key '%s' in [%s]" % (path, number, key, current.name))
            current.entries[key] = value.strip()
    return sections


def numbers(text, count, what):
    values = []
    for token in text.replace(";", " ").split():
        try:
            values.append(float(token))
        except ValueError:
            fail("%s: '%s' is not a number" % (what, token))
    if count is not None and len(values) != count:
        fail("%s: expected %d numbers, got %d" % (what, count, len(values)))
    return values


def boolean(text, what):
    t = text.strip().lower()
    if t in ("true", "yes", "on", "1"):
        return True
    if t in ("false", "no", "off", "0"):
        return False
    fail("%s: expected true or false, got '%s'" % (what, text))


# ---------------------------------------------------------------------------
# STL reading (facets in file order, coordinates as float32 like stl_reader)
# ---------------------------------------------------------------------------
def f32(x):
    return struct.unpack("<f", struct.pack("<f", float(x)))[0]


def read_stl(path):
    """Returns a list of (normal, (v0, v1, v2)) with float32-rounded values."""
    with open(path, "rb") as fh:
        data = fh.read()
    facets = []
    # ASCII when the file starts with "solid" and contains "facet" (stl_reader sniffs similarly)
    head = data[:512].lstrip()
    is_ascii = head.startswith(b"solid") and b"facet" in data[:4096] or b"vertex" in data[:4096]
    if is_ascii:
        text = data.decode("utf-8", errors="replace")
        normal = None
        verts = []
        for raw in text.splitlines():
            tokens = raw.split()
            if not tokens:
                continue
            if tokens[0] == "facet" and len(tokens) >= 5 and tokens[1] == "normal":
                normal = tuple(f32(t) for t in tokens[2:5])
                verts = []
            elif tokens[0] == "vertex" and len(tokens) >= 4:
                verts.append(tuple(f32(t) for t in tokens[1:4]))
            elif tokens[0] == "endfacet":
                if normal is None or len(verts) != 3:
                    fail("%s: malformed facet" % path)
                facets.append((normal, tuple(verts)))
                normal = None
                verts = []
    else:
        if len(data) < 84:
            fail("%s: binary STL is too short" % path)
        count = struct.unpack("<I", data[80:84])[0]
        if len(data) < 84 + 50 * count:
            fail("%s: binary STL is truncated" % path)
        for i in range(count):
            off = 84 + 50 * i
            values = struct.unpack("<12f", data[off:off + 48])
            facets.append((tuple(values[0:3]), (tuple(values[3:6]), tuple(values[6:9]), tuple(values[9:12]))))
    # stl_reader drops a facet whose corners are not three distinct points
    # (after float32 rounding) and keeps the numbering of the others.
    facets = [f for f in facets if len(set(f[1])) == 3]
    if not facets:
        fail("%s: no triangles" % path)
    return facets


def coplanar_ids(files):
    """Replicates Preprocessor::StlToGeometry + BfsCoplanarSurface.

    files: list of facet lists (one per STL file, in SBRSolver's order).
    Returns the coplanar id of every triangle in merged order.
    """
    tris = []        # (file index, unified vertex ids, normal)
    for file_index, facets in enumerate(files):
        ids = {}
        unique = []
        for normal, verts in facets:
            idx = []
            for v in verts:
                if v not in ids:
                    ids[v] = len(unique)
                    unique.append(v)
                idx.append(ids[v])
            tris.append((file_index, tuple(idx), normal))
    edge_map = {}
    adjacency = [[] for _ in tris]
    for t, (file_index, idx, normal) in enumerate(tris):
        for a, b in ((idx[0], idx[1]), (idx[1], idx[2]), (idx[0], idx[2])):
            key = (file_index, min(a, b), max(a, b))
            if key in edge_map:
                for other in edge_map[key]:
                    on = tris[other][2]
                    if on == normal or on == tuple(-c for c in normal):
                        adjacency[other].append(t)
                        adjacency[t].append(other)
                edge_map[key].append(t)
            else:
                edge_map[key] = [t]
    result = [-1] * len(tris)
    visited = [False] * len(tris)
    next_id = 0
    for start in range(len(tris)):
        if visited[start]:
            continue
        visited[start] = True
        queue = deque([start])
        while queue:
            u = queue.popleft()
            result[u] = next_id
            for v in adjacency[u]:
                if not visited[v]:
                    visited[v] = True
                    queue.append(v)
        next_id += 1
    return result


def read_faces_vtk(path):
    """Returns (coplanar ids, material ids) from SBRSolver's Faces.vtk."""
    fields = {}
    with open(path) as fh:
        lines = fh.read().split("\n")
    i = 0
    while i < len(lines):
        tokens = lines[i].split()
        if len(tokens) == 4 and tokens[0] in ("FaceID", "CoplanarID", "MaterialID") and tokens[3] == "int":
            n = int(tokens[2])
            values = []
            i += 1
            while len(values) < n and i < len(lines):
                values.extend(int(v) for v in lines[i].split())
                i += 1
            fields[tokens[0]] = values
            continue
        i += 1
    if "CoplanarID" not in fields or "MaterialID" not in fields:
        fail("%s: no CoplanarID/MaterialID fields (run sbrsolver with save_faces = true)" % path)
    return fields["CoplanarID"], fields["MaterialID"]


# ---------------------------------------------------------------------------
# Receivers (same arithmetic as Preprocessor::GenerateRXLine / GenerateRxPlane)
# ---------------------------------------------------------------------------
def receivers_line(start, end, resolution):
    direction = [e - s for s, e in zip(start, end)]
    length = math.sqrt(sum(d * d for d in direction))
    count = int(math.floor(length / resolution))
    direction = [d / length for d in direction]
    return [tuple(s + (i * resolution) * d for s, d in zip(start, direction)) for i in range(count)]


def receivers_plane(x_min, y_min, x_max, y_max, height, resolution):
    nx = int((x_max - x_min) / resolution)
    ny = int((y_max - y_min) / resolution)
    points = []
    for i in range(nx + 1):
        x = x_min + i * resolution
        for j in range(ny + 1):
            points.append((x, y_min + j * resolution, height))
    return points


def receivers_file(path):
    points = []
    with open(path) as fh:
        for raw in fh:
            tokens = raw.split()
            if len(tokens) != 3:
                continue
            try:
                points.append(tuple(float(t) for t in tokens))
            except ValueError:
                continue
    return points


# ---------------------------------------------------------------------------
# Antenna pattern tables
# ---------------------------------------------------------------------------
def pattern_amplitude(kind, theta):
    if kind == "isotropic":
        return 1.0
    if kind == "hertzian_dipole":
        return math.sqrt(DIRECTIVITY[kind]) * math.sin(theta)
    if kind == "half_wave_dipole":
        # cos(pi/2 cos t) / sin t evaluated as in FieldCompute::PatternAmplitude
        s = math.sin(theta)
        if s == 0.0:
            return 0.0
        h, c = math.sin(0.5 * theta), math.cos(0.5 * theta)
        return math.sqrt(DIRECTIVITY[kind]) * math.sin(0.5 * math.pi * 2.0 * min(h * h, c * c)) / s
    fail("unknown antenna type '%s'" % kind)


def write_pattern(path, kind, step_deg):
    """Writes a RayTubeTracing pattern table: header, 'nTheta nPhi', then
    'theta phi |Etheta|dB arg(Etheta) |Ephi|dB arg(Ephi)' rows (phi-major)."""
    n_theta = int(round(180.0 / step_deg)) + 1
    floor_db = -6000.0  # 10^(-300): effectively zero but still a finite table entry
    with open(path, "w", newline="\n") as fh:
        fh.write("Theta[deg] Phi[deg] Abs(Theta)[dB] Phase(Theta)[deg] Abs(Phi)[dB] Phase(Phi)[deg]\n")
        fh.write("%d 1\n" % n_theta)
        for i in range(n_theta):
            theta_deg = 180.0 * i / (n_theta - 1)
            amp = pattern_amplitude(kind, math.radians(theta_deg))
            db = 20.0 * math.log10(amp) if amp > 1e-300 else floor_db
            fh.write("%.17g 0 %.17g 0 %g 0\n" % (theta_deg, db, floor_db))


# ---------------------------------------------------------------------------
# Config.xml
# ---------------------------------------------------------------------------
PEC_CONDUCTIVITY = 1e30  # Fresnel limit of a perfect conductor: rPerp -> -1, rPar -> +1, no transmission


def material_xml(index, m, faces, perfect_conductor):
    lines = [
        "\t<material%d>" % index,
        "\t\t<!-- %s%s -->" % (m["name"], " (use_fresnel_coefficients = false: perfect conductor)" if perfect_conductor else ""),
        "\t\t<reflection_type>1</reflection_type>",
        "\t\t<transmission_type>1</transmission_type>",
        "\t\t<epsilon_setting_mode>1</epsilon_setting_mode>",
        "\t\t<complex_permittivity> <real></real><imag></imag> </complex_permittivity>",
        "\t\t<static_permittivity>%.17g</static_permittivity>" % (1.0 if perfect_conductor else m["rel_permittivity_re"]),
        "\t\t<conductivity>%.17g</conductivity>" % (PEC_CONDUCTIVITY if perfect_conductor else m["rel_conductivity"]),
        "\t\t<is_impedance_boundary>false</is_impedance_boundary>",
        "\t\t<thickness>0</thickness>",
        "\t\t<is_diffuse_reflect>false</is_diffuse_reflect>",
        "\t\t<is_diffuse_transmit>false</is_diffuse_transmit>",
    ]
    if index > 0:
        lines.append("\t\t<face_index>%s</face_index>" % ",".join(str(f) for f in faces))
    lines.append("\t</material%d>" % index)
    return "\n".join(lines)


def vec_xml(tag, v):
    return "\t<%s> <x>%.17g</x> <y>%.17g</y> <z>%.17g</z> </%s>" % (tag, v[0], v[1], v[2], tag)


def write_config(path, cfg):
    tx_z = cfg["tx_z"]
    tx_x = cfg["tx_x"]
    xml = []
    xml.append("<!-- RayTubeTracing configuration generated by tools/raytube/ini_to_raytube.py from %s -->" % cfg["ini_name"])
    xml.append("<frequency>%.17g</frequency>" % cfg["frequency"])
    xml.append("")
    xml.append("<!-- geometry -->")
    xml.append("\t<is_free_space>false</is_free_space>")
    xml.append("\t<geometry_filename>geometry.stl</geometry_filename>")
    xml.append("\t<num_material>%d</num_material>" % len(cfg["materials"]))
    for i, m in enumerate(cfg["materials"]):
        xml.append(material_xml(i, m, cfg["material_faces"][i], cfg["perfect_conductor"]))
    xml.append("")
    xml.append("\t<reflection_coef_filename></reflection_coef_filename>")
    xml.append("\t<transmission_coef_filename></transmission_coef_filename>")
    xml.append("\t<diffuse_reflection_filename></diffuse_reflection_filename>")
    xml.append("\t<diffuse_transmission_filename></diffuse_transmission_filename>")
    xml.append("\t<is_diffract_isolated_edges>false</is_diffract_isolated_edges>")
    xml.append("")
    xml.append("<!-- source -->")
    xml.append(vec_xml("source_position", cfg["tx_position"]))
    xml.append("\t<transmitting_power_dBm>%.17g</transmitting_power_dBm>" % cfg["power_dbm"])
    xml.append("\t<tx_antenna_pattern_file>tx_pattern.dat</tx_antenna_pattern_file>")
    xml.append("\t<tx_axis>")
    xml.append("\t\t<zx>%.17g</zx> <zy>%.17g</zy> <zz>%.17g</zz>" % tuple(tx_z))
    xml.append("\t\t<xx>%.17g</xx> <xy>%.17g</xy> <xz>%.17g</xz>" % tuple(tx_x))
    xml.append("\t</tx_axis>")
    xml.append("")
    xml.append("<!-- launching -->")
    xml.append("\t<ray_launch_mode>1</ray_launch_mode>")
    xml.append("\t<num_tesselate>%d</num_tesselate>" % cfg["tesselate"])
    xml.append("\t<num_tesselate_scat>1</num_tesselate_scat>")
    xml.append("\t<num_theta>2</num_theta>")
    xml.append("\t<num_phi>1</num_phi>")
    xml.append("\t<theta_phi_launch_zaxis> <x>0.0</x> <y>0.0</y> <z>1.0</z> </theta_phi_launch_zaxis>")
    xml.append("\t<theta_phi_launch_xaxis> <x>1.0</x> <y>0.0</y> <z>0.0</z> </theta_phi_launch_xaxis>")
    xml.append("\t<is_filter_rays_by_theta_angle>false</is_filter_rays_by_theta_angle>")
    xml.append("\t<is_filter_rays_by_phi_angle>false</is_filter_rays_by_phi_angle>")
    xml.append("\t<launch_angle_theta_start_deg></launch_angle_theta_start_deg>")
    xml.append("\t<launch_angle_theta_end_deg></launch_angle_theta_end_deg>")
    xml.append("\t<launch_angle_phi_start_deg></launch_angle_phi_start_deg>")
    xml.append("\t<launch_angle_phi_end_deg></launch_angle_phi_end_deg>")
    xml.append("\t<launch_filter_zaxis> <x>0.0</x> <y>0.0</y> <z>1.0</z> </launch_filter_zaxis>")
    xml.append("\t<launch_filter_xaxis> <x>1.0</x> <y>0.0</y> <z>0.0</z> </launch_filter_xaxis>")
    xml.append("\t<is_filter_rays_by_launch_id>false</is_filter_rays_by_launch_id>")
    xml.append("\t<ray_id_to_launch></ray_id_to_launch>")
    xml.append("")
    xml.append("<!-- tracing -->")
    xml.append("\t<num_reflection_max>%d</num_reflection_max>" % cfg["max_reflections"])
    xml.append("\t<num_transmission_max>%d</num_transmission_max>" % cfg["max_transmissions"])
    xml.append("\t<num_diffraction_max>0</num_diffraction_max>")
    xml.append("\t<num_diffuse_reflection_max>0</num_diffuse_reflection_max>")
    xml.append("\t<num_diffuse_transmission_max>0</num_diffuse_transmission_max>")
    xml.append("\t<num_ray_diffraction>1</num_ray_diffraction>")
    xml.append("\t<is_restrict_diffracting_edge>false</is_restrict_diffracting_edge>")
    xml.append("\t<diffracting_edge_index></diffracting_edge_index>")
    xml.append("\t<tracing_field_threshold>0</tracing_field_threshold>")
    xml.append("\t<diffuse_reflection_field_threshold>0</diffuse_reflection_field_threshold>")
    xml.append("\t<diffuse_transmission_field_threshold>0</diffuse_transmission_field_threshold>")
    xml.append("")
    xml.append("<!-- capturing -->")
    xml.append("\t<receiver_position_load_type>0</receiver_position_load_type>")
    xml.append("\t<receiver_position_filename>receivers.dat</receiver_position_filename>")
    xml.append("\t<is_limit_receivers>false</is_limit_receivers>")
    xml.append("\t<receiver_id_limit>0</receiver_id_limit>")
    xml.append("\t<rx_antenna_pattern_file>rx_pattern.dat</rx_antenna_pattern_file>")
    xml.append("\t<rx_axis>")
    xml.append("\t\t<zx>0.0</zx> <zy>0.0</zy> <zz>1.0</zz>")
    xml.append("\t\t<xx>1.0</xx> <xy>0.0</xy> <xz>0.0</xz>")
    xml.append("\t</rx_axis>")
    xml.append("\t<is_filter_nan_rays>false</is_filter_nan_rays>")
    xml.append("\t<is_filter_direct_rays>false</is_filter_direct_rays>")
    xml.append("\t<is_filter_low_power_rays>false</is_filter_low_power_rays>")
    xml.append("\t<power_truncation_threshold>0</power_truncation_threshold>")
    xml.append("\t<is_remove_duplicate_rays>true</is_remove_duplicate_rays>")
    xml.append("")
    xml.append("<!-- output -->")
    xml.append("\t<is_save_receivers_csv>true</is_save_receivers_csv>")
    xml.append("\t<is_save_electric_field_csv>true</is_save_electric_field_csv>")
    xml.append("\t<is_save_electric_field_matlab_csv>true</is_save_electric_field_matlab_csv>")
    xml.append("\t<is_save_facemodel_vtk>true</is_save_facemodel_vtk>")
    xml.append("\t<is_save_edgemodel_vtk>false</is_save_edgemodel_vtk>")
    xml.append("\t<is_save_transmitter_vtk>false</is_save_transmitter_vtk>")
    xml.append("\t<is_save_tessellated_icosahedron_vtk>false</is_save_tessellated_icosahedron_vtk>")
    xml.append("\t<is_save_raytube_vtk>false</is_save_raytube_vtk>")
    xml.append("\t<is_save_receivers_vtk>false</is_save_receivers_vtk>")
    xml.append("\t<is_save_raypath_vtk>%s</is_save_raypath_vtk>" % ("true" if cfg["save_paths"] else "false"))
    xml.append("\t<vtk_length_max>10.0</vtk_length_max>")
    xml.append("")
    xml.append("<!-- system -->")
    xml.append("\t<num_openmp_cpus>%d</num_openmp_cpus>" % cfg["threads"])
    xml.append("\t<num_log_interval>1000000</num_log_interval>")
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(xml) + "\n")


def write_stl(path, files):
    with open(path, "w", newline="\n") as fh:
        fh.write("solid merged\n")
        for facets in files:
            for normal, verts in facets:
                # 17 digits: std::stod in the reference then yields exactly the
                # float32 values SBRSolver works with (%.9g would not).
                fh.write("facet normal %.17g %.17g %.17g\n" % normal)
                fh.write("outer loop\n")
                for v in verts:
                    fh.write("vertex %.17g %.17g %.17g\n" % v)
                fh.write("endloop\n")
                fh.write("endfacet\n")
        fh.write("endsolid merged\n")


def write_receivers(path, points):
    with open(path, "w", newline="\n") as fh:
        fh.write("receiver positions\n")
        fh.write("%d\n" % len(points))
        for p in points:
            fh.write("%.17g %.17g %.17g\n" % p)


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ini")
    parser.add_argument("output_dir")
    parser.add_argument("--check-faces", metavar="FACES_VTK",
                        help="verify the per-face materials against the MaterialID field of SBRSolver's Faces.vtk")
    parser.add_argument("--tesselate", type=int, default=0,
                        help="num_tesselate for the reference (default: 2^tessellation of the INI)")
    parser.add_argument("--threads", type=int, default=1, help="num_openmp_cpus for the reference (default 1)")
    parser.add_argument("--pattern-step", type=float, default=0.0025,
                        help="theta spacing in degrees of the tabulated analytic pattern (default 0.0025)")
    parser.add_argument("--save-paths", action="store_true", help="let the reference write CapturedRays.vtk")
    args = parser.parse_args()

    ini = os.path.abspath(args.ini)
    base = os.path.dirname(ini)
    sections = read_ini(ini)

    def first(name):
        for s in sections:
            if s.name == name:
                return s
        return None

    def all_named(name):
        return [s for s in sections if s.name == name]

    scene = first("scene") or fail("%s: no [scene] section" % ini)
    tx = first("transmitter") or fail("%s: no [transmitter] section" % ini)
    rx = first("receivers") or fail("%s: no [receivers] section" % ini)
    tracing = first("ray_tracing")
    materials = all_named("material")
    if not materials:
        fail("%s: no [material] section" % ini)

    cfg = {"ini_name": os.path.basename(ini)}

    # materials
    names = {}
    mats = []
    for s in materials:
        name = s.require("name", ini)
        if name in names:
            fail("%s:%d: duplicate material name '%s'" % (ini, s.line, name))
        names[name] = len(mats)
        mats.append({
            "name": name,
            "rel_permittivity_re": numbers(s.require("rel_permittivity_re", ini), 1, "rel_permittivity_re")[0],
            "rel_conductivity": numbers(s.require("rel_conductivity", ini), 1, "rel_conductivity")[0],
        })
    cfg["materials"] = mats

    def material_index(text, where):
        if text not in names:
            fail("%s: unknown material '%s'" % (where, text))
        return names[text]

    # geometry files in SBRSolver's order with their per-file material
    sources = [(os.path.join(base, scene.require("stl_file", ini)),
                material_index(scene.get("material", mats[0]["name"]), "[scene] material"))]
    for g in all_named("geometry"):
        sources.append((os.path.join(base, g.require("file", ini)), material_index(g.require("material", ini), "[geometry] material")))
    files = [read_stl(path) for path, _ in sources]
    face_material = []
    for (path, material), facets in zip(sources, files):
        face_material.extend([material] * len(facets))
    ids = coplanar_ids(files)
    overrides = first("surface_materials")
    if overrides is not None:
        for key, value in overrides.entries.items():
            if not re.fullmatch(r"\d+", key):
                fail("%s:%d: [surface_materials] key '%s' is not a surface id" % (ini, overrides.line, key))
            surface = int(key)
            if surface > max(ids):
                fail("%s: [surface_materials] surface %d does not exist (ids 0 to %d)" % (ini, surface, max(ids)))
            material = material_index(value, "[surface_materials] %s" % key)
            for f in range(len(ids)):
                if ids[f] == surface:
                    face_material[f] = material
    if args.check_faces:
        vtk_ids, vtk_materials = read_faces_vtk(args.check_faces)
        if vtk_ids != ids or vtk_materials != face_material:
            bad = [f for f in range(min(len(ids), len(vtk_ids))) if vtk_ids[f] != ids[f] or vtk_materials[f] != face_material[f]]
            fail("per-face materials differ from %s (first differing faces: %s); the converter's coplanar grouping "
                 "no longer matches SBRSolver's" % (args.check_faces, bad[:10]))
    # RayTubeTracing gives material0 to every face and lets the other
    # materials override their face_index lists, each face at most once. A
    # material without faces would make it abort (empty face_index), so only
    # the materials that are actually used are written: the one with most
    # faces becomes material0, the others list their faces.
    used = sorted(set(face_material), key=lambda m: (-face_material.count(m), m))
    cfg["materials"] = [mats[m] for m in used]
    cfg["material_faces"] = [[f for f in range(len(face_material)) if face_material[f] == m] for m in used]

    # transmitter
    cfg["tx_position"] = numbers(tx.require("position", ini), 3, "position")
    axes = numbers(tx.get("coordinate_system", "1 0 0; 0 1 0; 0 0 1"), 9, "coordinate_system")
    cfg["tx_x"] = axes[0:3]
    cfg["tx_z"] = axes[6:9]
    power = numbers(tx.get("power_watts", "1"), 1, "power_watts")[0]
    cfg["power_dbm"] = 10.0 * math.log10(power * 1000.0)
    cfg["frequency"] = numbers(tx.require("frequency_hz", ini), 1, "frequency_hz")[0]
    tessellation = int(numbers(tx.get("tessellation", "2"), 1, "tessellation")[0])
    cfg["tesselate"] = args.tesselate if args.tesselate > 0 else 2 ** tessellation
    antenna = tx.get("antenna_type", "hertzian_dipole").lower()
    if antenna not in DIRECTIVITY:
        fail("%s: antenna_type '%s' is not supported" % (ini, antenna))
    if tx.get("pattern_file"):
        sys.stderr.write("ini_to_raytube: note: pattern_file is ignored, the analytic %s pattern is used\n" % antenna)

    # receivers
    layout = rx.get("layout", "line").lower()
    if layout == "line":
        points = receivers_line(numbers(rx.require("line_start", ini), 3, "line_start"),
                                numbers(rx.require("line_end", ini), 3, "line_end"),
                                numbers(rx.require("resolution", ini), 1, "resolution")[0])
    elif layout == "plane":
        points = receivers_plane(*(numbers(rx.require(k, ini), 1, k)[0]
                                   for k in ("x_min", "y_min", "x_max", "y_max", "height", "resolution")))
    elif layout == "file":
        points = receivers_file(os.path.join(base, rx.require("file", ini)))
    else:
        fail("%s: unknown receiver layout '%s'" % (ini, layout))
    if not points:
        fail("%s: the [receivers] section produced no receiver positions" % ini)

    cfg["max_reflections"] = int(numbers(tracing.get("max_reflections", "1"), 1, "max_reflections")[0]) if tracing else 1
    cfg["max_transmissions"] = int(numbers(tracing.get("max_transmissions", "0"), 1, "max_transmissions")[0]) if tracing else 0
    # use_fresnel_coefficients = false makes SBRSolver treat every surface as
    # a perfect conductor; the reference gets the same through the Fresnel
    # limit of an (almost) infinitely conductive material.
    cfg["perfect_conductor"] = bool(tracing) and not boolean(tracing.get("use_fresnel_coefficients", "true"), "use_fresnel_coefficients")
    cfg["threads"] = args.threads
    cfg["save_paths"] = args.save_paths

    os.makedirs(args.output_dir, exist_ok=True)
    write_stl(os.path.join(args.output_dir, "geometry.stl"), files)
    write_receivers(os.path.join(args.output_dir, "receivers.dat"), points)
    write_pattern(os.path.join(args.output_dir, "tx_pattern.dat"), antenna, args.pattern_step)
    write_pattern(os.path.join(args.output_dir, "rx_pattern.dat"), "isotropic", 90.0)
    write_config(os.path.join(args.output_dir, "Config.xml"), cfg)
    print("ini_to_raytube: wrote %s (%d triangles, %d receivers, %d materials, num_tesselate %d)"
          % (args.output_dir, len(face_material), len(points), len(mats), cfg["tesselate"]))


if __name__ == "__main__":
    main()
