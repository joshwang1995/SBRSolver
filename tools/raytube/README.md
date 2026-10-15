# Comparison with the RayTubeTracing reference solver

`RayTubeTracing.exe` is Takahiro Hashimoto's ray-tube tracer ((c) 2020, MSVC
build of 2020-11-04). It is the reference that SBRSolver's results are checked
against: the same scenes are simulated with both programs and the electric
field at every receiver must agree to every digit the reference prints.

The executable only runs on Windows (or under Wine). The comparison therefore
works in two modes:

* **with the executable** (`ctest` on Windows, or anywhere with
  `-DSBR_RAYTUBE_EXE=tools/raytube/RayTubeTracing.exe -DSBR_RAYTUBE_LAUNCHER=wine64`):
  the reference is run for every scene and its output compared with
  sbrsolver's, and also with the committed reference output below;
* **without it** (`ctest` on Linux and macOS): sbrsolver's output is compared
  with the committed reference output in `SBRSolver_rev2/data/reference/raytube/<scene>/`.

## Files

| file | purpose |
|---|---|
| `ini_to_raytube.py` | turns an SBRSolver INI into a RayTubeTracing run directory: `Config.xml`, the scene STL files merged into one, the receiver list, and the analytic antenna pattern tabulated at 0.0025° |
| `compare_with_raytube.py` | compares two `ElectricField.csv` files receiver by receiver (and the number of paths per receiver), tolerance `|a - b| <= 1e-9 max|E| + 1e-7 |b|` |
| `run_comparison.py` | the driver `ctest` runs: convert, run the executable or take the committed output, compare |
| `cases.json` | per-scene settings: the reference's `num_tesselate` and the receivers compared |

Regenerate the committed reference output (Windows, or Linux with Wine) after
changing a scene:

    cmake --build build && ./build/sbrsolver SBRSolver_rev2/data/office.ini build/out_office
    python3 -I tools/raytube/run_comparison.py SBRSolver_rev2/data/office.ini \
        --sbrsolver-output build/out_office --work build/raytube_office \
        --exe tools/raytube/RayTubeTracing.exe --launcher wine64 \
        --save-reference SBRSolver_rev2/data/reference/raytube/office

## What is compared, and how exactly

Both programs end up computing the same thing for every path: the exact
specular path through the image sources, the source field
`sqrt(P eta0 / 4 pi) * sqrt(gain) * theta_hat`, Fresnel reflection and
transmission of the TE/TM components, and `exp(-j k L) / L` over the unfolded
length `L`. The remaining differences are floating-point noise (about 1e-12
relative); the reference prints 8 significant digits, so the comparison
demands agreement to one unit in its last digit (`--rtol 1e-7`).

The number of paths per receiver is compared as well, so a lost or spurious
path is reported even when its field is small.

## Mapping of the INI settings

| SBRSolver | RayTubeTracing |
|---|---|
| `[transmitter] coordinate_system` rows | `tx_axis`: `z` = row 2, `x` = row 0 |
| `power_watts` | `transmitting_power_dBm` |
| `antenna_type` | `tx_antenna_pattern_file` tabulated from the analytic pattern (`pattern_file` is ignored by both) |
| `[material]` (`rel_permittivity_re`, `rel_conductivity`) | `epsilon_setting_mode` 1 (static permittivity and conductivity), Fresnel reflection and transmission (type 1) |
| `use_fresnel_coefficients = false` (perfect conductors) | every material with static permittivity 1 and conductivity 1e30, the Fresnel limit of a perfect conductor |
| STL files | merged into one ASCII file whose coordinates are the float32 values SBRSolver reads, printed with 17 digits; both solvers recompute the facet normals from the vertices |
| per-file materials and `[surface_materials]` | `face_index` lists; the converter recomputes the coplanar ids and `--check-faces` verifies them against sbrsolver's `Faces.vtk` |
| `tessellation` n | `num_tesselate` 2^n, overridden per scene by `cases.json` |
| `[receivers]` line / plane / file | `receivers.dat` with the same points in the same order |

## Limits of the reference

RayTubeTracing launches ray *tubes* (triangles of a tessellated icosahedron).
A tube can be split only at a single boundary edge of the surface it hits; a
tube whose corner rays land on several surfaces, or whose sides cross none
or several edges, is dropped as a "complex hit". Near surface edges and
corners the reference therefore misses paths that exist, and finer tubes
only shrink the affected region. `cases.json` sets `num_tesselate` per scene
to a value at which the reference finds every path (checked against an
independent image-method calculation), and for `ground.ini`, whose ground
strip is only 10 m wide, limits the comparison to the first kilometre of the
receiver line: beyond that the tubes no longer fit on the strip at any
affordable tessellation and the reference loses the ground reflection.

SBRSolver's launch has the same kind of limit (a ray instead of a tube must
hit the surface close enough to the exact path), which is why `ground.ini`
and `office.ini` use `tessellation = 3` and `bahen.ini` `tessellation = 4`.

A reference run without any field or path at the compared receivers fails
the comparison, so a scene cannot pass by checking zeros against zeros.
