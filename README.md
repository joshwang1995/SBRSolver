<a name="readme-top"></a>

# SBRSolver

<!-- ABOUT THE PROJECT -->
## About The Project

Implementation of the Shoot Bounce Relfect (SBR) method for electromagnetics simulation.

Features
* Transmitter Modelling
  - Ray generation via icosahedron sphere tessellation
  - Choose and apply theoretical antenna patterns (Hertzian Dipole, half-wave dipole, isotropic)
* Receiver Modelling
  - Reception sphere algorithm
  - Received field computation
* Scene Modelling
  - STL to triangle primitives
  - Bounding Volume Hierachy (BVH) implementation
  - Coplanar surface identification
* Ray Path Tracing
  - Ray math library
  - Ray launching algorithm
  - Ray intersection tests
  - Duplicate path removal
  - Image method path correction
* Field Calculation
  - TE/TM field decomposition
  - TE/TM field reflection
  - TE/TM field transmission

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

SBRSolver is a plain C++17 program built with [CMake](https://cmake.org). It compiles and runs
on Windows, Linux and macOS with any recent compiler (MSVC 2019+, GCC 8+, Clang 9+, Xcode 11+ on macOS 10.15 or newer).
The only dependency, [Eigen](https://eigen.tuxfamily.org), is bundled as a git submodule.

### Prerequisites

* Git
* CMake 3.16 or newer
* A C++17 compiler
  - Windows: Visual Studio 2019/2022 with the *Desktop development with C++* workload
    (the command-line *Build Tools* are enough), or MinGW-w64/Clang
  - Linux: `sudo apt install build-essential cmake git` (Debian/Ubuntu) or the equivalent
  - macOS 10.15 (Catalina) or newer: `xcode-select --install` and `brew install cmake`
* Optional: Python 3 to run the regression test

### Build

```sh
git clone --recurse-submodules https://github.com/joshwang1995/SBRSolver.git
cd SBRSolver
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
```

If you cloned without `--recurse-submodules`, fetch Eigen first with
`git submodule update --init --recursive`. Alternatively install Eigen 3.4 system-wide
(`apt install libeigen3-dev`, `brew install eigen`, `vcpkg install eigen3`) and CMake will pick it up.

The executable is written to `build/sbrsolver` on Linux/macOS. On Windows the Visual Studio generator
(CMake's default when Visual Studio is installed) puts it in `build\Release\sbrsolver.exe`; single-config
generators such as Ninja or MinGW Makefiles put it in `build\sbrsolver.exe`.
On Windows you can also open the folder in Visual Studio (*File | Open | Folder*) or generate a solution
with `cmake -S . -B build -G "Visual Studio 17 2022"`; both use the same `CMakeLists.txt`.

### Run

```sh
./build/sbrsolver                                                # runs SBRSolver_rev2/data/ground.ini
./build/sbrsolver SBRSolver_rev2/data/office.ini results/office  # another scene, into its own output directory
./build/sbrsolver my_case.ini                                    # custom config, results in output/ next to it
./build/sbrsolver --help
```

The first argument is a configuration file (see below), the second the output directory for the
VTK/CSV results. The configuration defaults to `ground.ini` in the data folder of the source tree
(the path is baked in at build time, so the program can be started from any working directory);
the output directory defaults to `output/` next to the configuration file and is created if missing.
`SBRSolver_rev2/data/output` holds the committed reference results that `ctest` compares against, so run
scenes other than `ground.ini` with an explicit output directory.

### Configuration file

All simulation settings live in an INI-style text file. `SBRSolver_rev2/data/ground.ini` is the
documented reference; `corner.ini`, `office.ini` and `bahen.ini` are further examples. File names
inside the configuration are relative to the configuration file. Lines starting with `#` or `;` are
comments, and `#` preceded by a space starts an inline comment. Numbers in a value are separated by
spaces (matrix rows by `;`) and always use `.` as the decimal mark. Unknown sections or keys, missing
required keys, values of the wrong type or range and missing files are reported with their line number.

| Section | Key | Meaning |
|---|---|---|
| `[scene]` | `stl_file` | scene geometry (STL, ASCII or binary) |
| | `save_edges`, `save_faces` | write `Edges.vtk` / `Faces.vtk` (default `true`) |
| `[transmitter]` | `position` | ray origin, `x y z` in metres |
| | `coordinate_system` | three rows of the TX frame, e.g. `0 1 0; 0 0 1; 1 0 0` (default identity) |
| | `power_watts`, `frequency_hz` | transmit power (default `1`) and frequency |
| | `tessellation` | icosahedron subdivision level 0 to 8; launches 10·4ⁿ+2 rays (default `2`, 162 rays) |
| | `antenna_type` | `isotropic`, `hertzian_dipole` (default) or `half_wave_dipole` |
| | `pattern_file` | optional measured pattern; read and validated, not applied yet |
| `[receivers]` | `layout = line` | `line_start`, `line_end`, `resolution` (generated layouts are limited to 100 million receiver × ray combinations, about 3 GB) |
| | `layout = plane` | `x_min`, `y_min`, `x_max`, `y_max`, `height`, `resolution` |
| | `layout = file` | `file`, an `RX_*.dat` location file |
| `[ray_tracing]` | `max_reflections`, `max_transmissions` | bounce limits (default `1` and `0`) |
| | `use_fresnel_coefficients` | default `true` |
| `[material]` | `name`, `reflection_loss`, `transmission_loss`, `rel_conductivity`, `rel_permittivity_re`, `rel_permittivity_im`, `frequency_hz` | one section per material; the first is material 0, applied to every surface. `frequency_hz` defaults to the transmitter frequency |
| `[outputs]` | `ray_paths`, `receivers`, `field_csv`, `icosahedron` | which result files to write (defaults `true`, `true`, `true`, `false`) |

The scene, receiver layout and material parameters are no longer set in the source code.

### Test

```sh
cd build && ctest -C Release --output-on-failure
```

This runs `ground.ini` and compares the resulting electric field with the reference results
committed in `SBRSolver_rev2/data/output/ElectricField.csv`. The same build and test matrix runs in
GitHub Actions on Ubuntu (GCC and Clang), macOS and Windows for every pull request and every push to `master`.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Contributing

**_Git is used for version control in this entire project. Any use of Git without fully understanding the consequence of each action can be detrimental. Read the following resource to get started with using Git: https://git-scm.com/book/en/v2_**

1. Do not commit directly to the `master` branch. It should only contain verified working code.
2. Create your own branch from `master`: `git switch -c my-feature master`
3. Commit frequently so you can go back if something goes wrong. Think of it as quick save in video games:
   `git add -A && git commit -m "Short, descriptive message"`
4. Push your branch to GitHub: `git push -u origin my-feature`
5. Once your fix/implementation is finished, open a pull request to merge it into `master`.
   **Make sure the build and `ctest` pass on your machine, and get someone else to review it for you.**

Any Git client (command line, Visual Studio, VS Code, GitKraken, ...) works for the steps above.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- TROUBLESHOOTING -->
## Troubleshooting

* **CMake reports that Eigen was not found.** The submodule folder `eigen/` is empty. Run
  ```sh
  git submodule update --init --recursive
  ```
  in the repository root and configure again.
* **`Error: invalid configuration`** at start-up. The message lists each offending line of the
  configuration file (unknown key, bad value, missing file). File names inside the configuration are
  relative to the configuration file, not to the current directory.
* **OpenMP.** Nothing in the solver uses OpenMP yet, so it is off by default. To link it (for example
  before re-enabling the `#pragma omp` loops in `RTSolver.cpp`) configure with `-DSBR_ENABLE_OPENMP=ON`;
  on macOS install it first with `brew install libomp` and add `-DOpenMP_ROOT=$(brew --prefix libomp)`.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
