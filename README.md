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
./build/sbrsolver                        # uses SBRSolver_rev2/data, writes to SBRSolver_rev2/data/output
./build/sbrsolver path/to/data out_dir   # custom data and output directories
./build/sbrsolver --help
```

The first argument is the data directory (containing `stl_files/` and the `RX_*.dat` receiver files),
the second the output directory for the VTK/CSV results. The data directory defaults to `SBRSolver_rev2/data`
in the source tree (the path is baked in at build time, so the program can be started from any working
directory); the output directory defaults to `<data_dir>/output` and is created if missing. The scene, receiver layout and material
parameters are currently set in `SBRSolver_rev2/SBRSolver_rev2.cpp`.

### Test

```sh
cd build && ctest -C Release --output-on-failure
```

This runs the default scene and compares the resulting electric field with the reference results
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
* **`STL file not found`** at start-up. Pass the data directory explicitly:
  `sbrsolver /path/to/SBRSolver/SBRSolver_rev2/data`.
* **OpenMP.** Nothing in the solver uses OpenMP yet, so it is off by default. To link it (for example
  before re-enabling the `#pragma omp` loops in `RTSolver.cpp`) configure with `-DSBR_ENABLE_OPENMP=ON`;
  on macOS install it first with `brew install libomp` and add `-DOpenMP_ROOT=$(brew --prefix libomp)`.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
