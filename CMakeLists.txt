cmake_minimum_required(VERSION 3.16)

project(SBRSolver
    VERSION 2.0
    DESCRIPTION "Shoot-Bounce-Reflect (SBR) electromagnetic ray tracer"
    LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(SBR_WARNINGS "Enable the compiler's standard warning set" ON)

# Default to a Release build for single-config generators (Makefiles, Ninja).
# Multi-config generators (Visual Studio, Xcode) pick the config at build time.
get_property(_sbr_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _sbr_multi_config AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Eigen (header only). Prefer the git submodule in ./eigen, otherwise fall
# back to an installed Eigen3 package (apt, brew, vcpkg, conda, ...).
# ---------------------------------------------------------------------------
set(_sbr_eigen_submodule "${CMAKE_CURRENT_SOURCE_DIR}/eigen")
if(EXISTS "${_sbr_eigen_submodule}/Eigen/Core")
    add_library(sbr_eigen INTERFACE)
    target_include_directories(sbr_eigen SYSTEM INTERFACE "${_sbr_eigen_submodule}")
    message(STATUS "SBRSolver: using Eigen from git submodule ${_sbr_eigen_submodule}")
else()
    find_package(Eigen3 3.4 CONFIG QUIET)
    if(Eigen3_FOUND)
        add_library(sbr_eigen INTERFACE)
        target_link_libraries(sbr_eigen INTERFACE Eigen3::Eigen)
        message(STATUS "SBRSolver: using installed Eigen3 ${Eigen3_VERSION}")
    else()
        message(FATAL_ERROR
            "Eigen was not found.\n"
            "Either fetch the bundled submodule:\n"
            "    git submodule update --init --recursive\n"
            "or install Eigen 3.4+ (apt install libeigen3-dev, brew install eigen, "
            "vcpkg install eigen3) and re-run CMake.")
    endif()
endif()

# ---------------------------------------------------------------------------
# Threads (the solver parallelises with std::thread, see common/ParallelFor.h)
# ---------------------------------------------------------------------------
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# FindThreads succeeds on every Windows toolchain, but libstdc++ only provides
# std::thread / std::mutex when built with gthreads (the MinGW-w64 "win32"
# thread model before GCC 13 has neither), so probe for them explicitly.
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES Threads::Threads)
check_cxx_source_compiles("
#include <mutex>
#include <thread>
int main() { std::mutex m; std::lock_guard<std::mutex> l(m); std::thread t([] {}); t.join(); return 0; }
" SBR_HAS_STD_THREAD)
unset(CMAKE_REQUIRED_LIBRARIES)
if(NOT SBR_HAS_STD_THREAD)
    message(FATAL_ERROR
        "SBRSolver: this toolchain's standard library has no std::thread/std::mutex. "
        "On Windows use MSVC, or a MinGW-w64 build with the posix or mcf thread model "
        "(for example the MSYS2 UCRT64/MINGW64 packages, or GCC 13+ for the win32 model).")
endif()

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
set(SBR_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/SBRSolver_rev2")

add_executable(sbrsolver
    ${SBR_SRC_DIR}/SBRSolver_rev2.cpp
    ${SBR_SRC_DIR}/FieldCompute.cpp
    ${SBR_SRC_DIR}/Icosahedron.cpp
    ${SBR_SRC_DIR}/Preprocessor.cpp
    ${SBR_SRC_DIR}/Ray.cpp
    ${SBR_SRC_DIR}/RTSolver.cpp
    ${SBR_SRC_DIR}/SimulationConfig.cpp
    ${SBR_SRC_DIR}/Triangle.cpp
    ${SBR_SRC_DIR}/timer/Timer.cpp
    # headers listed so IDE generators show them
    ${SBR_SRC_DIR}/FieldCompute.h
    ${SBR_SRC_DIR}/Icosahedron.h
    ${SBR_SRC_DIR}/Preprocessor.h
    ${SBR_SRC_DIR}/Ray.h
    ${SBR_SRC_DIR}/RTSolver.h
    ${SBR_SRC_DIR}/SimulationConfig.h
    ${SBR_SRC_DIR}/StlReader.h
    ${SBR_SRC_DIR}/Triangle.h
    ${SBR_SRC_DIR}/BVH/BVH.h
    ${SBR_SRC_DIR}/BVH/BVHNode.h
    ${SBR_SRC_DIR}/common/Constants.h
    ${SBR_SRC_DIR}/common/DataStructures.h
    ${SBR_SRC_DIR}/common/ParallelFor.h
    ${SBR_SRC_DIR}/common/VecMatDef.h
    ${SBR_SRC_DIR}/timer/Timer.h
)

target_include_directories(sbrsolver PRIVATE "${SBR_SRC_DIR}")
target_link_libraries(sbrsolver PRIVATE sbr_eigen Threads::Threads)

# Where the shipped configuration files, scenes and receiver files live, so the
# binary can be run from any working directory. The configuration file to use
# is the first command-line argument (default: <data>/ground.ini).
target_compile_definitions(sbrsolver PRIVATE
    SBR_DEFAULT_DATA_DIR="${SBR_SRC_DIR}/data")

if(WIN32)
    # Keep <windows.h> (if anything ever includes it) from defining min/max macros.
    target_compile_definitions(sbrsolver PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
endif()
if(MSVC)
    target_compile_definitions(sbrsolver PRIVATE _CRT_SECURE_NO_WARNINGS)
    if(SBR_WARNINGS)
        # /wd4267: size_t -> int narrowing, previously silenced with #pragma warning
        target_compile_options(sbrsolver PRIVATE /W4 /wd4267)
    endif()
elseif(SBR_WARNINGS)
    target_compile_options(sbrsolver PRIVATE -Wall -Wextra)
endif()

# std::filesystem lives in a separate library on older toolchains (libstdc++ 8:
# stdc++fs, libc++ 7/8: c++fs). Probe by linking instead of guessing from the
# compiler id, because Clang can be paired with either standard library.
set(_sbr_fs_src "#include <filesystem>\nint main() { return std::filesystem::exists(std::filesystem::current_path()) ? 0 : 1; }")
check_cxx_source_compiles("${_sbr_fs_src}" SBR_FILESYSTEM_BUILTIN)
if(NOT SBR_FILESYSTEM_BUILTIN)
    set(CMAKE_REQUIRED_LIBRARIES stdc++fs)
    check_cxx_source_compiles("${_sbr_fs_src}" SBR_FILESYSTEM_NEEDS_STDCXXFS)
    set(CMAKE_REQUIRED_LIBRARIES c++fs)
    check_cxx_source_compiles("${_sbr_fs_src}" SBR_FILESYSTEM_NEEDS_CXXFS)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(SBR_FILESYSTEM_NEEDS_STDCXXFS)
        target_link_libraries(sbrsolver PRIVATE stdc++fs)
    elseif(SBR_FILESYSTEM_NEEDS_CXXFS)
        target_link_libraries(sbrsolver PRIVATE c++fs)
    else()
        message(FATAL_ERROR
            "SBRSolver: this compiler/standard library has no usable std::filesystem. "
            "Needs GCC 8+, Clang 7+ with libstdc++ 8+ or libc++ 7+, MSVC 2017 15.7+, "
            "or Xcode 11+ targeting macOS 10.15 or newer.")
    endif()
endif()

# Visual Studio: start the debugger in the source tree like the old project did.
set_target_properties(sbrsolver PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${SBR_SRC_DIR}")

# ---------------------------------------------------------------------------
# Tests (ctest): run the reference configuration (ground.ini) from the build
# tree and, when a Python interpreter is available, compare the resulting field
# with the reference results committed under SBRSolver_rev2/data/output.
# ---------------------------------------------------------------------------
enable_testing()

set(_sbr_test_output "${CMAKE_CURRENT_BINARY_DIR}/test_output")
add_test(NAME sbrsolver_smoke
         COMMAND sbrsolver "${SBR_SRC_DIR}/data/ground.ini" "${_sbr_test_output}")
set_tests_properties(sbrsolver_smoke PROPERTIES TIMEOUT 600)

# A second scene exercising per-file and per-surface materials, compared with
# its reference field as well (data/reference/office_ElectricField.csv).
set(_sbr_test_output_office "${CMAKE_CURRENT_BINARY_DIR}/test_output_office")
add_test(NAME sbrsolver_office_smoke
         COMMAND sbrsolver "${SBR_SRC_DIR}/data/office.ini" "${_sbr_test_output_office}")
set_tests_properties(sbrsolver_office_smoke PROPERTIES TIMEOUT 600 FIXTURES_SETUP sbr_office_run)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME sbrsolver_reference_field
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_fields.py"
                     "${SBR_SRC_DIR}/data/output/ElectricField.csv"
                     "${_sbr_test_output}/ElectricField.csv")
    # Fixtures (not DEPENDS) so the smoke run always precedes the comparison,
    # even when only the comparison test is selected with ctest -R.
    set_tests_properties(sbrsolver_smoke PROPERTIES FIXTURES_SETUP sbr_default_run)
    set_tests_properties(sbrsolver_reference_field PROPERTIES FIXTURES_REQUIRED sbr_default_run)

    add_test(NAME sbrsolver_office_reference_field
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_fields.py"
                     "${SBR_SRC_DIR}/data/reference/office_ElectricField.csv"
                     "${_sbr_test_output_office}/ElectricField.csv")
    set_tests_properties(sbrsolver_office_reference_field PROPERTIES FIXTURES_REQUIRED sbr_office_run)
else()
    message(STATUS "SBRSolver: Python3 not found, skipping the reference field comparison test")
endif()
