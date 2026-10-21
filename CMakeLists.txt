cmake_minimum_required(VERSION 3.16)

project(SBRSolver
    VERSION 2.0
    DESCRIPTION "Shoot-Bounce-Reflect (SBR) electromagnetic ray tracer"
    LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(SBR_WARNINGS "Enable the compiler's standard warning set" ON)
option(SBR_BUILD_VIEWER "Build the built-in 3D viewer (Polyscope, --view / [outputs] viewer = true)" ON)
option(SBR_ENABLE_TRACY "Instrument the solver for the Tracy profiler (connect the Tracy GUI while it runs)" OFF)
option(SBR_TRACY_ON_DEMAND "Tracy records only once the GUI connects; OFF records from start-up (run with TRACY_NO_EXIT=1 to wait for the GUI at exit)" ON)

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
# Third-party libraries (git submodules under third_party/)
# ---------------------------------------------------------------------------
set(_sbr_third_party "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

if(SBR_BUILD_VIEWER)
    if(NOT EXISTS "${_sbr_third_party}/polyscope/CMakeLists.txt" OR NOT EXISTS "${_sbr_third_party}/polyscope/deps/glfw/CMakeLists.txt")
        message(FATAL_ERROR
            "SBRSolver: the viewer needs the Polyscope submodule and its own submodules. Run\n"
            "    git submodule update --init --recursive\n"
            "or configure with -DSBR_BUILD_VIEWER=OFF.")
    endif()
    # Polyscope bundles GLFW, Dear ImGui, glm and stb. The mock backend lets the
    # viewer run without a display (used by the tests); EGL is picked up on Linux
    # when available for headless rendering.
    set(POLYSCOPE_BACKEND_OPENGL3_GLFW ON CACHE BOOL "" FORCE)
    set(POLYSCOPE_BACKEND_OPENGL_MOCK ON CACHE BOOL "" FORCE)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    add_subdirectory("${_sbr_third_party}/polyscope" EXCLUDE_FROM_ALL)
    # Polyscope's headers are not warning-clean under -Wextra; treat them as
    # system headers so our own warnings stay visible.
    foreach(_sbr_dep polyscope imgui glm glfw)
        if(TARGET ${_sbr_dep})
            set_target_properties(${_sbr_dep} PROPERTIES
                INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "$<TARGET_PROPERTY:${_sbr_dep},INTERFACE_INCLUDE_DIRECTORIES>")
        endif()
    endforeach()
    message(STATUS "SBRSolver: viewer enabled (Polyscope)")
endif()

if(SBR_ENABLE_TRACY)
    if(NOT EXISTS "${_sbr_third_party}/tracy/CMakeLists.txt")
        message(FATAL_ERROR
            "SBRSolver: the Tracy submodule is missing. Run\n"
            "    git submodule update --init --recursive\n"
            "or configure with -DSBR_ENABLE_TRACY=OFF.")
    endif()
    # The client records once the Tracy GUI connects (on demand). Until then
    # the zones cost nothing measurable, but the client itself adds about
    # 0.2 s of start-up (timer calibration), a few threads and a listening
    # socket (TCP 8086) to every run, which is why the option is off by
    # default and CI builds with it on.
    set(TRACY_ENABLE ON CACHE BOOL "" FORCE)
    set(TRACY_ON_DEMAND ${SBR_TRACY_ON_DEMAND} CACHE BOOL "" FORCE)
    set(TRACY_STATIC ON CACHE BOOL "" FORCE)
    add_subdirectory("${_sbr_third_party}/tracy" EXCLUDE_FROM_ALL)
    message(STATUS "SBRSolver: Tracy profiling enabled (on demand)")
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
    ${SBR_SRC_DIR}/viewer/SceneData.cpp
    ${SBR_SRC_DIR}/viewer/FieldGrid.cpp
    # headers listed so IDE generators show them
    ${SBR_SRC_DIR}/viewer/SceneData.h
    ${SBR_SRC_DIR}/viewer/FieldGrid.h
    ${SBR_SRC_DIR}/viewer/Viewer.h
    ${SBR_SRC_DIR}/common/Profiling.h
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

if(SBR_BUILD_VIEWER)
    target_sources(sbrsolver PRIVATE ${SBR_SRC_DIR}/viewer/Viewer.cpp)
    target_link_libraries(sbrsolver PRIVATE polyscope)
    target_compile_definitions(sbrsolver PRIVATE SBR_HAS_VIEWER=1)
endif()
if(SBR_ENABLE_TRACY)
    target_link_libraries(sbrsolver PRIVATE Tracy::TracyClient)
    target_compile_definitions(sbrsolver PRIVATE SBR_HAS_TRACY=1)
endif()

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
        set(_sbr_fs_lib stdc++fs)
    elseif(SBR_FILESYSTEM_NEEDS_CXXFS)
        set(_sbr_fs_lib c++fs)
    else()
        message(FATAL_ERROR
            "SBRSolver: this compiler/standard library has no usable std::filesystem. "
            "Needs GCC 8+, Clang 7+ with libstdc++ 8+ or libc++ 7+, MSVC 2017 15.7+, "
            "or Xcode 11+ targeting macOS 10.15 or newer.")
    endif()
    target_link_libraries(sbrsolver PRIVATE ${_sbr_fs_lib})
endif()

# Visual Studio: start the debugger in the source tree like the old project did.
set_target_properties(sbrsolver PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${SBR_SRC_DIR}")

# ---------------------------------------------------------------------------
# Tests (ctest): run the four scenes in SBRSolver_rev2/data from the build tree
# and, when a Python interpreter is available, compare the fields with the
# committed self-references (data/output, data/reference) and with the
# RayTubeTracing reference solver (tools/raytube).
# ---------------------------------------------------------------------------
enable_testing()

# Every scene is run once (test_output_<scene>), then compared with its
# committed self-reference (cross-platform determinism, tools/compare_fields.py)
# and with the RayTubeTracing reference solver (physics, tools/raytube).
set(_sbr_scenes ground corner office bahen)
foreach(_scene IN LISTS _sbr_scenes)
    add_test(NAME sbrsolver_run_${_scene}
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/${_scene}.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}")
    # Fixtures (not DEPENDS) so the run always precedes its comparisons, even
    # when only a comparison test is selected with ctest -R.
    set_tests_properties(sbrsolver_run_${_scene} PROPERTIES TIMEOUT 600 FIXTURES_SETUP sbr_run_${_scene})
endforeach()

# Unit tests of the viewer's data model and field interpolation (no display).
add_executable(sbr_viewer_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/viewer_tests.cpp
    ${SBR_SRC_DIR}/viewer/SceneData.cpp
    ${SBR_SRC_DIR}/viewer/FieldGrid.cpp
    ${SBR_SRC_DIR}/Triangle.cpp)
target_include_directories(sbr_viewer_tests PRIVATE "${SBR_SRC_DIR}")
target_link_libraries(sbr_viewer_tests PRIVATE sbr_eigen)
if(MSVC)
    target_compile_definitions(sbr_viewer_tests PRIVATE _CRT_SECURE_NO_WARNINGS)
    if(SBR_WARNINGS)
        target_compile_options(sbr_viewer_tests PRIVATE /W4 /wd4267)
    endif()
elseif(SBR_WARNINGS)
    target_compile_options(sbr_viewer_tests PRIVATE -Wall -Wextra)
endif()
add_test(NAME viewer_unit_tests COMMAND sbr_viewer_tests)

# Smoke test of the viewer itself: a synthetic run goes through every overlay
# with Polyscope's mock backend, so it needs no display. With a display (or
# xvfb-run) `sbr_viewer_smoke --screenshot shot.png` renders it for real.
if(SBR_BUILD_VIEWER)
    add_executable(sbr_viewer_smoke
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/viewer_smoke.cpp
        ${SBR_SRC_DIR}/viewer/Viewer.cpp
        ${SBR_SRC_DIR}/viewer/SceneData.cpp
        ${SBR_SRC_DIR}/viewer/FieldGrid.cpp
        ${SBR_SRC_DIR}/Triangle.cpp)
    target_include_directories(sbr_viewer_smoke PRIVATE "${SBR_SRC_DIR}")
    target_link_libraries(sbr_viewer_smoke PRIVATE sbr_eigen polyscope ${_sbr_fs_lib})   # Viewer.cpp uses std::filesystem
    if(MSVC)
        target_compile_definitions(sbr_viewer_smoke PRIVATE _CRT_SECURE_NO_WARNINGS)
    elseif(SBR_WARNINGS)
        target_compile_options(sbr_viewer_smoke PRIVATE -Wall -Wextra)
    endif()
    add_test(NAME viewer_smoke_mock COMMAND sbr_viewer_smoke --backend openGL_mock --frames 2 --runs 2)
    # Three receiver layers: the field volume (hex mesh) and the slice plane exist.
    add_test(NAME viewer_smoke_mock_volume COMMAND sbr_viewer_smoke --backend openGL_mock --frames 2 --receiver-layers 3)
    set_tests_properties(viewer_smoke_mock viewer_smoke_mock_volume PROPERTIES TIMEOUT 300)
endif()

if(SBR_BUILD_VIEWER)
    # The viewer on a finished run without any rendering (mock backend): every
    # overlay is registered and two frames are ticked. Runs on every platform.
    add_test(NAME viewer_mock_office
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_viewer"
                     --viewer-backend openGL_mock --viewer-frames 2)
    set_tests_properties(viewer_mock_office PROPERTIES TIMEOUT 600)
    # Malformed viewer options are rejected before the run starts.
    add_test(NAME viewer_cli_rejects_frames
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_viewer"
                     --viewer-backend openGL_mock --viewer-frames 2x)
    add_test(NAME viewer_cli_rejects_screenshot_type
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_viewer"
                     --viewer-backend openGL_mock --screenshot "${CMAKE_CURRENT_BINARY_DIR}/test_output_viewer/office.txt")
    set_tests_properties(viewer_cli_rejects_frames viewer_cli_rejects_screenshot_type PROPERTIES WILL_FAIL TRUE TIMEOUT 60)
    # A real rendering into a PNG where a display (or Xvfb) is available.
    find_program(SBR_XVFB_RUN xvfb-run)
    if(SBR_XVFB_RUN)
        add_test(NAME viewer_screenshot_office
                 COMMAND "${SBR_XVFB_RUN}" -a "$<TARGET_FILE:sbrsolver>" "${SBR_SRC_DIR}/data/office.ini"
                         "${CMAKE_CURRENT_BINARY_DIR}/test_output_viewer_shot"
                         --screenshot "${CMAKE_CURRENT_BINARY_DIR}/test_output_viewer_shot/office.png")
        # The viewer removes a stale PNG first and fails when none was written.
        set_tests_properties(viewer_screenshot_office PROPERTIES TIMEOUT 600)
    endif()
endif()

# The comparisons need a Python interpreter. They are skipped without one,
# unless SBR_REQUIRE_REFERENCE_TESTS is set (as CI does), which makes the
# missing interpreter a configure error instead of a silently smaller suite.
option(SBR_REQUIRE_REFERENCE_TESTS "Fail to configure when the reference comparison tests cannot be registered" OFF)
if(SBR_REQUIRE_REFERENCE_TESTS)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
else()
    find_package(Python3 COMPONENTS Interpreter QUIET)
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME sbrsolver_reference_ground
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_fields.py"
                     "${SBR_SRC_DIR}/data/output/ElectricField.csv"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_ground/ElectricField.csv")
    set_tests_properties(sbrsolver_reference_ground PROPERTIES FIXTURES_REQUIRED sbr_run_ground)
    add_test(NAME sbrsolver_reference_office
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_fields.py"
                     "${SBR_SRC_DIR}/data/reference/office_ElectricField.csv"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office/ElectricField.csv")
    set_tests_properties(sbrsolver_reference_office PROPERTIES FIXTURES_REQUIRED sbr_run_office)

    # Comparison with RayTubeTracing. With SBR_RAYTUBE_EXE the reference solver
    # is run (on Windows the committed executable is used by default); without
    # it the committed reference output under data/reference/raytube is used.
    set(_sbr_raytube_default "")
    if(WIN32 AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/raytube/RayTubeTracing.exe")
        set(_sbr_raytube_default "${CMAKE_CURRENT_SOURCE_DIR}/tools/raytube/RayTubeTracing.exe")
    endif()
    set(SBR_RAYTUBE_EXE "${_sbr_raytube_default}" CACHE FILEPATH
        "RayTubeTracing.exe to run for the raytube_* tests (empty: compare with the committed reference output)")
    set(SBR_RAYTUBE_LAUNCHER "" CACHE STRING
        "Program that runs SBR_RAYTUBE_EXE on this platform, e.g. wine64 (empty: run it directly)")
    string(REPLACE "\\" "/" SBR_RAYTUBE_EXE "${SBR_RAYTUBE_EXE}")
    if(SBR_RAYTUBE_EXE AND NOT EXISTS "${SBR_RAYTUBE_EXE}")
        message(FATAL_ERROR "SBRSolver: SBR_RAYTUBE_EXE does not exist: ${SBR_RAYTUBE_EXE}")
    endif()
    foreach(_scene IN LISTS _sbr_scenes)
        set(_args "${SBR_SRC_DIR}/data/${_scene}.ini"
                  --sbrsolver-output "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}"
                  --work "${CMAKE_CURRENT_BINARY_DIR}/raytube_${_scene}"
                  --reference "${SBR_SRC_DIR}/data/reference/raytube/${_scene}")
        if(SBR_RAYTUBE_EXE)
            list(APPEND _args --exe "${SBR_RAYTUBE_EXE}")
            if(SBR_RAYTUBE_LAUNCHER)
                list(APPEND _args --launcher "${SBR_RAYTUBE_LAUNCHER}")
            endif()
        endif()
        add_test(NAME raytube_${_scene}
                 COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/raytube/run_comparison.py" ${_args})
        set_tests_properties(raytube_${_scene} PROPERTIES FIXTURES_REQUIRED sbr_run_${_scene} TIMEOUT 1800)
    endforeach()
    if(SBR_RAYTUBE_EXE)
        message(STATUS "SBRSolver: raytube_* tests run ${SBR_RAYTUBE_EXE}")
    else()
        message(STATUS "SBRSolver: raytube_* tests compare with the committed RayTubeTracing output (set SBR_RAYTUBE_EXE to run it)")
    endif()
else()
    message(STATUS "SBRSolver: Python3 not found, skipping the reference comparison tests")
endif()
