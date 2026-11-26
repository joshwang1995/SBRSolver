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
# The GPU tracer runs on Vulkan. macOS has no Vulkan driver with the
# features it needs, so the option is off there.
if(APPLE)
    set(_sbr_gpu_default OFF)
else()
    set(_sbr_gpu_default ON)
endif()
option(SBR_BUILD_GPU "Build the Vulkan GPU tracer (--device gpu)" ${_sbr_gpu_default})
option(SBR_BUILD_GLSLANG "Build the bundled glslang shader compiler even when one is installed" OFF)

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
    # By default the client records once the Tracy GUI connects (on demand);
    # with SBR_TRACY_ON_DEMAND=OFF it records from start-up (and waits for
    # the GUI at exit when TRACY_NO_EXIT=1 is set), which suits short runs.
    # Until the GUI connects the zones cost nothing measurable, but the
    # client itself adds about 0.2 s of start-up (timer calibration), a few
    # threads and a listening socket (TCP 8086) to every run, which is why
    # SBR_ENABLE_TRACY is off by default and CI builds with it on.
    set(TRACY_ENABLE ON CACHE BOOL "" FORCE)
    set(TRACY_ON_DEMAND ${SBR_TRACY_ON_DEMAND} CACHE BOOL "" FORCE)
    set(TRACY_STATIC ON CACHE BOOL "" FORCE)
    add_subdirectory("${_sbr_third_party}/tracy" EXCLUDE_FROM_ALL)
    if(SBR_TRACY_ON_DEMAND)
        message(STATUS "SBRSolver: Tracy profiling enabled (on demand)")
    else()
        message(STATUS "SBRSolver: Tracy profiling enabled (recording from start-up; run with TRACY_NO_EXIT=1 to wait for the GUI at exit)")
    endif()
endif()

if(SBR_BUILD_GPU)
    foreach(_sbr_dep Vulkan-Headers volk glslang)
        if(NOT EXISTS "${_sbr_third_party}/${_sbr_dep}/CMakeLists.txt")
            message(FATAL_ERROR
                "SBRSolver: the GPU tracer needs the ${_sbr_dep} submodule. Run\n"
                "    git submodule update --init --recursive\n"
                "or configure with -DSBR_BUILD_GPU=OFF.")
        endif()
    endforeach()
    # The shaders are compiled to SPIR-V at build time, by an installed
    # glslang (apt install glslang-tools, or the Vulkan SDK) or, without one,
    # by the bundled glslang, which takes a few minutes to build.
    if(NOT SBR_BUILD_GLSLANG)
        find_program(SBR_GLSLANG NAMES glslang glslangValidator)
    endif()
    if(SBR_GLSLANG AND NOT SBR_BUILD_GLSLANG)
        set(_sbr_glslang_command "${SBR_GLSLANG}")
        set(_sbr_glslang_target "")
        message(STATUS "SBRSolver: GPU tracer enabled, shaders compiled with ${SBR_GLSLANG}")
    else()
        set(ENABLE_OPT OFF CACHE BOOL "" FORCE)
        set(ENABLE_HLSL OFF CACHE BOOL "" FORCE)
        set(ENABLE_GLSLANG_BINARIES ON CACHE BOOL "" FORCE)
        set(BUILD_EXTERNAL OFF CACHE BOOL "" FORCE)
        set(GLSLANG_TESTS OFF CACHE BOOL "" FORCE)
        set(GLSLANG_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
        add_subdirectory("${_sbr_third_party}/glslang" EXCLUDE_FROM_ALL)
        set(_sbr_glslang_command "$<TARGET_FILE:glslang-standalone>")
        set(_sbr_glslang_target glslang-standalone)
        message(STATUS "SBRSolver: GPU tracer enabled, shaders compiled with the bundled glslang")
    endif()
    # The Vulkan headers and the volk loader. volk loads the driver's
    # library at run time, so the solver starts without one and falls back
    # to the CPU tracer.
    add_library(sbr_vulkan INTERFACE)
    target_include_directories(sbr_vulkan SYSTEM INTERFACE
        "${_sbr_third_party}/Vulkan-Headers/include" "${_sbr_third_party}/volk")
    target_compile_definitions(sbr_vulkan INTERFACE VK_NO_PROTOTYPES)
    # Every shader is compiled twice, in single and in double precision, and
    # embedded in the executable as a header (cmake/EmbedSpirv.cmake).
    set(SBR_SHADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/SBRSolver_rev2/gpu/shaders")
    set(SBR_SHADER_OUT "${CMAKE_CURRENT_BINARY_DIR}/shaders")
    file(MAKE_DIRECTORY "${SBR_SHADER_OUT}")
    set(SBR_SHADER_HEADERS "")
    function(sbr_add_shader name)
        foreach(_precision single double)
            if(_precision STREQUAL "double")
                set(_double 1)
            else()
                set(_double 0)
            endif()
            set(_spv "${SBR_SHADER_OUT}/${name}_${_precision}.spv")
            set(_header "${SBR_SHADER_OUT}/${name}_${_precision}.spv.h")
            add_custom_command(OUTPUT "${_spv}"
                COMMAND ${_sbr_glslang_command} -V --target-env vulkan1.1 "-I${SBR_SHADER_DIR}"
                        -DSBR_DOUBLE=${_double} -o "${_spv}" "${SBR_SHADER_DIR}/${name}.comp"
                DEPENDS "${SBR_SHADER_DIR}/${name}.comp" "${SBR_SHADER_DIR}/common.glsl" ${_sbr_glslang_target}
                COMMENT "Compiling shader ${name} (${_precision})")
            add_custom_command(OUTPUT "${_header}"
                COMMAND ${CMAKE_COMMAND} "-DINPUT=${_spv}" "-DOUTPUT=${_header}" "-DNAME=${name}_${_precision}"
                        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedSpirv.cmake"
                DEPENDS "${_spv}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedSpirv.cmake"
                COMMENT "Embedding shader ${name} (${_precision})")
            list(APPEND SBR_SHADER_HEADERS "${_header}")
        endforeach()
        set(SBR_SHADER_HEADERS "${SBR_SHADER_HEADERS}" PARENT_SCOPE)
    endfunction()
    sbr_add_shader(probe)
    sbr_add_shader(trace)
    sbr_add_shader(capture)
    add_custom_target(sbr_shaders DEPENDS ${SBR_SHADER_HEADERS})
endif()

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
set(SBR_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/SBRSolver_rev2")

add_executable(sbrsolver
    ${SBR_SRC_DIR}/SBRSolver_rev2.cpp
    ${SBR_SRC_DIR}/CpuTracer.cpp
    ${SBR_SRC_DIR}/FieldCompute.cpp
    ${SBR_SRC_DIR}/Icosahedron.cpp
    ${SBR_SRC_DIR}/Preprocessor.cpp
    ${SBR_SRC_DIR}/ReceiverGrid.cpp
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
    ${SBR_SRC_DIR}/gpu/GpuTypes.h
    ${SBR_SRC_DIR}/CpuTracer.h
    ${SBR_SRC_DIR}/FieldCompute.h
    ${SBR_SRC_DIR}/Icosahedron.h
    ${SBR_SRC_DIR}/Preprocessor.h
    ${SBR_SRC_DIR}/Ray.h
    ${SBR_SRC_DIR}/ReceiverGrid.h
    ${SBR_SRC_DIR}/RTSolver.h
    ${SBR_SRC_DIR}/SimulationConfig.h
    ${SBR_SRC_DIR}/StlReader.h
    ${SBR_SRC_DIR}/Tracer.h
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
if(SBR_BUILD_GPU)
    target_sources(sbrsolver PRIVATE
        ${SBR_SRC_DIR}/gpu/GpuContext.cpp
        ${SBR_SRC_DIR}/gpu/GpuTracer.cpp
        ${SBR_SRC_DIR}/gpu/VolkImpl.cpp
        ${SBR_SRC_DIR}/gpu/GpuContext.h
        ${SBR_SRC_DIR}/gpu/GpuRecords.h
        ${SBR_SRC_DIR}/gpu/GpuTracer.h)
    target_include_directories(sbrsolver PRIVATE "${SBR_SHADER_OUT}")
    target_link_libraries(sbrsolver PRIVATE sbr_vulkan ${CMAKE_DL_LIBS})
    target_compile_definitions(sbrsolver PRIVATE SBR_HAS_GPU=1)
    add_dependencies(sbrsolver sbr_shaders)
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
    # Every run writes Paths.csv for the GPU comparisons; the office run
    # also writes Segments.csv for the comparison tool tests.
    set(_sbr_run_options --paths-csv)
    if(_scene STREQUAL "office")
        set(_sbr_run_options --paths-csv --segments-csv)
    endif()
    add_test(NAME sbrsolver_run_${_scene}
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/${_scene}.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}" ${_sbr_run_options})
    # Fixtures (not DEPENDS) so the run always precedes its comparisons, even
    # when only a comparison test is selected with ctest -R.
    set_tests_properties(sbrsolver_run_${_scene} PROPERTIES TIMEOUT 600 FIXTURES_SETUP sbr_run_${_scene})
endforeach()

# The dense receiver layouts of data/fine run as well; the field file must
# then hold one finite field per receiver (checked by tools/check_fields.py
# when a Python interpreter is available).
foreach(_scene IN LISTS _sbr_scenes)
    add_test(NAME sbrsolver_run_${_scene}_fine
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/fine/${_scene}.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_fine" --paths-csv)
    set_tests_properties(sbrsolver_run_${_scene}_fine PROPERTIES TIMEOUT 600 FIXTURES_SETUP sbr_run_${_scene}_fine)
endforeach()

# The large-office plan scene: both access points radiate at once.
add_test(NAME sbrsolver_run_large_office
         COMMAND sbrsolver "${SBR_SRC_DIR}/data/large_office/large_office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_large_office" --paths-csv)
set_tests_properties(sbrsolver_run_large_office PROPERTIES TIMEOUT 1800 FIXTURES_SETUP sbr_run_large_office)

# The office scene with two transmitters radiating at once.
add_test(NAME sbrsolver_run_office_two_tx
         COMMAND sbrsolver "${SBR_SRC_DIR}/data/office_two_tx.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx" --paths-csv --segments-csv)
set_tests_properties(sbrsolver_run_office_two_tx PROPERTIES TIMEOUT 600 FIXTURES_SETUP sbr_run_office_two_tx)
add_test(NAME sbrsolver_run_office_two_tx_power
         COMMAND sbrsolver "${SBR_SRC_DIR}/data/office_two_tx_power.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx_power")
set_tests_properties(sbrsolver_run_office_two_tx_power PROPERTIES TIMEOUT 600 FIXTURES_SETUP sbr_run_office_two_tx_power)

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

# Unit tests of the capture stage: the receiver grid must find the same
# candidates as a test of every receiver against every segment.
add_executable(sbr_tracer_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/tracer_tests.cpp
    ${SBR_SRC_DIR}/CpuTracer.cpp
    ${SBR_SRC_DIR}/ReceiverGrid.cpp
    ${SBR_SRC_DIR}/Triangle.cpp)
target_include_directories(sbr_tracer_tests PRIVATE "${SBR_SRC_DIR}")
target_link_libraries(sbr_tracer_tests PRIVATE sbr_eigen Threads::Threads)
if(MSVC)
    target_compile_definitions(sbr_tracer_tests PRIVATE _CRT_SECURE_NO_WARNINGS)
    if(SBR_WARNINGS)
        target_compile_options(sbr_tracer_tests PRIVATE /W4 /wd4267)
    endif()
elseif(SBR_WARNINGS)
    target_compile_options(sbr_tracer_tests PRIVATE -Wall -Wextra)
endif()
add_test(NAME tracer_unit_tests COMMAND sbr_tracer_tests)

# Unit tests of the Vulkan host layer with the probe shader. Without a
# Vulkan device the test exits with 77 and ctest reports it as skipped.
if(SBR_BUILD_GPU)
    add_executable(sbr_gpu_tests
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/gpu_tests.cpp
        ${SBR_SRC_DIR}/gpu/GpuContext.cpp
        ${SBR_SRC_DIR}/gpu/VolkImpl.cpp)
    target_include_directories(sbr_gpu_tests PRIVATE "${SBR_SRC_DIR}" "${SBR_SHADER_OUT}")
    target_link_libraries(sbr_gpu_tests PRIVATE sbr_vulkan ${CMAKE_DL_LIBS})
    add_dependencies(sbr_gpu_tests sbr_shaders)
    if(WIN32)
        target_compile_definitions(sbr_gpu_tests PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    if(MSVC)
        target_compile_definitions(sbr_gpu_tests PRIVATE _CRT_SECURE_NO_WARNINGS)
        if(SBR_WARNINGS)
            target_compile_options(sbr_gpu_tests PRIVATE /W4 /wd4267)
        endif()
    elseif(SBR_WARNINGS)
        target_compile_options(sbr_gpu_tests PRIVATE -Wall -Wextra)
    endif()
    add_test(NAME gpu_unit_tests COMMAND sbr_gpu_tests)
    set_tests_properties(gpu_unit_tests PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 300)

    # The scenes on the GPU tracer, in both precisions, each compared with
    # the CPU run of the same scene below. Without a Vulkan device the runs
    # exit with 77 and every test that depends on them skips too.
    set(_sbr_gpu_scenes ground corner office bahen office_two_tx)
    foreach(_scene IN LISTS _sbr_gpu_scenes)
        foreach(_precision double single)
            add_test(NAME gpu_run_${_scene}_${_precision}
                     COMMAND sbrsolver "${SBR_SRC_DIR}/data/${_scene}.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_gpu_${_precision}"
                             --device gpu --gpu-precision ${_precision} --paths-csv)
            set_tests_properties(gpu_run_${_scene}_${_precision} PROPERTIES
                                 SKIP_RETURN_CODE 77 TIMEOUT 600 FIXTURES_SETUP sbr_gpu_${_scene}_${_precision})
        endforeach()
    endforeach()
    foreach(_scene office bahen)
        add_test(NAME gpu_run_${_scene}_fine_double
                 COMMAND sbrsolver "${SBR_SRC_DIR}/data/fine/${_scene}.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_fine_gpu_double"
                         --device gpu --gpu-precision double --paths-csv)
        set_tests_properties(gpu_run_${_scene}_fine_double PROPERTIES
                             SKIP_RETURN_CODE 77 TIMEOUT 1800 FIXTURES_SETUP sbr_gpu_${_scene}_fine_double)
    endforeach()
    # The plan scene with its two access points, 74,181 segments and 3.8
    # million candidates: about 11 s on lavapipe with four cores.
    foreach(_precision double single)
        add_test(NAME gpu_run_large_office_${_precision}
                 COMMAND sbrsolver "${SBR_SRC_DIR}/data/large_office/large_office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_large_office_gpu_${_precision}"
                         --device gpu --gpu-precision ${_precision} --paths-csv)
        set_tests_properties(gpu_run_large_office_${_precision} PROPERTIES
                             SKIP_RETURN_CODE 77 TIMEOUT 1800 FIXTURES_SETUP sbr_gpu_large_office_${_precision})
    endforeach()
    # A second single-precision office run, to compare with the first.
    add_test(NAME gpu_run_office_single_repeat
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_gpu_single_repeat"
                     --device gpu --gpu-precision single --paths-csv)
    set_tests_properties(gpu_run_office_single_repeat PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 600 FIXTURES_SETUP sbr_gpu_office_single_repeat)
    # The office in small batches with a small candidate buffer: 642
    # launches of at most 2 segments of 96 bytes at 40,000 bytes per batch
    # give four batches, and 500 pairs for 1,687 candidates force reruns.
    add_test(NAME gpu_run_office_batches
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_gpu_batches"
                     --device gpu --gpu-precision double --paths-csv)
    set_tests_properties(gpu_run_office_batches PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 600 FIXTURES_SETUP sbr_gpu_office_batches
                         ENVIRONMENT "SBR_GPU_BATCH_BYTES=40000;SBR_GPU_CANDIDATE_CAPACITY=500")
endif()

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
    # Two transmitters: the overlays are built from a scene with two
    # transmitters and paths tagged with them.
    add_test(NAME viewer_smoke_mock_two_tx COMMAND sbr_viewer_smoke --backend openGL_mock --frames 2 --transmitters 2)
    set_tests_properties(viewer_smoke_mock viewer_smoke_mock_volume viewer_smoke_mock_two_tx PROPERTIES TIMEOUT 300)
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
    add_test(NAME viewer_cli_rejects_backend
             COMMAND sbrsolver "${SBR_SRC_DIR}/data/office.ini" "${CMAKE_CURRENT_BINARY_DIR}/test_output_viewer"
                     --viewer-backend openGL_bogus --viewer-frames 1)
    set_tests_properties(viewer_cli_rejects_frames viewer_cli_rejects_screenshot_type viewer_cli_rejects_backend
                         PROPERTIES WILL_FAIL TRUE TIMEOUT 60)
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
    # tools/compare_paths.py matches the paths of two runs by receiver and
    # key. The office run against itself passes; against the two-transmitter
    # run it must fail, because the second transmitter's paths exist on one
    # side only. The failing run is the proof that the tool sees a difference.
    add_test(NAME paths_tool_pass
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office")
    set_tests_properties(paths_tool_pass PROPERTIES FIXTURES_REQUIRED sbr_run_office)
    add_test(NAME paths_tool_fail
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx")
    set_tests_properties(paths_tool_fail PROPERTIES FIXTURES_REQUIRED "sbr_run_office;sbr_run_office_two_tx" WILL_FAIL TRUE)
    # The same for tools/compare_segments.py: the two-transmitter run has
    # more launches, so the comparison must fail.
    add_test(NAME segments_tool_pass
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_segments.py"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office/Segments.csv" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office/Segments.csv")
    set_tests_properties(segments_tool_pass PROPERTIES FIXTURES_REQUIRED sbr_run_office)
    add_test(NAME segments_tool_fail
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_segments.py"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office/Segments.csv" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx/Segments.csv")
    set_tests_properties(segments_tool_fail PROPERTIES FIXTURES_REQUIRED "sbr_run_office;sbr_run_office_two_tx" WILL_FAIL TRUE)
    if(SBR_BUILD_GPU)
        # The GPU runs against the CPU runs. In double precision the paths
        # must be the same at every receiver; in single precision a bounded
        # number of receivers may differ. The single bounds are twice the
        # counts measured on lavapipe (Mesa 25.2.8), at least 5.
        set(_sbr_single_bound_ground 5)
        set(_sbr_single_bound_corner 5)
        set(_sbr_single_bound_office 5)
        set(_sbr_single_bound_bahen 5)
        set(_sbr_single_bound_office_two_tx 5)
        foreach(_scene IN LISTS _sbr_gpu_scenes)
            add_test(NAME gpu_paths_${_scene}_double
                     COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                             "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_gpu_double")
            set_tests_properties(gpu_paths_${_scene}_double PROPERTIES SKIP_RETURN_CODE 77
                                 FIXTURES_REQUIRED "sbr_gpu_${_scene}_double;sbr_run_${_scene}")
            add_test(NAME gpu_paths_${_scene}_single
                     COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                             "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_gpu_single"
                             --max-differing ${_sbr_single_bound_${_scene}})
            set_tests_properties(gpu_paths_${_scene}_single PROPERTIES SKIP_RETURN_CODE 77
                                 FIXTURES_REQUIRED "sbr_gpu_${_scene}_single;sbr_run_${_scene}")
        endforeach()
        foreach(_scene office bahen)
            add_test(NAME gpu_paths_${_scene}_fine_double
                     COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                             "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_fine" "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_fine_gpu_double")
            set_tests_properties(gpu_paths_${_scene}_fine_double PROPERTIES SKIP_RETURN_CODE 77
                                 FIXTURES_REQUIRED "sbr_gpu_${_scene}_fine_double;sbr_run_${_scene}_fine")
        endforeach()
        add_test(NAME gpu_paths_large_office_double
                 COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                         "${CMAKE_CURRENT_BINARY_DIR}/test_output_large_office" "${CMAKE_CURRENT_BINARY_DIR}/test_output_large_office_gpu_double")
        set_tests_properties(gpu_paths_large_office_double PROPERTIES SKIP_RETURN_CODE 77
                             FIXTURES_REQUIRED "sbr_gpu_large_office_double;sbr_run_large_office")
        add_test(NAME gpu_paths_large_office_single
                 COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                         "${CMAKE_CURRENT_BINARY_DIR}/test_output_large_office" "${CMAKE_CURRENT_BINARY_DIR}/test_output_large_office_gpu_single"
                         --max-differing 5)
        set_tests_properties(gpu_paths_large_office_single PROPERTIES SKIP_RETURN_CODE 77
                             FIXTURES_REQUIRED "sbr_gpu_large_office_single;sbr_run_large_office")
        # Two single-precision runs give the same file, and the batched run
        # gives the same paths and fields as the run in one batch.
        add_test(NAME gpu_repeat_office
                 COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                         "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_gpu_single" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_gpu_single_repeat"
                         --field-rtol 0)
        set_tests_properties(gpu_repeat_office PROPERTIES SKIP_RETURN_CODE 77
                             FIXTURES_REQUIRED "sbr_gpu_office_single;sbr_gpu_office_single_repeat")
        add_test(NAME gpu_batches_office
                 COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_paths.py"
                         "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_gpu_double" "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_gpu_batches"
                         --field-rtol 0)
        set_tests_properties(gpu_batches_office PROPERTIES SKIP_RETURN_CODE 77
                             FIXTURES_REQUIRED "sbr_gpu_office_double;sbr_gpu_office_batches")
    endif()
    # Dense layouts: one finite field per receiver. Every receiver of the
    # ground line and the corner floor sees the transmitter. In the office
    # the room closed by the glass panel and the rooms without a line of
    # sight (658 receivers) and in Bahen a few door recesses (7) get no field
    # with one reflection and no transmission, so those allow zeros.
    set(_sbr_fine_rx_ground RX_Ground_fine.dat)
    set(_sbr_fine_rx_corner RX_Corner_fine.dat)
    set(_sbr_fine_rx_office RX_Office_fine.dat)
    set(_sbr_fine_rx_bahen RX_Bahen_fine.dat)
    set(_sbr_fine_check_ground --all-nonzero)
    set(_sbr_fine_check_corner --all-nonzero)
    set(_sbr_fine_check_office --max-zero 700)
    set(_sbr_fine_check_bahen --max-zero 10)
    # The plan scene: one finite field per receiver. The rooms that neither
    # access point reaches within six wall transmissions get no field:
    # 1003 receivers.
    add_test(NAME sbrsolver_fields_large_office
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_fields.py"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_large_office/ElectricField.csv"
                     "${SBR_SRC_DIR}/data/large_office/RX_large_office.dat" --max-zero 1100)
    set_tests_properties(sbrsolver_fields_large_office PROPERTIES FIXTURES_REQUIRED sbr_run_large_office)
    # Two transmitters radiating at once: the field of the full run is the
    # sum of the two single-transmitter runs (which the check runs itself)
    # and every path carries its transmitter. The second transmitter is at
    # the other end of the room, so the total field is zero only where
    # neither one reaches.
    add_test(NAME sbrsolver_two_tx_superposition
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_superposition.py"
                     "$<TARGET_FILE:sbrsolver>" "${SBR_SRC_DIR}/data/office_two_tx.ini"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx_single" --timeout 250)
    add_test(NAME sbrsolver_fields_office_two_tx
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_fields.py"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx/ElectricField.csv"
                     "${SBR_SRC_DIR}/data/RX_Office.dat" --max-zero 30)
    set_tests_properties(sbrsolver_two_tx_superposition sbrsolver_fields_office_two_tx
                         PROPERTIES FIXTURES_REQUIRED sbr_run_office_two_tx TIMEOUT 600)
    # The same with the powers of the two transmitters added: the field of
    # the full run is the root of the summed squares of the single fields.
    add_test(NAME sbrsolver_two_tx_power_sum
             COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_superposition.py"
                     "$<TARGET_FILE:sbrsolver>" "${SBR_SRC_DIR}/data/office_two_tx_power.ini"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx_power"
                     "${CMAKE_CURRENT_BINARY_DIR}/test_output_office_two_tx_power_single" --timeout 250)
    set_tests_properties(sbrsolver_two_tx_power_sum PROPERTIES FIXTURES_REQUIRED sbr_run_office_two_tx_power TIMEOUT 600)
    foreach(_scene IN LISTS _sbr_scenes)
        add_test(NAME sbrsolver_fields_${_scene}_fine
                 COMMAND "${Python3_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_fields.py"
                         "${CMAKE_CURRENT_BINARY_DIR}/test_output_${_scene}_fine/ElectricField.csv"
                         "${SBR_SRC_DIR}/data/fine/${_sbr_fine_rx_${_scene}}" ${_sbr_fine_check_${_scene}})
        set_tests_properties(sbrsolver_fields_${_scene}_fine PROPERTIES FIXTURES_REQUIRED sbr_run_${_scene}_fine)
    endforeach()

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
