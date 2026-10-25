cmake_minimum_required(VERSION 3.14)

# The plain .cpp files are compiled by the host C++ compiler (g++, clang++);
# only src/engine.cpp goes to nvcc (CMAKE_CUDA_COMPILER). If CXX (environment
# or -DCMAKE_CXX_COMPILER) points to nvcc, as in some Kokkos setups, use that
# nvcc as the CUDA compiler and the host compiler g++ for C++.
if(DEFINED CMAKE_CXX_COMPILER)
  set(GEMPIC_REQUESTED_CXX "${CMAKE_CXX_COMPILER}")
else()
  set(GEMPIC_REQUESTED_CXX "$ENV{CXX}")
endif()
if(GEMPIC_REQUESTED_CXX)
  get_filename_component(GEMPIC_REQUESTED_CXX_NAME "${GEMPIC_REQUESTED_CXX}" NAME)
  if(GEMPIC_REQUESTED_CXX_NAME STREQUAL "nvcc")
    find_program(GEMPIC_HOST_CXX NAMES g++ c++ clang++)
    if(NOT GEMPIC_HOST_CXX)
      message(FATAL_ERROR "CXX is nvcc, and no host C++ compiler (g++) was found: set CXX=g++")
    endif()
    message(STATUS "CXX is nvcc: C++ files -> ${GEMPIC_HOST_CXX}, CUDA -> ${GEMPIC_REQUESTED_CXX}")
    set(CMAKE_CXX_COMPILER "${GEMPIC_HOST_CXX}" CACHE FILEPATH "C++ compiler" FORCE)
    if(NOT DEFINED CMAKE_CUDA_COMPILER AND NOT DEFINED ENV{CUDACXX})
      find_program(GEMPIC_NVCC_FULL "${GEMPIC_REQUESTED_CXX}")
      set(CMAKE_CUDA_COMPILER "${GEMPIC_NVCC_FULL}" CACHE FILEPATH "CUDA compiler")
    endif()
  endif()
endif()

project(gempic_rz VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # -std=c++17, not gnu++17 (also accepted by nvcc)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
option(GEMPIC_NATIVE "optimize for the build machine (-march=native)" ON)
option(GEMPIC_OPENMP "CPU build: parallel engine with OpenMP threads (backend = openmp)" ON)
option(GEMPIC_CUDA "build the GPU engine for NVIDIA GPUs (backend = cuda)" OFF)
option(GEMPIC_HDF5 "HDF5 output (output.format = hdf5), if the HDF5 library is found" ON)
option(GEMPIC_MPI "MPI: domain decomposition along zeta (pipeline), if MPI is found" ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wall> $<$<COMPILE_LANGUAGE:CXX>:-Wextra>)
  if(GEMPIC_NATIVE)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-march=native>)
  endif()
endif()

if(GEMPIC_CUDA)
  # GPU architectures: e.g. -DCMAKE_CUDA_ARCHITECTURES=native, or 80 (A100), 90 (H100)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89 90)
  endif()
  enable_language(CUDA)
  # the CUDA runtime is linked by CMake itself (statically: no LD_LIBRARY_PATH
  # needed at run time); find_package(CUDAToolkit) is deliberately not used, it
  # needs FindThreads, which fails on some systems
  if(NOT DEFINED CMAKE_CUDA_RUNTIME_LIBRARY)
    set(CMAKE_CUDA_RUNTIME_LIBRARY Static)
  endif()
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set_source_files_properties(src/engine.cpp PROPERTIES LANGUAGE CUDA)
endif()

# library with everything except main()
add_library(gemcore STATIC
  src/units.cpp src/deck.cpp src/config.cpp src/atomic_data.cpp src/ionization.cpp src/mapping.cpp
  src/grid.cpp src/solver.cpp src/smoothing.cpp src/simulation.cpp src/diagnostics.cpp src/output_hdf5.cpp src/engine.cpp src/comm.cpp src/pipeline.cpp src/eigen.cpp src/expr.cpp src/beam_import.cpp src/laser_file.cpp)
target_include_directories(gemcore PUBLIC src)

# HDF5 output (C library)
set(GEMPIC_HAVE_HDF5 OFF)
if(GEMPIC_HDF5)
  find_package(HDF5 COMPONENTS C)
  if(HDF5_FOUND)
    set(GEMPIC_HAVE_HDF5 ON)
    message(STATUS "HDF5 ${HDF5_VERSION}: output.format = hdf5 available")
  else()
    message(WARNING "HDF5 not found: output.format = hdf5 is disabled (install e.g. libhdf5-dev)")
  endif()
endif()
function(gempic_use_hdf5 target)
  if(GEMPIC_HAVE_HDF5)
    target_compile_definitions(${target} PUBLIC GEMPIC_HAVE_HDF5 ${HDF5_DEFINITIONS})
    target_include_directories(${target} PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(${target} PUBLIC ${HDF5_LIBRARIES})
  endif()
endfunction()
gempic_use_hdf5(gemcore)

# MPI (domain decomposition along zeta)
set(GEMPIC_HAVE_MPI OFF)
if(GEMPIC_MPI)
  find_package(MPI COMPONENTS CXX)
  if(MPI_CXX_FOUND)
    set(GEMPIC_HAVE_MPI ON)
    message(STATUS "MPI ${MPI_CXX_VERSION}: run with  mpirun -np <P> gempic_rz <deck>")
  else()
    message(STATUS "MPI not found: single-process build")
  endif()
endif()
function(gempic_use_mpi target)
  if(GEMPIC_HAVE_MPI)
    target_compile_definitions(${target} PUBLIC GEMPIC_HAVE_MPI)
    target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
  endif()
endfunction()
gempic_use_mpi(gemcore)
if(GEMPIC_CUDA)
  target_compile_definitions(gemcore PUBLIC GEMPIC_HAVE_CUDA)
  target_compile_options(gemcore PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda --expt-relaxed-constexpr>)
elseif(GEMPIC_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(gemcore PUBLIC OpenMP::OpenMP_CXX)
    # test build: the GPU code paths emulated on CPU threads (see portable.hpp)
    add_executable(gempic_rz_gpuemu src/main.cpp src/units.cpp src/deck.cpp src/config.cpp src/atomic_data.cpp
      src/ionization.cpp src/mapping.cpp src/grid.cpp src/solver.cpp src/smoothing.cpp src/simulation.cpp src/diagnostics.cpp
      src/output_hdf5.cpp src/engine.cpp src/comm.cpp src/pipeline.cpp src/eigen.cpp src/expr.cpp src/beam_import.cpp src/laser_file.cpp)
    gempic_use_hdf5(gempic_rz_gpuemu)
    gempic_use_mpi(gempic_rz_gpuemu)
    target_include_directories(gempic_rz_gpuemu PRIVATE src)
    target_compile_definitions(gempic_rz_gpuemu PRIVATE GEM_EMULATE_GPU)
    target_link_libraries(gempic_rz_gpuemu PRIVATE OpenMP::OpenMP_CXX)
  else()
    message(WARNING "OpenMP not found: backend = openmp runs on one thread")
  endif()
endif()

add_executable(gempic_rz src/main.cpp)
target_link_libraries(gempic_rz gemcore)

# tests
enable_testing()
add_executable(test_deck tests/test_deck.cpp)
target_link_libraries(test_deck gemcore)
target_compile_definitions(test_deck PRIVATE GEM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_test(NAME deck_parser COMMAND test_deck)

add_executable(compare_columns tools/compare_columns.cpp)
add_executable(check_physics tests/check_physics.cpp)
target_link_libraries(check_physics gemcore)
# physics and regression tests for every backend of this build
set(GEMPIC_TEST_BACKENDS serial)
if(GEMPIC_CUDA)
  list(APPEND GEMPIC_TEST_BACKENDS cuda)       # needs a GPU when ctest runs
elseif(GEMPIC_OPENMP)
  list(APPEND GEMPIC_TEST_BACKENDS openmp)
endif()
foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
  if(b STREQUAL "serial")
    set(sfx "")
  else()
    set(sfx "_${b}")
  endif()
  add_test(NAME physics${sfx}
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:gempic_rz> -DCHK=$<TARGET_FILE:check_physics> -DBACKEND=${b}
            -DSRC=${CMAKE_SOURCE_DIR} -P ${CMAKE_SOURCE_DIR}/tests/physics.cmake)
  add_test(NAME regression${sfx}
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:gempic_rz> -DCMP=$<TARGET_FILE:compare_columns> -DBACKEND=${b}
            -DSRC=${CMAKE_SOURCE_DIR} -P ${CMAKE_SOURCE_DIR}/tests/regression.cmake)
endforeach()
if(TARGET gempic_rz_gpuemu)
  # the GPU code paths, emulated on CPU threads, against the serial reference
  add_test(NAME regression_gpu_emulation
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:gempic_rz_gpuemu> -DCMP=$<TARGET_FILE:compare_columns>
            -DBACKEND=openmp -DSRC=${CMAKE_SOURCE_DIR} -P ${CMAKE_SOURCE_DIR}/tests/regression.cmake)
endif()

# MPI pipeline (2-4 ranks; MPIEXEC_PREFLAGS e.g. for a batch system)
if(GEMPIC_HAVE_MPI AND GEMPIC_HAVE_HDF5)
  add_test(NAME mpi_pipeline
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/mpi.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:compare_columns>
            ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS})
endif()

# HDF5 output and checkpoint / restart
if(GEMPIC_HAVE_HDF5)
  add_executable(h5_extract tools/h5_extract.cpp)
  gempic_use_hdf5(h5_extract)
  foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
    add_test(NAME restart_${b}
      COMMAND sh ${CMAKE_SOURCE_DIR}/tests/restart.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:compare_columns> ${b})
  endforeach()
  add_test(NAME hdf5_output
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/hdf5.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:h5_extract>
            $<TARGET_FILE:compare_columns>)
  # free-hand zeta and z grids: backend against serial, restart
  foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
    add_test(NAME flexgrid_${b}
      COMMAND sh ${CMAKE_SOURCE_DIR}/tests/flexgrid.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:compare_columns> ${b})
  endforeach()
  if(TARGET gempic_rz_gpuemu)
    add_test(NAME flexgrid_gpu_emulation
      COMMAND sh ${CMAKE_SOURCE_DIR}/tests/flexgrid.sh $<TARGET_FILE:gempic_rz_gpuemu> $<TARGET_FILE:compare_columns>
              openmp)
  endif()
  # profiles given as functions (function parser) against the built-in profiles; undulator taper
  foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
    add_test(NAME functions_${b}
      COMMAND sh ${CMAKE_SOURCE_DIR}/tests/functions.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:compare_columns> ${b})
  endforeach()
  # ion blocks with their own density profile (dopants)
  foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
    add_test(NAME dopant_${b}
      COMMAND sh ${CMAKE_SOURCE_DIR}/tests/dopant.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:compare_columns> ${b})
  endforeach()
  # radial smoothing of the plasma sources (plasma.smooth_length)
  foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
    add_test(NAME smoothing_${b}
      COMMAND sh ${CMAKE_SOURCE_DIR}/tests/smoothing.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:compare_columns> ${b})
  endforeach()
  # transport elements: quadrupole, plasma lens, solenoid, chicane against linear optics
  foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
    add_test(NAME transport_${b}
      COMMAND sh ${CMAKE_SOURCE_DIR}/tests/transport.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:h5_extract> ${b})
  endforeach()
  # beams from files (openPMD, native), selection, resampling, Fawley shot noise
  add_test(NAME beam_import
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/beam_import.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:compare_columns>)
  # openPMD output: validator (if installed), fields and particles against the native files
  add_test(NAME openpmd
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/openpmd.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:h5_extract>)
  set_tests_properties(openpmd PROPERTIES SKIP_RETURN_CODE 77)
  add_test(NAME undulator
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/undulator.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:h5_extract>)
  # beam push in a planar undulator at dz = lambda_u / 16: numerics.beam_field = mean against level
  foreach(b IN LISTS GEMPIC_TEST_BACKENDS)
    add_test(NAME planar_step_${b} COMMAND sh ${CMAKE_SOURCE_DIR}/tests/planar_step.sh $<TARGET_FILE:gempic_rz> ${b})
  endforeach()
  # the parallel backend of this build (cuda in a CUDA build, which rejects openmp; else openmp or serial)
  list(GET GEMPIC_TEST_BACKENDS -1 GEMPIC_PAR_BACKEND)
  # bunch.axis_refine: more macroparticles near the axis at the same charge density
  add_test(NAME axis_refine
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/axis_refine.sh $<TARGET_FILE:gempic_rz> ${GEMPIC_PAR_BACKEND})
  # bunch.distribution_r = function, ux_function, uy_function: a helical ring on its matched orbits
  add_test(NAME bunch_functions
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/bunch_functions.sh $<TARGET_FILE:gempic_rz> ${GEMPIC_PAR_BACKEND})
  # flying focus: luminal Bessel pulse fed through the radial absorber layer (pattern kept, moves at c)
  add_test(NAME flying_focus
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/flying_focus.sh $<TARGET_FILE:gempic_rz>)
  # lasers from and to LASY envelope files: import, exact propagation, boundary injection, export
  add_test(NAME lasy COMMAND sh ${CMAKE_SOURCE_DIR}/tests/lasy.sh $<TARGET_FILE:gempic_rz>)
  set_tests_properties(lasy PROPERTIES SKIP_RETURN_CODE 77)
  # energy spread dE/E of bunches.txt against the particle data (gamma = 1000, dE/E down to 1e-6)
  add_test(NAME energy_spread COMMAND sh ${CMAKE_SOURCE_DIR}/tests/energy_spread.sh $<TARGET_FILE:gempic_rz>)
  set_tests_properties(energy_spread PROPERTIES SKIP_RETURN_CODE 77)
  # GUI (tools/gempic_gui.py): runs started through its HTTP API; fields against the code's gather
  add_test(NAME gui COMMAND sh ${CMAKE_SOURCE_DIR}/tests/gui.sh $<TARGET_FILE:gempic_rz>)
  set_tests_properties(gui PROPERTIES SKIP_RETURN_CODE 77)
  add_test(NAME cherenkov
    COMMAND sh ${CMAKE_SOURCE_DIR}/tests/cherenkov.sh $<TARGET_FILE:gempic_rz> $<TARGET_FILE:h5_extract>)
endif()
