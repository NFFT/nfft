cmake_minimum_required(VERSION 3.20)

project(NFFT3
  VERSION 3.5.4
  DESCRIPTION "Nonequispaced Fast Fourier Transforms"
  LANGUAGES C)

set(NFFT_VERSION_TYPE "alpha" CACHE STRING "Version type suffix (e.g. alpha)")

list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
option(BUILD_SHARED_LIBS "Build shared libraries (DLLs on Windows)" ON)

# Options
option(NFFT_ENABLE_FLOAT       "Build single (float) precision"   OFF)
option(NFFT_ENABLE_LONG_DOUBLE "Build long double precision"      OFF)
option(NFFT_ENABLE_OPENMP      "Build the OpenMP (_omp) library"  OFF)
option(NFFT_ENABLE_MAXOPT      "Use aggressive compiler flags"    ON)
option(NFFT_ENABLE_TESTS       "Build the CUnit test suite"       OFF)

# Benchmarks: a single knob that both selects the CodSpeed measurement mode and
# (when not "off") enables the benchmark build. The mode is baked in at build time
# and forwarded to codspeed-cpp's CODSPEED_MODE.
#   off        - benchmarks not built (default).
#   simulation - deterministic instruction count via callgrind (needs valgrind to
#                measure); the mode CI gates on.
#   walltime   - wall-clock timing; the binary writes a local JSON of per-benchmark
#                stats to $CODSPEED_PROFILE_FOLDER (offline, no runner/upload).
set(NFFT_BENCHMARK_MODE "off" CACHE STRING
    "Benchmark build + CodSpeed mode: off|simulation|walltime")
set_property(CACHE NFFT_BENCHMARK_MODE PROPERTY STRINGS off simulation walltime)

# If tests are enabled, allow opt-in to set of more exhaustive tests.
include(CMakeDependentOption)
cmake_dependent_option(NFFT_ENABLE_EXHAUSTIVE_UNIT_TESTS
  "Enable larger/slower exhaustive unit tests (requires NFFT_ENABLE_TESTS=ON)"
  OFF "NFFT_ENABLE_TESTS" OFF)

option(NFFT_ENABLE_EXAMPLES      "Build the example programs"     ON)
option(NFFT_ENABLE_APPLICATIONS  "Build the application programs" ON)

set(NFFT_WINDOW "kaiserbessel" CACHE STRING
    "Window function: kaiserbessel|gaussian|bspline|sinc|dirac")

# Language interfaces
option(NFFT_ENABLE_JULIA "Build Julia interface (double + shared only)" OFF)
option(NFFT_WITH_OCTAVE  "Build Octave mex interface"             OFF)
set(NFFT_WITH_MATLAB "" CACHE PATH
    "Path to a MATLAB installation root; non-empty builds the mex interface against MATLAB")

# If building a mex interface, allow opt-in to building with OpenMP kernel. 
cmake_dependent_option(NFFT_ENABLE_MATLAB_THREADS
  "Build mex interface with the threaded (OpenMP) kernel"
  "${NFFT_ENABLE_OPENMP}" "NFFT_WITH_OCTAVE OR NFFT_WITH_MATLAB" OFF)

# Guards

# Precision mutually exclusivity.
if(NFFT_ENABLE_FLOAT AND NFFT_ENABLE_LONG_DOUBLE)
  message(FATAL_ERROR "NFFT_ENABLE_FLOAT and NFFT_ENABLE_LONG_DOUBLE are mutually exclusive")
endif()

# Julia interface only supports double.
if(NFFT_ENABLE_JULIA AND (NFFT_ENABLE_FLOAT OR NFFT_ENABLE_LONG_DOUBLE))
  message(FATAL_ERROR "NFFT_ENABLE_JULIA requires double precision.")
endif()

# Julia interface requires shared libraries.
if(NFFT_ENABLE_JULIA AND NOT BUILD_SHARED_LIBS)
  message(FATAL_ERROR "NFFT_ENABLE_JULIA requires BUILD_SHARED_LIBS=ON.")
endif()

# Julia wrappers resolve library symbols by name at runtime, so each Julia module 
# must link the interface-kernel archive whole, supported since CMake 3.24.
if(NFFT_ENABLE_JULIA AND CMAKE_VERSION VERSION_LESS 3.24)
  message(FATAL_ERROR "NFFT_ENABLE_JULIA requires CMake >= 3.24 "
                      "(whole-archive link of the interface kernel).")
endif()

# Octave mex only supports double or single.
if(NFFT_WITH_OCTAVE AND NFFT_ENABLE_LONG_DOUBLE)
  message(FATAL_ERROR "The Octave mex interface does not support long-double precision.")
endif()

# MATLAB mex only supports double or single.
if(NFFT_WITH_MATLAB AND NFFT_ENABLE_LONG_DOUBLE)
  message(FATAL_ERROR "The MATLAB mex interface does not support long-double precision.")
endif()

# Matlab/Octave mex interfaces are mutually exclusive.
if(NFFT_WITH_OCTAVE AND NFFT_WITH_MATLAB)
  message(FATAL_ERROR "NFFT_WITH_OCTAVE and NFFT_WITH_MATLAB are mutually exclusive "
                      "(both build the same <mod>mex modules in matlab/).")
endif()

# OpenMP required for threaded Matlab interface.
if(NFFT_ENABLE_MATLAB_THREADS AND NOT NFFT_ENABLE_OPENMP)
  message(FATAL_ERROR "NFFT_ENABLE_MATLAB_THREADS requires NFFT_ENABLE_OPENMP=ON.")
endif()

# Precision fixed per build tree.
if(NFFT_ENABLE_FLOAT)
  set(NFFT_PREC_SUFFIX "f")
  set(NFFT_PRECISION_MACRO "NFFT_PRECISION_SINGLE")
elseif(NFFT_ENABLE_LONG_DOUBLE)
  set(NFFT_PREC_SUFFIX "l")
  set(NFFT_PRECISION_MACRO "NFFT_PRECISION_LONG_DOUBLE")
else()
  set(NFFT_PREC_SUFFIX "")
  set(NFFT_PRECISION_MACRO "NFFT_PRECISION_DOUBLE")
endif()

message(STATUS "NFFT3 ${PROJECT_VERSION}${NFFT_VERSION_TYPE}: precision suffix='${NFFT_PREC_SUFFIX}', window=${NFFT_WINDOW}")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

include(nfft_config)

# Generate config.h into build dir.
configure_file(${PROJECT_SOURCE_DIR}/cmake.config.h.in
               ${PROJECT_BINARY_DIR}/config.h @ONLY)

# Generate ticks.h from existing template. TICKS_PER_SECOND only affects
# MEASURE_TIME timing (off by default); fallback is 1.
set(TICKS_PER_SECOND 1)
configure_file(${PROJECT_SOURCE_DIR}/include/ticks.h.in
               ${PROJECT_BINARY_DIR}/ticks.h @ONLY)

# Include search order: build dir (generated headers) BEFORE source include/.
set(NFFT_GENERATED_INCLUDE_DIR ${PROJECT_BINARY_DIR})

# FFTW3 required dependency.
find_package(FFTW3 REQUIRED)

# Define HAVE_FFTW_THREADS to true if an _omp (or _threads) FFTW3 lib exists.
if(NFFT_ENABLE_OPENMP AND (FFTW3_OMP_FOUND OR FFTW3_THREADS_FOUND))
  set(HAVE_FFTW_THREADS ON)
  # Regenerate config.h now that HAVE_FFTW_THREADS is known.
  configure_file(${PROJECT_SOURCE_DIR}/cmake.config.h.in
                 ${PROJECT_BINARY_DIR}/config.h @ONLY)
endif()

if(NFFT_ENABLE_OPENMP AND NOT FFTW3_OMP_FOUND)
  if(FFTW3_THREADS_FOUND)
    message(WARNING "Only generic multi-threaded FFTW3 library found. This may be less "
      "efficient because NFFT will use OpenMP threads while FFTW3 uses its own threading mechanism.")
  else()
    message(WARNING "Only single-threaded FFTW3 library found. Multi-threading will only "
      "be used in NFFT code, not in FFTW3 computations.")
  endif()
endif()

if(NFFT_ENABLE_OPENMP)
  # On macOS, two libomp.dylib copies at different paths can both load into one
  # process and abort with "OMP: Error #15" (e.g. Homebrew's llvm ships its own,
  # Homebrew's fftw uses the libomp formula). Preset FindOpenMP's library to the
  # libomp that FFTW3's OpenMP library uses. OpenMP_omp_LIBRARY is used for
  # clang, OpenMP_libomp_LIBRARY for AppleClang.
  if(APPLE AND FFTW3_OMP_FOUND AND NOT OpenMP_omp_LIBRARY AND NOT OpenMP_libomp_LIBRARY)
    find_program(NFFT_OTOOL otool)
    mark_as_advanced(NFFT_OTOOL)
    if(NFFT_OTOOL)
      execute_process(COMMAND ${NFFT_OTOOL} -L ${FFTW3_OMP_LIBRARY}
                      OUTPUT_VARIABLE _nfft_fftw_deps ERROR_QUIET)
      string(REGEX MATCH "[^ \t\n]*/libomp\\.dylib" _nfft_fftw_omp "${_nfft_fftw_deps}")
      if(_nfft_fftw_omp AND EXISTS "${_nfft_fftw_omp}")
        set(OpenMP_omp_LIBRARY "${_nfft_fftw_omp}" CACHE FILEPATH "OpenMP runtime used by FFTW3")
        set(OpenMP_libomp_LIBRARY "${_nfft_fftw_omp}" CACHE FILEPATH "OpenMP runtime used by FFTW3")
        message(STATUS "Using the OpenMP runtime of FFTW3: ${_nfft_fftw_omp}")
      endif()
      unset(_nfft_fftw_deps)
      unset(_nfft_fftw_omp)
    endif()
  endif()

  find_package(OpenMP COMPONENTS C REQUIRED)
  include(CheckCSourceCompiles)

  # Link-test `#pragma omp atomic` on the SELECTED precision's real type, and
  # discover whether -latomic is needed. The test type is chosen by the same 
  # NFFT_SINGLE/NFFT_LDOUBLE defines the library compiles with.
  set(_nfft_atomic_test "
    #include <omp.h>
    #include <stdint.h>
    int main(void){
    #if defined(NFFT_SINGLE)
      float r = 0.0f, increment = 2.0f, expected_min = 150.0f;
    #elif defined(NFFT_LDOUBLE)
      long double r = 0.0L, increment = 2.0L, expected_min = 150.0L;
    #else
      double r = 0.0, increment = 2.0, expected_min = 150.0;
    #endif
      int i;
      volatile int loop_count = (int)((uintptr_t)&r % 50) + 75;
      #pragma omp parallel for
      for (i = 0; i < loop_count; i++) {
        #pragma omp atomic
        r += increment;
      }
      return (r < expected_min) ? 1 : 0;
    }")

  if(NFFT_ENABLE_FLOAT)
    set(CMAKE_REQUIRED_DEFINITIONS "-DNFFT_SINGLE=1")
  elseif(NFFT_ENABLE_LONG_DOUBLE)
    set(CMAKE_REQUIRED_DEFINITIONS "-DNFFT_LDOUBLE=1")
  else()
    set(CMAKE_REQUIRED_DEFINITIONS "")
  endif()

  # Try with no extra lib first, then with -latomic.
  set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_C)
  unset(NFFT_OPENMP_ATOMIC_NONE CACHE)
  check_c_source_compiles("${_nfft_atomic_test}" NFFT_OPENMP_ATOMIC_NONE)
  if(NFFT_OPENMP_ATOMIC_NONE)
    set(NFFT_OPENMP_ATOMIC_LIBS "" CACHE INTERNAL "Extra libs for OpenMP atomic on the real type")
  else()
    set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_C atomic)
    unset(NFFT_OPENMP_ATOMIC_LATOMIC CACHE)
    check_c_source_compiles("${_nfft_atomic_test}" NFFT_OPENMP_ATOMIC_LATOMIC)
    if(NFFT_OPENMP_ATOMIC_LATOMIC)
      set(NFFT_OPENMP_ATOMIC_LIBS "atomic" CACHE INTERNAL "Extra libs for OpenMP atomic on the real type")
    else()
      message(FATAL_ERROR
        "OpenMP atomic on the real type is not supported (tried with and without -latomic)")
    endif()
  endif()

  unset(CMAKE_REQUIRED_DEFINITIONS)
  unset(CMAKE_REQUIRED_LIBRARIES)
  unset(_nfft_atomic_test)
  message(STATUS "OpenMP atomic on real type: extra libs = '${NFFT_OPENMP_ATOMIC_LIBS}'")

  # Start NFFT's and FFTW3's OpenMP code in one process. Fails if two OpenMP
  # runtimes get loaded.
  if(FFTW3_OMP_FOUND AND NOT CMAKE_CROSSCOMPILING)
    include(CheckCSourceRuns)
    set(_f "fftw${NFFT_PREC_SUFFIX}")
    set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_C FFTW3::fftw3_omp FFTW3::fftw3)
    unset(NFFT_OPENMP_FFTW3_RUNS CACHE)
    check_c_source_runs("
      #include <omp.h>
      #include <fftw3.h>
      int main(void){
        int n = 0;
        ${_f}_complex *x;
        ${_f}_plan p;
        #pragma omp parallel reduction(+:n)
        n += 1;
        if (n < 1 || !${_f}_init_threads())
          return 1;
        ${_f}_plan_with_nthreads(2);
        x = ${_f}_alloc_complex(64);
        p = ${_f}_plan_dft_1d(64, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
        ${_f}_execute(p);
        ${_f}_destroy_plan(p);
        ${_f}_free(x);
        ${_f}_cleanup_threads();
        return 0;
      }" NFFT_OPENMP_FFTW3_RUNS)
    unset(CMAKE_REQUIRED_LIBRARIES)
    unset(_f)
    if(NOT NFFT_OPENMP_FFTW3_RUNS)
      message(FATAL_ERROR
        "A program that uses both OpenMP and the OpenMP FFTW3 library failed to run. "
        "This is most likely because two different OpenMP runtimes were loaded, e.g. "
        "libomp from Homebrew's llvm and from Homebrew's libomp formula. Set "
        "OpenMP_omp_LIBRARY (clang) or OpenMP_libomp_LIBRARY (AppleClang) to the "
        "OpenMP runtime that FFTW3 uses.")
    endif()
  endif()
endif()

include(nfft_maxopt)   # Aggressive compiler options.
include(nfft_programs)  # nfft_add_serial_program / nfft_add_omp_program helpers

# Kernel always built.
add_subdirectory(kernel)

# Optional benchmarks (built unless NFFT_BENCHMARK_MODE is "off").
if(NOT NFFT_BENCHMARK_MODE STREQUAL "off")
  include(nfft_benchmarks) # parse NFFT_AGNOSTIC_BENCHMARKS into gating vars
  add_subdirectory(benchmarks)
endif()

# Optional tests, gracefully skipped if CUnit not found.
if(NFFT_ENABLE_TESTS)
  find_package(CUnit)
  if(CUnit_FOUND)
    enable_testing()
    add_subdirectory(tests)
  else()
    message(STATUS "CUnit not found — skipping the test suite (set NFFT_ENABLE_TESTS=OFF to silence).")
  endif()
endif()

# Optional example programs, never installed.
if(NFFT_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

# Optional applications programs, never installed.
if(NFFT_ENABLE_APPLICATIONS)
  add_subdirectory(applications)
endif()

# Optional language interfaces.
if(NFFT_ENABLE_JULIA OR NFFT_WITH_OCTAVE OR NFFT_WITH_MATLAB)
  include(nfft_interface_kernel)
  nfft_define_interface_kernels()
  include(nfft_interfaces)
endif()

if(NFFT_ENABLE_JULIA)
  add_subdirectory(julia)
endif()

if(NFFT_WITH_OCTAVE OR NFFT_WITH_MATLAB)
  add_subdirectory(matlab)
endif()

include(nfft_install)
