# Gate bench_nfft_direct on BENCHMARK_AGNOSTIC_WINDOW (mirrors benchmarks/Makefile.am:4).
# Returning early also skips the codspeed FetchContent, so non-kaiserbessel cells
# pay no benchmark cost.
if(NOT BENCHMARK_AGNOSTIC_WINDOW)
  message(STATUS "bench_nfft_direct disabled (BENCHMARK_AGNOSTIC_WINDOW=OFF)")
  return()
endif()

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# Codspeed version.
set(CODSPEED_VERSION "v2.3.0" CACHE STRING "codspeed-cpp tag")

# Build the CodSpeed google_benchmark fork with our toolchain.
set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_TESTING        OFF CACHE BOOL "" FORCE)
# google_benchmark builds with -Werror by default. Clang's -Wsign-conversion flags a
# signedness conversion in codspeed-cpp's own core/include/measurement.hpp (which gcc
# does not), breaking the codspeed library build in every clang cell. Disable Werror
# for the FetchContent sub-build only; this does not touch NFFT's own warnings.
set(BENCHMARK_ENABLE_WERROR         OFF CACHE BOOL "" FORCE)

FetchContent_Declare(codspeed
  GIT_REPOSITORY https://github.com/CodSpeedHQ/codspeed-cpp
  GIT_TAG        ${CODSPEED_VERSION}
  GIT_SUBMODULES_RECURSE TRUE
  EXCLUDE_FROM_ALL                            # keep codspeed's install() rules
                                              # (which hardcode DESTINATION /usr/local)
                                              # out of NFFT's `cmake --install`
  SOURCE_SUBDIR  google_benchmark)            # add only the google_benchmark dir
# -DFETCHCONTENT_SOURCE_DIR_CODSPEED=/path/to/codspeed-cpp uses the provided checkout 
# instead of downloading. Useful in CI to avoid network issues.
#
# Force a static build for this dependency since we default to BUILD_SHARED_LIBS=ON 
# globally. This would yield a shared ibbenchmark with unresolved core symbols.
set(_nfft_saved_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF)
FetchContent_MakeAvailable(codspeed)
set(BUILD_SHARED_LIBS ${_nfft_saved_BUILD_SHARED_LIBS})

# Pin the CodSpeed "<file>::" name prefix deterministically. codspeed-cpp derives
# it from std::filesystem::relative(__FILE__, CODSPEED_ROOT_DIR); with the absolute
# __FILE__ CMake passes and CODSPEED_ROOT_DIR at the repo root, it resolves to
# "benchmarks/bench_nfft_direct.cpp" regardless of the run-time working directory.
# CODSPEED_ROOT_DIR is read inside the codspeed core, so it must be set on the
# `codspeed` target, not on bench_nfft_direct.
target_compile_definitions(codspeed PRIVATE
  "CODSPEED_ROOT_DIR=\"${PROJECT_SOURCE_DIR}\"")

add_executable(bench_nfft_direct bench_nfft_direct.cpp util.h)
target_include_directories(bench_nfft_direct PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${NFFT_GENERATED_INCLUDE_DIR})              # config.h (with BENCHMARKS_PREFIX)
target_link_libraries(bench_nfft_direct PRIVATE nfft3 FFTW3::fftw3 benchmark::benchmark)
if(UNIX)
  target_link_libraries(bench_nfft_direct PRIVATE m)
endif()
