cmake_minimum_required(VERSION 3.20)
project(nfft3_consumer_fc LANGUAGES C)

include(FetchContent)

# Point at the in-repo source tree (passed via -DNFFT3_SOURCE_DIR=...).
# Turn off benchmarks so the sub-build doesn't FetchContent codspeed-cpp,
# and off OpenMP to keep the smoke test minimal.
set(NFFT_ENABLE_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(NFFT_ENABLE_OPENMP     OFF CACHE BOOL "" FORCE)

FetchContent_Declare(nfft3 SOURCE_DIR "${NFFT3_SOURCE_DIR}")
FetchContent_MakeAvailable(nfft3)

add_executable(use_nfft3_fc ${CMAKE_CURRENT_SOURCE_DIR}/../find_package/use_nfft3.c)
target_link_libraries(use_nfft3_fc PRIVATE NFFT3::nfft3)
