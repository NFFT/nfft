cmake_minimum_required(VERSION 3.20)
project(nfft3_precision_probe C)

# Driven by a single -DNFFT_PROBE_SUFFIX=<f|l|> (empty = double).
set(NFFT_PROBE_SUFFIX "" CACHE STRING "Precision suffix to probe (f, l, or empty)")

set(_pkg "NFFT3${NFFT_PROBE_SUFFIX}")
set(_tgt "NFFT3::nfft3${NFFT_PROBE_SUFFIX}")
set(_pfx "nfft${NFFT_PROBE_SUFFIX}_")

find_package(${_pkg} REQUIRED)
if(NOT TARGET ${_tgt})
  message(FATAL_ERROR "Expected imported target ${_tgt} from package ${_pkg}, but it is missing")
endif()
message(STATUS "Probe: ${_pkg} found, target ${_tgt} present, API prefix '${_pfx}'")

add_executable(use_nfft3_precision use_nfft3_precision.c)
target_compile_definitions(use_nfft3_precision PRIVATE NFFT_PREFIX=${_pfx})
target_link_libraries(use_nfft3_precision PRIVATE ${_tgt})

if(NFFT_PROBE_SUFFIX STREQUAL "f")
  set(_prec NFFT_PRECISION_SINGLE)
elseif(NFFT_PROBE_SUFFIX STREQUAL "l")
  set(_prec NFFT_PRECISION_LONG_DOUBLE)
else()
  set(_prec NFFT_PRECISION_DOUBLE)
endif()
add_executable(use_nfft3_nocomplex ../find_package/use_nfft3_nocomplex.c)
target_compile_definitions(use_nfft3_nocomplex PRIVATE ${_prec})
target_link_libraries(use_nfft3_nocomplex PRIVATE ${_tgt})
