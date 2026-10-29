# Copyright (c) 2026 The NFFT Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/>.

# nfft_benchmarks.cmake 
# Parses NFFT_AGNOSTIC_BENCHMARKS, a comma-separated
# list of param:flag pairs (param in {window,openmp,precision}, flag in {0,1}),
# into three booleans used to gate which benchmark variants are built:
#   BENCHMARK_AGNOSTIC_WINDOW / BENCHMARK_AGNOSTIC_OPENMP / BENCHMARK_AGNOSTIC_PRECISION
# Empty string  -> all ON (default)
# Non-empty      -> all OFF, then set per flag
# Unknown flag   -> error

set(NFFT_AGNOSTIC_BENCHMARKS "" CACHE STRING
    "Comma-separated param:flag pairs (window/openmp/precision : 0/1) selecting agnostic benchmark variants")

if(NFFT_AGNOSTIC_BENCHMARKS STREQUAL "")
  set(BENCHMARK_AGNOSTIC_WINDOW    ON)
  set(BENCHMARK_AGNOSTIC_OPENMP    ON)
  set(BENCHMARK_AGNOSTIC_PRECISION ON)
else()
  set(BENCHMARK_AGNOSTIC_WINDOW    OFF)
  set(BENCHMARK_AGNOSTIC_OPENMP    OFF)
  set(BENCHMARK_AGNOSTIC_PRECISION OFF)
  string(REPLACE "," ";" _nfft_agn_list "${NFFT_AGNOSTIC_BENCHMARKS}")
  foreach(_pair IN LISTS _nfft_agn_list)
    if(_pair STREQUAL "")
      # skip empty token (mirrors the '') case in configure.ac:290
    elseif(_pair STREQUAL "window:1")
      set(BENCHMARK_AGNOSTIC_WINDOW ON)
    elseif(_pair STREQUAL "window:0")
      set(BENCHMARK_AGNOSTIC_WINDOW OFF)
    elseif(_pair STREQUAL "openmp:1")
      set(BENCHMARK_AGNOSTIC_OPENMP ON)
    elseif(_pair STREQUAL "openmp:0")
      set(BENCHMARK_AGNOSTIC_OPENMP OFF)
    elseif(_pair STREQUAL "precision:1")
      set(BENCHMARK_AGNOSTIC_PRECISION ON)
    elseif(_pair STREQUAL "precision:0")
      set(BENCHMARK_AGNOSTIC_PRECISION OFF)
    else()
      message(FATAL_ERROR
        "Unknown agnostic benchmark flag \"${_pair}\". Valid flags are: precision:0/1, window:0/1, openmp:0/1")
    endif()
  endforeach()
endif()

message(STATUS "Agnostic benchmarks: window=${BENCHMARK_AGNOSTIC_WINDOW} "
               "openmp=${BENCHMARK_AGNOSTIC_OPENMP} precision=${BENCHMARK_AGNOSTIC_PRECISION}")
