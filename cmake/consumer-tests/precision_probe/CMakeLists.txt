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
