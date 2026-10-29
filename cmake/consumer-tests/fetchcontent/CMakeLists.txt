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
project(nfft3_consumer_fc LANGUAGES C)

include(FetchContent)

# Point at the in-repo source tree (passed via -DNFFT3_SOURCE_DIR=...).
# Keep benchmarks off (the default) so the sub-build doesn't FetchContent
# codspeed-cpp, and off OpenMP to keep the smoke test minimal.
set(NFFT_BENCHMARK_MODE off CACHE STRING "" FORCE)
set(NFFT_ENABLE_OPENMP  OFF CACHE BOOL   "" FORCE)

FetchContent_Declare(nfft3 SOURCE_DIR "${NFFT3_SOURCE_DIR}")
FetchContent_MakeAvailable(nfft3)

add_executable(use_nfft3_fc ${CMAKE_CURRENT_SOURCE_DIR}/../find_package/use_nfft3.c)
target_link_libraries(use_nfft3_fc PRIVATE NFFT3::nfft3)
