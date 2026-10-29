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

# FindFFTW3 — locate FFTW3 for the selected precision (NFFT_PREC_SUFFIX) plus
# the optional _omp and _threads variants. Honors FFTW3_ROOT (var or env) and
# FFTW3_INCLUDEDIR / FFTW3_LIBDIR.
#
# Provides imported targets:
#   FFTW3::fftw3          (required)
#   FFTW3::fftw3_omp      (if found)
#   FFTW3::fftw3_threads  (if found)
# and sets FFTW3_FOUND, FFTW3_OMP_FOUND, FFTW3_THREADS_FOUND.

set(_s "${NFFT_PREC_SUFFIX}")

find_path(FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ${FFTW3_ROOT} ENV FFTW3_ROOT ${FFTW3_INCLUDEDIR}
  PATH_SUFFIXES include)

find_library(FFTW3_LIBRARY
  NAMES fftw3${_s} fftw3${_s}-3 libfftw3${_s}
  HINTS ${FFTW3_ROOT} ENV FFTW3_ROOT ${FFTW3_LIBDIR}
  PATH_SUFFIXES lib lib64)

find_library(FFTW3_OMP_LIBRARY
  NAMES fftw3${_s}_omp fftw3${_s}_omp-3
  HINTS ${FFTW3_ROOT} ENV FFTW3_ROOT ${FFTW3_LIBDIR}
  PATH_SUFFIXES lib lib64)

find_library(FFTW3_THREADS_LIBRARY
  NAMES fftw3${_s}_threads fftw3${_s}_threads-3
  HINTS ${FFTW3_ROOT} ENV FFTW3_ROOT ${FFTW3_LIBDIR}
  PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
  REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR)

if(FFTW3_FOUND AND NOT TARGET FFTW3::fftw3)
  add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
  set_target_properties(FFTW3::fftw3 PROPERTIES
    IMPORTED_LOCATION "${FFTW3_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}")
endif()

if(FFTW3_OMP_LIBRARY AND NOT TARGET FFTW3::fftw3_omp)
  set(FFTW3_OMP_FOUND TRUE)
  add_library(FFTW3::fftw3_omp UNKNOWN IMPORTED)
  set_target_properties(FFTW3::fftw3_omp PROPERTIES
    IMPORTED_LOCATION "${FFTW3_OMP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}")
endif()

if(FFTW3_THREADS_LIBRARY AND NOT TARGET FFTW3::fftw3_threads)
  set(FFTW3_THREADS_FOUND TRUE)
  add_library(FFTW3::fftw3_threads UNKNOWN IMPORTED)
  set_target_properties(FFTW3::fftw3_threads PROPERTIES
    IMPORTED_LOCATION "${FFTW3_THREADS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}")
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY FFTW3_OMP_LIBRARY FFTW3_THREADS_LIBRARY)
