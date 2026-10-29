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

# CUnit test suite.

# Shared source list.
set(NFFT_TEST_SOURCES
  check.c
  util.c    util.h
  accuracy_log.c accuracy_log.h
  reflect.c reflect.h
  bspline.c bspline.h
  bessel.c  bessel.h
  window.c  window.h
  kbpoly.c  kbpoly.h
  nfft.c    nfft.h
  nfct.c    nfct.h
  nfst.c    nfst.h)

  # Non-threaded test program.
add_executable(checkall ${NFFT_TEST_SOURCES})

# Generated includes first (build/), then source includes (include/).
target_include_directories(checkall PRIVATE
  ${NFFT_GENERATED_INCLUDE_DIR}      # config.h, ticks.h
  ${PROJECT_SOURCE_DIR}/include)     # nfft3.h, infft.h
target_link_libraries(checkall PRIVATE nfft3 FFTW3::fftw3 CUnit::CUnit)

if(UNIX)
# Add libm.
  target_link_libraries(checkall PRIVATE m)
endif()

add_test(NAME checkall COMMAND checkall)

# Optional threaded test program.
if(NFFT_ENABLE_OPENMP)
  add_executable(checkall_threads ${NFFT_TEST_SOURCES})
  target_include_directories(checkall_threads PRIVATE
    ${NFFT_GENERATED_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(checkall_threads PRIVATE
    nfft3_omp OpenMP::OpenMP_C FFTW3::fftw3 ${NFFT_OPENMP_ATOMIC_LIBS} CUnit::CUnit)
  if(FFTW3_OMP_FOUND)
    target_link_libraries(checkall_threads PRIVATE FFTW3::fftw3_omp)
  elseif(FFTW3_THREADS_FOUND)
    target_link_libraries(checkall_threads PRIVATE FFTW3::fftw3_threads)
  endif()
  if(UNIX)
    target_link_libraries(checkall_threads PRIVATE m)
  endif()
  add_test(NAME checkall_threads COMMAND checkall_threads)
endif()
