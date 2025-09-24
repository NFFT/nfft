# Copyright (c) 2025 Jens Keiner, Stefan Kunis, Daniel Potts
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# @synopsis NFFT_OPENMP_ATOMIC_FLOAT
# @summary Check whether using "#pragma omp atomic" with the selected floating-point data type requires linking additional libraries.
# @category C
#
# @version 2025-09-23
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>
#
#  Set shell variable `OPENMP_ATOMIC_FLOAT_LIBS' to the additional libraries needed for OpenMP
#  atomic operations on floating point type to work. If no valid configuration is found, set
#  `OPENMP_ATOMIC_FLOAT_LIBS' to "unknown".
AC_DEFUN([NFFT_OPENMP_ATOMIC_FLOAT],
[
  AC_MSG_CHECKING([which additional libraries are needed for OpenMP atomic operations on floating point type to work])
  
  # Save current flags.
  ac_save_CFLAGS="$CFLAGS"
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_LIBS="$LIBS"

  # Add OpenMP flags.
  CFLAGS="$CFLAGS $OPENMP_CFLAGS"
  LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"

  # Add OpenMP flags.
  #if test "x$CC" = "xclang"; then
  #  CFLAGS="$CFLAGS -g -Rpass=.* -Rpass-missed=.* -fsave-optimization-record -foptimization-record-file=opt.txt $OPENMP_CFLAGS"
  #  LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"
  #else
  #  CFLAGS="$CFLAGS -g -fopt-info-all=opt.txt $OPENMP_CFLAGS"
  #  LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"
  #fi
  
  # Define test program as a macro for reuse.
  m4_define([OMP_ATOMIC_TEST_PROGRAM], [
    AC_LANG_PROGRAM([
      #include <omp.h>
      #include <stdint.h>
    ], [
      #if defined(NFFT_LDOUBLE)
      long double r = 0.0L;
      long double increment = 2.0L;
      long double expected_min = 150.0L;
      #elif defined(NFFT_SINGLE)
      float r = 0.0f;
      float increment = 2.0f;
      float expected_min = 150.0f;
      #else
      double r = 0.0;
      double increment = 2.0;
      double expected_min = 150.0;
      #endif
      int i;
      // Make loop count unpredictable to compiler by using address arithmetic.
      volatile int loop_count = (int)((uintptr_t)&r % 50) + 75;
      // Make increment unpredictable by using address-based modification.
      increment += (((uintptr_t)&i % 3) == 0) ? 0.0 : 0.0;

      #pragma omp parallel for
      for (i = 0; i < loop_count; i++) {
        #pragma omp atomic
        r += increment;
      }
      
      if (r < expected_min) {
        return 1;
      }
      
      return 0;
    ])
  ])

  # Try different LDFLAGS options for atomic operations.
  OPENMP_ATOMIC_FLOAT_LIBS="unknown"
  for atomic_flag in "" "-latomic"; do
    if test "x$OPENMP_ATOMIC_FLOAT_LIBS" = "xunknown"; then
      LIBS_bak=$LIBS
      LIBS="$LIBS $atomic_flag"
      
      AC_LINK_IFELSE([
        OMP_ATOMIC_TEST_PROGRAM
      ], [
        OPENMP_ATOMIC_FLOAT_LIBS="$atomic_flag"
      ], [
        # Test failed, continue to next iteration
        :
      ])

      # Restore flags.
      LIBS="$LIBS_bak"
    fi
  done

  if test "x$OPENMP_ATOMIC_FLOAT_LIBS" = "xunknown"; then
    AC_MSG_RESULT([unknown])
  else
    if test "x$OPENMP_ATOMIC_FLOAT_LIBS" = "x"; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([$OPENMP_ATOMIC_FLOAT_LIBS])
    fi
  fi
  
  # Restore original flags.
  CFLAGS="$ac_save_CFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"
])
