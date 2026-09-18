# Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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
# @synopsis NFFT_OPENMP_MATCH_FFTW3
# @summary Make NFFT link the same OpenMP runtime as the OpenMP FFTW3 library.
# @category C
#
# @version 2026-09-18
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>
#
#  On macOS, two libomp.dylib copies at different paths can both load into one
#  process, which aborts with "OMP: Error #15". Homebrew's llvm ships its own
#  libomp while Homebrew's fftw links the separate libomp formula. If the libomp
#  we would link differs from the one libfftw3_omp uses, prepend the directory of
#  the latter to LDFLAGS. Then run a program that starts both runtimes and fail
#  if it aborts. Expects fftw3_LDFLAGS, fftw3_LIBS and fftw3_LIBS_omp to be set.
AC_DEFUN([NFFT_OPENMP_MATCH_FFTW3],
[
  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  CFLAGS="$CFLAGS $OPENMP_CFLAGS $fftw3_CPPFLAGS"
  LIBS="$fftw3_LIBS_omp $fftw3_LIBS $OPENMP_LIBS $LIBS"

  m4_define([NFFT_OPENMP_FFTW3_TEST_PROGRAM], [
    AC_LANG_PROGRAM([
      #include <omp.h>
      #include <fftw3.h>
    ], [
      int n = 0;
      fftw${PREC_SUFFIX}_complex *x;
      fftw${PREC_SUFFIX}_plan p;
      #pragma omp parallel reduction(+:n)
      n += 1;
      if (n < 1 || !fftw${PREC_SUFFIX}_init_threads())
        return 1;
      fftw${PREC_SUFFIX}_plan_with_nthreads(2);
      x = fftw${PREC_SUFFIX}_alloc_complex(64);
      p = fftw${PREC_SUFFIX}_plan_dft_1d(64, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw${PREC_SUFFIX}_execute(p);
      fftw${PREC_SUFFIX}_destroy_plan(p);
      fftw${PREC_SUFFIX}_free(x);
      fftw${PREC_SUFFIX}_cleanup_threads();
      return 0;
    ])
  ])

  case $host_os in
    darwin*)
      AC_CHECK_PROG([OTOOL], [otool], [otool])
      if test "x$OTOOL" != "x"; then
        AC_MSG_CHECKING([whether NFFT and FFTW3 link the same libomp])
        nfft_omp_ours=""
        nfft_omp_fftw=""
        LDFLAGS_bak="$LDFLAGS"
        LDFLAGS="$OPENMP_CFLAGS $fftw3_LDFLAGS $LDFLAGS"
        AC_LINK_IFELSE([NFFT_OPENMP_FFTW3_TEST_PROGRAM], [
          nfft_omp_ours=`$OTOOL -L conftest$ac_exeext | sed -n 's|^[[	 ]]*\([[^ ]]*/libomp\.dylib\) .*|\1|p'`
          nfft_fftw_omp_lib=`$OTOOL -L conftest$ac_exeext | sed -n "s|^[[	 ]]*\([[^ ]]*libfftw3${PREC_SUFFIX}_omp[[^ /]]*\.dylib\) .*|\1|p"`
          if test -f "$nfft_fftw_omp_lib"; then
            nfft_omp_fftw=`$OTOOL -L "$nfft_fftw_omp_lib" | sed -n 's|^[[	 ]]*\([[^ ]]*/libomp\.dylib\) .*|\1|p'`
          fi
        ])
        LDFLAGS="$LDFLAGS_bak"
        if test "x$nfft_omp_ours" = "x" -o "x$nfft_omp_fftw" = "x"; then
          AC_MSG_RESULT([unknown])
        elif test "x$nfft_omp_ours" = "x$nfft_omp_fftw"; then
          AC_MSG_RESULT([yes])
        elif test -f "$nfft_omp_fftw"; then
          nfft_omp_fftw_dir=`dirname "$nfft_omp_fftw"`
          LDFLAGS="-L$nfft_omp_fftw_dir $LDFLAGS"
          AC_MSG_RESULT([no, prepending -L$nfft_omp_fftw_dir to LDFLAGS])
        else
          AC_MSG_RESULT([no, FFTW3 uses $nfft_omp_fftw])
        fi
      fi
      ;;
  esac

  AC_MSG_CHECKING([whether NFFT and FFTW3 OpenMP code can run in one process])
  LDFLAGS_bak="$LDFLAGS"
  LDFLAGS="$OPENMP_CFLAGS $fftw3_LDFLAGS $LDFLAGS"
  AC_RUN_IFELSE([NFFT_OPENMP_FFTW3_TEST_PROGRAM], [AC_MSG_RESULT([yes])], [
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([A program that uses both OpenMP and the OpenMP FFTW3 library failed to run.]
[This is most likely because two different OpenMP runtimes were loaded, e.g. libomp from]
[Homebrew's llvm and from Homebrew's libomp formula. Make sure the linker finds the OpenMP]
[runtime that FFTW3 uses first, e.g. by putting its directory first in LDFLAGS. See config.log.])
  ], [AC_MSG_RESULT([skipped, cross-compiling])])
  LDFLAGS="$LDFLAGS_bak"

  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"
])
