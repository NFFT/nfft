# Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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
# @synopsis NFFT_LIB_FFTW3
# @summary Check configure options and assign variables related to the fftw3 library.
# @category C
#
# @version 2008-12-07
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>
#
#  If we find the library, set the shell variable `ax_have_fftw3' to `yes'.
#  Otherwise, set `ax_have_fftw3' to `no'.

AC_DEFUN([NFFT_LIB_FFTW3],
[
  AC_REQUIRE([AX_OPENMP])

  saved_CPPFLAGS="$CPPFLAGS"
  saved_CFLAGS="$CFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  nfft_fftw3_have_lib="no"
  nfft_fftw3_CPPFLAGS=""
  nfft_fftw3_LDLFAGS=""
  nfft_fftw3_LIBS=""
  nfft_fftw3_have_lib_threads="no"
  nfft_fftw3_LIBS_threads=""
  nfft_fftw3_have_lib_omp="no"
  nfft_fftw3_LIBS_omp=""

  AC_ARG_WITH(fftw3, [AS_HELP_STRING([--with-fftw3=DIR],
  [compile with fftw3 in root directory DIR])], with_fftw3=$withval, with_fftw3="yes")

  AC_ARG_WITH(fftw3-libdir, [AS_HELP_STRING([--with-fftw3-libdir=DIR],
  [compile with fftw3 library directory DIR])], fftw3_lib_dir=$withval, 
    fftw3_lib_dir="yes")

  AC_ARG_WITH(fftw3-includedir, [AS_HELP_STRING([--with-fftw3-includedir=DIR],
  [compile with fftw3 include directory DIR])], fftw3_include_dir=$withval, 
    fftw3_include_dir="yes")

  if test "x$with_fftw3" != "xyes"; then
    if test "x${fftw3_include_dir}" = "xyes"; then
      if test -d "$with_fftw3/include"; then
        fftw3_include_dir="$with_fftw3/include"
      elif test -d "$with_fftw3/api"; then
        fftw3_include_dir="$with_fftw3/api"
      else
        fftw3_include_dir="$with_fftw3"
      fi
    fi
    if test "x${fftw3_lib_dir}" = "xyes"; then 
      if test -d "$with_fftw3/lib"; then
        fftw3_lib_dir="$with_fftw3/lib"
      elif test -d "$with_fftw3/.libs"; then
        fftw3_lib_dir="$with_fftw3/.libs"
      else
        fftw3_lib_dir="$with_fftw3"
      fi
    fi
  fi

  if test "x${fftw3_include_dir}" != "xyes"; then
    AX_CHECK_DIR([${fftw3_include_dir}],[],
      [AC_MSG_ERROR([The FFTW3 include directory "${fftw3_include_dir}" does not exist.])])
    nfft_fftw3_CPPFLAGS="-I$fftw3_include_dir"
  fi

  if test "x${fftw3_lib_dir}" != "xyes"; then 
    AX_CHECK_DIR([${fftw3_lib_dir}],[],
      [AC_MSG_ERROR([The FFTW3 library directory "${fftw3_lib_dir}" does not exist.])])
    nfft_fftw3_LDFLAGS="-L$fftw3_lib_dir"
  fi

  CPPFLAGS="$nfft_fftw3_CPPFLAGS $CPPFLAGS"

  # Check if header is present and usable.
  AC_CHECK_HEADER([fftw3.h], [nfft_lib_fftw3_have_header=yes], [nfft_lib_fftw3_have_header=no])

  if test "x$nfft_lib_fftw3_have_header" = "xyes"; then
    LDFLAGS="$nfft_fftw3_LDFLAGS $saved_LDFLAGS"

    # Check if uniprocessor FFTW3 library is present and usable.
    LIBS="$saved_LIBS"
    AS_UNSET([ac_cv_search_fftw${PREC_SUFFIX}_execute])
    AC_SEARCH_LIBS([fftw${PREC_SUFFIX}_execute], [fftw3${PREC_SUFFIX} fftw3${PREC_SUFFIX}-3], [], [], [-lm])

    if test "x$ac_res" != "xno"; then
      nfft_fftw3_have_lib="yes"
      nfft_fftw3_LIBS="-lm"
      if test "x$ac_res" != "xnone required"; then
        nfft_fftw3_LIBS="$ac_res $nfft_fftw3_LIBS"
      fi
    fi

    if test "x$nfft_fftw3_have_lib" = "xyes"; then
      # Check if generic multi-threaded FFTW3 library is present and usable without linking any specific threads implementation.
      LIBS="${nfft_fftw3_LIBS} $saved_LIBS"
      AS_UNSET([ac_cv_search_fftw${PREC_SUFFIX}_init_threads])
      AC_SEARCH_LIBS([fftw${PREC_SUFFIX}_init_threads], [fftw3${PREC_SUFFIX}_threads fftw3${PREC_SUFFIX}_threads-3], [], [], [])

      if test "x$ac_res" != "xno"; then
        nfft_fftw3_have_lib_threads="yes"
        if test "x$ac_res" != "xnone required"; then
          nfft_fftw3_LIBS_threads="$ac_res"
        fi
      fi

      if test "x$nfft_fftw3_have_lib_threads" = "xno"; then
        # Check if generic multi-threaded FFTW3 library is present and usable when linking with -lpthread.
        LIBS="${nfft_fftw3_LIBS} $saved_LIBS"
        AS_UNSET([ac_cv_search_fftw${PREC_SUFFIX}_init_threads])
        AC_SEARCH_LIBS([fftw${PREC_SUFFIX}_init_threads], [fftw3${PREC_SUFFIX}_threads fftw3${PREC_SUFFIX}_threads-3], [], [], [-lpthread])

        if test "x$ac_res" != "xno"; then
          nfft_fftw3_have_lib_threads="yes"
          nfft_fftw3_LIBS_threads="-lpthread"
          if test "x$ac_res" != "xnone required"; then
            nfft_fftw3_LIBS_threads="$ac_res $nfft_fftw3_LIBS_threads"
          fi
        fi
      fi

      # Check if multi-threaded FFTW3 library with OpenMP is present and usable.
      LIBS="${nfft_fftw3_LIBS} $saved_LIBS"
      CFLAGS="$OPENMP_CFLAGS $saved_CFLAGS"
      LDFLAGS="$OPENMP_CFLAGS $saved_LDFLAGS"
      AS_UNSET([ac_cv_search_fftw${PREC_SUFFIX}_init_threads])
      AC_SEARCH_LIBS([fftw${PREC_SUFFIX}_init_threads], [fftw3${PREC_SUFFIX}_omp fftw3${PREC_SUFFIX}_omp-3], [], [], [])

      if test "x$ac_res" != "xno"; then
        nfft_fftw3_have_lib_omp="yes"
        if test "x$ac_res" != "xnone required"; then
          nfft_fftw3_LIBS_omp="$ac_res"
        fi
      fi
    fi
  fi

  # Kept for compatibility with depdent macros (for now),
  fftw3_CPPFLAGS=$nfft_fftw3_CPPFLAGS
  fftw3_LDFLAGS=$nfft_fftw3_LDFLAGS
  fftw3_LIBS=$nfft_fftw3_LIBS
  fftw3_threads_LIBS="$nfft_fftw3_LIBS_threads $nfft_fftw3_LIBS"

  # Restore saved flags.
  CPPFLAGS="$saved_CPPFLAGS"
  CFLAGS="$saved_CFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"
])
