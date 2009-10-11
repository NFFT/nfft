# $Id$
#
# Copyright (c) 2005, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
################################################################################
# Process this file with autoconf to produce a configure script.
################################################################################

# autoconf initialization
AC_INIT([NFFT],[3.1.2],[potts@mathematik.tu-chemnitz.de])

# revision id from svn
AC_REVISION([$Id$])

# copyright notice
AC_COPYRIGHT([2003, 2009, Jens Keiner, Stefan Kunis, Daniel Potts])

# m4 macros go here
AC_CONFIG_MACRO_DIR([m4])

# where to put auxilliary files
AC_CONFIG_AUX_DIR([config])

# how to recognize the source directory
AC_CONFIG_SRCDIR([include/nfft3.h])

# header to create
AC_CONFIG_HEADERS([include/config.h:include/config.h.in])

# canonical host system type string
AC_CANONICAL_HOST

# substitute abs_srcdir in generated Makefiles
AC_SUBST([abs_srcdir])

# minimum required libtool version
LT_PREREQ([2.2.0])

# libtool initialization
LT_INIT([win32-dll])

# substitute LIBTOOL_DEPS variable in generated Makefiles
AC_SUBST([LIBTOOL_DEPS])

# version information for shared library
SHARED_VERSION_INFO="0:0:0"

# substitute SHARED_VERSION_INFO in generated Makefiles
AC_SUBST(SHARED_VERSION_INFO)

# check for make
AC_PROG_MAKE_SET

# automake initialization
AM_INIT_AUTOMAKE(1.10)

################################################################################
# options for customizing the build process
################################################################################

# maintainer mode option
AM_MAINTAINER_MODE

# enable or disable parts of NFFT

# whether we need the fpt module
need_fpt="no"

# build all modules by default in maintainer mode or if option is given
AC_ARG_ENABLE(all, [AC_HELP_STRING([--enable-all],[build all modules])], 
  ok=$enableval, ok=no)
if test "x$ok" = "xyes" -o "x$USE_MAINTAINER_MODE" = "xyes"; then
  nfft_module_default="yes"
else
  nfft_module_default="no"
fi

# options for modules
AX_NFFT_MODULE([nfct],[NFCT],[nonequispaced fast cosine transform])
AX_NFFT_MODULE([nfst],[NFST],[nonequispaced fast sine transform])
AX_NFFT_MODULE([nfsft],[NFSFT],[nonequispaced fast spherical Fourier transform],
  [need_fpt="yes"])
AX_NFFT_MODULE([nfsoft],[NFSOFT],[nonequispaced fast SO(3) Fourier transform],
  [need_fpt="yes"])
AX_NFFT_MODULE([nnfft],[NNFFT],[nonequispaced fast Fourier transform -- ]#
  [nonequispaced in both domains])
AX_NFFT_MODULE([nsfft],[NSFFT],[nonequispaced sparse fast Fourier transform])
AX_NFFT_MODULE([mri],[MRI],[magnet resonance imaging])
AX_NFFT_MODULE([fpt],[FPT],[fast polynomial transform],[],[],
  ["x$ok" = "xyes" -o "x$need_fpt" = "xyes"])

# multithreaded code
AC_ARG_ENABLE(threads, [AC_HELP_STRING([--enable-threads],
  [enable multithreaded code])], enable_threads=$enableval, enable_threads=no)
AM_CONDITIONAL(HAVE_THREADS, test "x$enable_threads" = "xyes" )

# debug mode
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],
  [compile with extra runtime checks for debugging])], enable_debug=$enableval,
  enable_debug=no)
if test "x$enable_debug" = "xyes"; then
  AC_DEFINE(NFFT_DEBUG,1,[Define to enable extra debugging code.])
fi

# runtime time measurements
AC_ARG_ENABLE(measure-time, [AC_HELP_STRING([--enable-measure-time],
  [measure time during execution])], ok=$enableval, ok=no)
if test "x$ok" = "xyes"; then
  AC_DEFINE(MEASURE_TIME,1,[Define to enable runtime time measurements.])
fi

# runtime time measurements for FFTW part
AC_ARG_ENABLE(measure-time-fftw, [AC_HELP_STRING([--enable-measure-time-fftw],
  [measure time of FFTW transforms during execution])], ok=$enableval, ok=no)
if test "x$ok" = "xyes"; then
  AC_DEFINE(MEASURE_TIME_FFTW,1,[Define to enable time measurements for FFTW]
  [transforms.])
fi

# select window function
AC_ARG_WITH(window, [AC_HELP_STRING([--with-window=ARG],[choose window function
  (ARG can be one of: kaiserbessel (default), gaussian, bspline, sinc,
  dirac)])],
  window=$withval, window="kaiserbessel")

case "$window" in
  gaussian)
    AC_DEFINE(GAUSSIAN,1,[Define to enable Gaussian window function.]);;
  bspline)
    AC_DEFINE(B_SPLINE,1,[Define to enable B-spline window function.]);;
  sinc)
    AC_DEFINE(SINC_POWER,1,[Define to enable sinc power window function.]);;
  delta)
    AC_DEFINE(DIRAC_DELTA,1,[Define to enable dirac delta window function.]);;
  *)
    AC_DEFINE(KAISER_BESSEL,1,[Define to enable Kaiser-Bessel window function.]);;
esac

################################################################################
# compiler characteristis
################################################################################

# select programming language
AC_LANG(C)

# compiler vendor
AX_COMPILER_VENDOR

# check for C99 compliant mode (possibly with GNU extensions)
AC_PROG_CC_C99

# enable "const" keyword
AC_C_CONST

# enable "restrict" keyword
AC_C_RESTRICT

# enable "inline" keyword
AC_C_INLINE

# check if compiler supports blocks
AX_CC_BLOCKS

# reset CC variable (the C99 option is added back below)
CC="$ac_save_CC"

################################################################################
# program checks
################################################################################

# C preprocessor
AC_PROG_CPP_WERROR

# assembler
AM_PROG_AS

# BSD-compatible install
AC_PROG_INSTALL

# whether ln -s works
AC_PROG_LN_S

################################################################################
# 3rd party libraries
################################################################################

# FFTW
AX_LIB_FFTW3

if test "x$ax_lib_fftw3" = "xno"; then
  AC_MSG_ERROR([You don't seem to have the FFTW 3 library installed. You can ]#
     [download it from http://www.fftw.org. If you have installed FFTW 3, ]#
     [make sure that this configure script can find it. See ./configure ]#
     [--help for more information.])
fi

if  test "x$enable_threads" = "xyes" -a "x$ax_lib_fftw3_threads" = "xno"; then
  AC_MSG_ERROR([You don't seem to have the threaded FFTW 3 library installed.])
fi

# libdispatch
AX_LIB_DISPATCH

# Matlab
AX_PROG_MATLAB

################################################################################
# compiler options
################################################################################

# Try to choose good compiler options.
if test "x$ac_test_CFLAGS" != "xset"; then
  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"
  CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"
  LIBS="$LIBS $fftw3_LIBS"
  LDFLAGS="$LDFLAGS $fftw3_LDFLAGS"
  AX_CC_MAXOPT
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"
fi

# override CFLAGS selection when debugging
if test "x${enable_debug}" = "xyes"; then
  CFLAGS="-g"
  if test "x$ax_cv_c_compiler_vendor" = "xapple"; then
    CFLAGS="$CFLAGS $ax_cv_apple_gcc_archflag"
  fi
fi

# add gcc warnings, in debug/maintainer mode only
if test "x${enable_debug}" = "xyes" || test "x$USE_MAINTAINER_MODE" = "xyes"; 
then
  if test "x$ac_cv_c_compiler_gnu" = "xyes"; then
    CFLAGS="$CFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
    CFLAGS="$CFLAGS -Wno-long-long -Wshadow -Wbad-function-cast -Wwrite-strings"
    CFLAGS="$CFLAGS -Wstrict-prototypes -Wredundant-decls -Wnested-externs"
    CFLAGS="$CFLAGS -Wundef -Wconversion -Wmissing-prototypes "
    CFLAGS="$CFLAGS -Wmissing-declarations"
  fi
fi

# option to accept C99
CFLAGS="$CFLAGS $ac_cv_prog_cc_c99" 

CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"

# add Matlab CFLAGS
if test "x$ax_prog_matlab" = "xyes"; then
  CFLAGS="$CFLAGS $matlab_CFLAGS"
fi

################################################################################
# header files/data types/compiler characteristics
################################################################################

AC_CHECK_HEADERS([math.h stdio.h stdlib.h sys/resource.h time.h complex.h \
  string.h float.h limits.h stdarg.h stddef.h sys/types.h stdint.h inttypes.h \
  stdbool.h sys/time.h malloc.h])

AC_HEADER_TIME

AC_TYPE_SIZE_T
AC_CHECK_TYPE([long double],
  [AC_DEFINE(HAVE_LONG_DOUBLE, 1, [Define to 1 if the compiler supports]
  [`long double'])],[])
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t]
  [is defined in <sys/time.h>])],,
  [
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
  ])
AC_CHECK_TYPES(uintptr_t, [], [AC_CHECK_SIZEOF(void *)], [$ac_includes_default
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#endif])

AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(unsigned long long)
AC_CHECK_SIZEOF(size_t)
AC_CHECK_SIZEOF(ptrdiff_t)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long double)

# library functions
AC_FUNC_ALLOCA
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
saved_LIBS=$LIBS
AC_CHECK_LIB(m, sin)
LIBS=$saved_LIBS

AC_CHECK_FUNCS([BSDgettimeofday gettimeofday gethrtime read_real_time])
AC_CHECK_FUNCS([time_base_to_time drand48 memset posix_memalign memalign])
AC_CHECK_FUNCS([_mm_malloc _mm_free clock_gettime mach_absolute_time sysctl])
AC_CHECK_FUNCS([abort snprintf sqrt tgamma lgamma log sin cos acos pow exp])
AC_CHECK_FUNCS([log2 log10 log1p sqrt tgamma lgamma log sin cos acos pow exp copysign lrint cexp])
AC_CHECK_FUNCS([log2f log10f log1pf sqrtf tgammaf lgammaf logf sinf cosf acosf powf expf copysignf lrintf cexpf])
AC_CHECK_FUNCS([log2l log10l log1pl sqrtl tgammal lgammal logl sinl cosl acosl powl expl copysignl lrintl cexpl])
AC_CHECK_FUNCS([nanosleep sleep drand48 srand48])

AC_CHECK_DECLS([memalign, posix_memalign])
AC_CHECK_DECLS([log2, log10, log1p, sqrt, tgamma, lgamma, log, sin, cos, acos, pow, exp, copysign, lrint],[],[],[#include <math.h>])
AC_CHECK_DECLS([log2f, log10f, log1pf, sqrtf, tgammaf, lgammaf, logf, sinf, cosf, acosf, powf, expf, copysignf, lrintf],[],[],[#include <math.h>])
AC_CHECK_DECLS([log2l, log10l, log1pl, sqrtl, tgammal, lgammal, logl, sinl, cosl, acosl, powl, expl, copysignl, lrintl],[],[],[#include <math.h>])
AC_CHECK_DECLS([cexp],[],[],[#include <complex.h>])
AC_CHECK_DECLS([cexpf],[],[],[#include <complex.h>])
AC_CHECK_DECLS([cexpl],[],[],[#include <complex.h>])
AC_CHECK_DECLS([nanosleep],[],[],[#include <time.h>])
AC_CHECK_DECLS([drand48],[],[],[#include <stdlib.h>])
AC_CHECK_DECLS([srand48],[],[],[#include <stdlib.h>])
AC_CHECK_DECLS([sleep],[],[],[#include <unistd.h>])

# The output files to be generated
AC_CONFIG_FILES( \
  Makefile \
  include/Makefile \
  include/nfft3conf.h \
  3rdparty/Makefile \
  3rdparty/cstripack/Makefile \
  util/Makefile \
  kernel/Makefile \
  kernel/fpt/Makefile \
  kernel/mri/Makefile \
  kernel/nfct/Makefile \
  kernel/nfft/Makefile \
  kernel/nfsft/Makefile \
  kernel/nfsoft/Makefile \
  kernel/nfst/Makefile \
  kernel/nnfft/Makefile \
  kernel/nsfft/Makefile \
  kernel/solver/Makefile \
  examples/Makefile \
  examples/fpt/Makefile \
  examples/mri/Makefile \
  examples/nfct/Makefile \
  examples/nfft/Makefile \
  examples/nfsft/Makefile \
  examples/nfsoft/Makefile \
  examples/nfst/Makefile \
  examples/nnfft/Makefile \
  examples/nsfft/Makefile \
  examples/solver/Makefile \
  applications/Makefile \
  applications/fastgauss/Makefile \
  applications/fastsum/Makefile \
  applications/fastsumS2/Makefile \
  applications/quadratureS2/Makefile \
  applications/mri/Makefile \
  applications/mri/mri2d/Makefile \
  applications/mri/mri3d/Makefile \
  applications/polarFFT/Makefile \
  applications/radon/Makefile \
  matlab/Makefile \
  matlab/nfsft/Makefile \
  matlab/nfft/Makefile \
  matlab/nfsft/@f_hat/Makefile 
)

if test "x$USE_MAINTAINER_MODE" = "xyes"; then
AC_CONFIG_FILES( \
  tests/Makefile \
  VERSION
  doxygen/doxygen.Doxyfile
  tests/kernel/Makefile \
  tests/kernel/fpt/Makefile \
  tests/lowlevel/Makefile \
  tests/lowlevel/nfsft/Makefile \
  applications/iterS2/Makefile
)
fi

AC_OUTPUT
