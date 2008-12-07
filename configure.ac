# $Id: nfsftmex.c 1548 2007-04-18 07:34:23Z keiner $
#
# Copyright (c) 2005, 2008 Jens Keiner, Stefan Kunis, Daniel Potts
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

# initialize autoconf
AC_INIT([NFFT],[3.1.0-alpha],[potts@mathematik.tu-chemnitz.de])

# m4 macros
AC_CONFIG_MACRO_DIR([m4])

# copyright notice
AC_COPYRIGHT([2003, 2008, Jens Keiner, Stefan Kunis, Daniel Potts])

# directory for auxilliary files
AC_CONFIG_AUX_DIR([config])

# minimum libtool version
LT_PREREQ([2.2.0])

# initialize libtool
LT_INIT([win32-dll])

AC_SUBST([LIBTOOL_DEPS])

# revision number
AC_REVISION([$Id$])

# Version information for library file
SHARED_VERSION_INFO="0:0:0"
AC_SUBST(SHARED_VERSION_INFO)

# how to recognize configure source directory
AC_CONFIG_SRCDIR(include/nfft3.h)

# config headers
AC_CONFIG_HEADERS([include/config.h:include/config.h.in])

# check for make
AC_PROG_MAKE_SET

# initialise automake
AM_INIT_AUTOMAKE(1.10)

# host system type
AC_CANONICAL_HOST

# make abs_srcdir available in Makefiles
AC_SUBST([abs_srcdir])

# ############################################################################
# options
# ############################################################################

# enable or disable parts of nfft

need_fpt="no"

AC_ARG_ENABLE(nfct, [AC_HELP_STRING([--enable-nfct],[compile with nfct module (nonequispaced fast cosine transform)])], ok=$enableval, ok=yes)
if test "x$ok" = "xyes"; then
  HAVE_NFCT="#define HAVE_NFCT 1"
  AM_CONDITIONAL(HAVE_NFCT, test "yes" = "yes" )
else
  HAVE_NFCT="#undef HAVE_NFCT"
  AM_CONDITIONAL(HAVE_NFCT, test "yes" = "no" )
fi
AC_SUBST(HAVE_NFCT)

AC_ARG_ENABLE(nfst, [AC_HELP_STRING([--enable-nfst],[compile with nfst module (nonequispaced fast sine transform)])], ok=$enableval, ok=yes)
if test "x$ok" = "xyes"; then
  HAVE_NFST="#define HAVE_NFST 1"
  AM_CONDITIONAL(HAVE_NFST, test "yes" = "yes" )
else
  HAVE_NFST="#undef HAVE_NFST"
  AM_CONDITIONAL(HAVE_NFST, test "yes" = "no" )
fi
AC_SUBST(HAVE_NFST)

AC_ARG_ENABLE(nfsft, [AC_HELP_STRING([--enable-nfsft],[compile with nfsft module (nonequispaced fast spherical Fourier transform; needs fpt module)])], ok=$enableval, ok=yes)
if test "x$ok" = "xyes"; then
  need_fpt="yes"
  HAVE_NFSFT="#define HAVE_NFSFT 1"
  AM_CONDITIONAL(HAVE_NFSFT, test "yes" = "yes" )
else
  HAVE_NFSFT="#undef HAVE_NFSFT"
  AM_CONDITIONAL(HAVE_NFSFT, test "yes" = "no" )
fi
AC_SUBST(HAVE_NFSFT)

AC_ARG_ENABLE(nnfft, [AC_HELP_STRING([--enable-nnfft],[compile with nnfft module (nonequispaced fast Fourier transform in time and frequency)])], ok=$enableval, ok=yes)
if test "x$ok" = "xyes"; then
  HAVE_NNFFT="#define HAVE_NNFFT 1"
  AM_CONDITIONAL(HAVE_NNFFT, test "yes" = "yes" )
else
  HAVE_NNFFT="#undef HAVE_NNFFT"
  AM_CONDITIONAL(HAVE_NNFFT, test "yes" = "no" )
fi
AC_SUBST(HAVE_NNFFT)

AC_ARG_ENABLE(nsfft, [AC_HELP_STRING([--enable-nsfft],[compile with nsfft module (nonequispaced sparse fast Fourier transform)])], ok=$enableval, ok=yes)
if test "x$ok" = "xyes"; then
  HAVE_NSFFT="#define HAVE_NSFFT 1"
  AM_CONDITIONAL(HAVE_NSFFT, test "yes" = "yes" )
else
  HAVE_NSFFT="#undef HAVE_NSFFT"
  AM_CONDITIONAL(HAVE_NSFFT, test "yes" = "no" )
fi
AC_SUBST(HAVE_NSFFT)

AC_ARG_ENABLE(mri, [AC_HELP_STRING([--enable-mri],[compile with mri module (magnet resonance imaging)])], ok=$enableval, ok=yes)
if test "x$ok" = "xyes"; then
  HAVE_MRI="#define HAVE_MRI 1"
  AM_CONDITIONAL(HAVE_MRI, test "yes" = "yes" )
else
  HAVE_MRI="#undef HAVE_MRI"
  AM_CONDITIONAL(HAVE_MRI, test "yes" = "no" )
fi
AC_SUBST(HAVE_MRI)

AC_ARG_ENABLE(fpt, [AC_HELP_STRING([--enable-fpt],[compile with fpt module (fast polynomial transform)])], ok=$enableval, ok=yes)
if test "x$ok" = "xyes"; then
  have_fpt="yes"
  HAVE_FPT="#define HAVE_FPT 1"
  AM_CONDITIONAL(HAVE_FPT, test "yes" = "yes" )
else
  have_fpt="no"
  HAVE_FPT="#undef HAVE_FPT"
  AM_CONDITIONAL(HAVE_FPT, test "yes" = "no" )
fi
AC_SUBST(HAVE_FPT)

if test "x$need_fpt" = "xyes" -a "x$have_fpt" = "xno"; then
  AC_MSG_ERROR([You must enable the FPT module for this configuration (--enable-fpt)])
fi

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile nfft with extra runtime checks for debugging])], ok=$enableval, ok=no)
if test "x$ok" = "xyes"; then
  enable_debug="yes"
  AC_DEFINE(NFFT_DEBUG,1,[Define to enable extra NFFT debugging code.])
fi

AC_ARG_ENABLE(measure-time, [AC_HELP_STRING([--enable-measure-time],
  [measure time during execution])], ok=$enableval, ok=no)
if test "x$ok" = "xyes"; then
  AC_DEFINE(MEASURE_TIME,1,[Define to enable time measurements during execution.])
fi

AC_ARG_ENABLE(measure-time-fftw, [AC_HELP_STRING([--enable-measure-time-fftw],
  [measure time of FFTW transforms during execution])], ok=$enableval, ok=no)
if test "x$ok" = "xyes"; then
  AC_DEFINE(MEASURE_TIME_FFTW,1,[Define to enable time measurements for FFTW transforms.])
fi

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

# maintainer mode option
AM_MAINTAINER_MODE

# ############################################################################
# Checks for programs.
# ############################################################################

# Set programming language to C.
AC_LANG(C)

# Search for C compiler. Try to select C99 compliant mode.
AC_PROG_CC_C99

# Reset CC variable.
CC="$ac_save_CC"

if test "x$ac_cv_prog_cc_c99" = "x-std=gnu99"; then
  ac_cv_prog_cc_c99="-std=c99"
fi

# Make sure "const" keyword can be used.
AC_C_CONST

# Make sure "restrict" keyword can be used.
AC_C_RESTRICT

# Make sure "inline" keyword can be used.
AC_C_INLINE

# Determine compiler vendor (thanks fftw)
AX_COMPILER_VENDOR

# Check for assembler (thanks fftw).
AMX_PROG_AS # AM_PROG_AS breaks at least one version of libtool.

# Search for the C preprocessor.
AC_PROG_CPP_WERROR

# Check for BSD-compatible install program. Enforces inclusion of own install
# script.
AC_PROG_INSTALL

# Check if ln -s works.
AC_PROG_LN_S

# Check for fftw3.
AX_LIB_FFTW3

if test "x$ax_lib_fftw3" = "xno"; then
  AC_MSG_ERROR(
    [The fftw3 library could not be found. You can download it from
     http://www.fftw.org. If you have installed fftw3, make sure that
     this configure script can find it. See ./configure --help for
     more information.])
fi

# Check for MATLAB.
AX_PROG_MATLAB
AM_CONDITIONAL(HAVE_MATLAB, test "x$ax_prog_matlab" = "xyes" )
AC_SUBST(matlab_CPPFLAGS)
AC_SUBST(matlab_LIBS)
AC_SUBST(matlab_LDFLAGS)
AC_SUBST(matlab_mexext)

# Try to choose good compiler options.
if test "x$ac_test_CFLAGS" != "xset"; then
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"
  LIBS="$LIBS $fftw3_LIBS $matlab_LIBADD"
  LDFLAGS="$LDFLAGS $fftw3_LDFLAGS $matlab_LDFLAGS"
  AX_CC_MAXOPT
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"
fi

if test "x$ax_lib_fftw3_libtool" = "xyes"; then
  fftw3_LIBS="$fftw3_lib_dir/libfftw3.la"
fi

AC_SUBST(fftw3_LIBS)
AC_SUBST(fftw3_LDFLAGS)

CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"

# override CFLAGS selection when debugging
if test "x${enable_debug}" = "xyes"; then
  CFLAGS="-g3"
fi

# add gcc warnings, in debug/maintainer mode only
if test "x$enable_debug" = "xyes" || test "x$USE_MAINTAINER_MODE" = "xyes"; then
  if test "x$ac_cv_c_compiler_gnu" = "xyes"; then
    CFLAGS="$CFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wno-long-long -Wshadow -Wbad-function-cast -Wwrite-strings -Wstrict-prototypes -Wredundant-decls -Wnested-externs -Wundef -Wconversion -Wmissing-prototypes -Wmissing-declarations"
  fi
fi

# Manually add option for C99 support to compiler flags.
AX_CHECK_COMPILER_FLAGS([$ac_cv_prog_cc_c99],
  [CFLAGS="$CFLAGS $ac_cv_prog_cc_c99"],[AC_MSG_ERROR([Compiler does not accept option to accept C99.])])

# add Matlab CFLAGS
if test "x$ax_prog_matlab" = "xyes"; then
  CFLAGS="$CFLAGS $matlab_CFLAGS"
fi

# header files
AC_HEADER_STDC
AC_CHECK_HEADERS([stdbool.h libintl.h malloc.h stddef.h stdlib.h string.h \
  strings.h sys/time.h unistd.h limits.h c_asm.h intrinsics.h stdint.h \
  mach/mach_time.h sys/sysctl.h c_asm.h intrinsics.h sys/sysctl])
# c_asm.h: Header file for enabling asm() on Digital Unix
# intrinsics.h: cray unicos
# sys/sysctl.h: MacOS X altivec detection

# typedefs, structures, and compiler characteristics
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_CHECK_TYPE([long double],
  [AC_DEFINE(HAVE_LONG_DOUBLE, 1, [Define to 1 if the compiler supports `long double'])],[])

AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,
[
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
])

AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(unsigned long long)
AC_CHECK_SIZEOF(size_t)
AC_CHECK_SIZEOF(ptrdiff_t)

AC_CHECK_TYPES(uintptr_t, [], [AC_CHECK_SIZEOF(void *)], [$ac_includes_default
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#endif])

AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long double)

# library functions
AC_FUNC_ALLOCA
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_LIB(m, sin)

AC_CHECK_FUNCS([BSDgettimeofday gettimeofday gethrtime read_real_time])
AC_CHECK_FUNCS([time_base_to_time drand48 memset posix_memalign memalign])
AC_CHECK_FUNCS([_mm_malloc _mm_free clock_gettime mach_absolute_time sysctl])
AC_CHECK_FUNCS([abort snprintf sqrt tgamma lgamma log sin cos acos pow exp])
AC_CHECK_FUNCS([log1p sqrt tgamma lgamma log sin cos acos pow exp copysign lrint])
AC_CHECK_FUNCS([log1pf sqrtf tgammaf lgammaf logf sinf cosf acosf powf expf copysignf lrintf])
AC_CHECK_FUNCS([log1pl sqrtl tgammal lgammal logl sinl cosl acosl powl expl copysignl lrintl])
AC_CHECK_FUNCS([nanosleep sleep drand48])

AC_CHECK_DECLS([memalign, posix_memalign])
AC_CHECK_DECLS([log1p, sqrt, tgamma, lgamma, log, sin, cos, acos, pow, exp, copysign, lrint],[],[],[#include <math.h>])
AC_CHECK_DECLS([log1pf, sqrtf, tgammaf, lgammaf, logf, sinf, cosf, acosf, powf, expf, copysignf, lrintf],[],[],[#include <math.h>])
AC_CHECK_DECLS([log1pl, sqrtl, tgammal, lgammal, logl, sinl, cosl, acosl, powl, expl, copysignl, lrintl],[],[],[#include <math.h>])
AC_CHECK_DECLS([nanosleep],[],[],[#include <time.h>])
AC_CHECK_DECLS([drand48],[],[],[#include <stdlib.h>])
AC_CHECK_DECLS([sleep],[],[],[#include <unistd.h>])

# The output files to be generated
AC_OUTPUT(Makefile \
          include/Makefile \
          include/nfftconf.h \
          3rdparty/Makefile \
          3rdparty/cstripack/Makefile \
          util/Makefile \
          kernel/Makefile \
          kernel/fpt/Makefile \
          kernel/mri/Makefile \
          kernel/nfct/Makefile \
          kernel/nfft/Makefile \
          kernel/nfsft/Makefile \
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
          examples/nfst/Makefile \
          examples/nnfft/Makefile \
          examples/nnfft/accuracy/Makefile \
          examples/nnfft/simple_test/Makefile \
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
          matlab/nfsft/Makefile
          matlab/nfsft/@f_hat/Makefile
)

if test "x$USE_MAINTAINER_MODE" = "xyes"; then
AC_OUTPUT(tests/Makefile \
          VERSION
          doxygen/doxygen.Doxyfile
          tests/kernel/Makefile \
          tests/kernel/fpt/Makefile \
          tests/lowlevel/Makefile \
          tests/lowlevel/nfsft/Makefile \
          applications/iterS2/Makefile \
          applications/texture/Makefile \
          applications/texture/examplesTexture/Makefile \
          applications/texture/examplesTexture/correctness/Makefile \
          applications/texture/examplesTexture/data/Makefile \
          applications/texture/examplesTexture/data/data_publication/Makefile \
          applications/texture/examplesTexture/data/util/Makefile \
          applications/texture/examplesTexture/data_creator/Makefile \
          applications/texture/examplesTexture/interpolation/Makefile \
          applications/texture/examplesTexture/matrix/Makefile \
          applications/texture/examplesTexture/paper1/Makefile \
          applications/texture/examplesTexture/paper1/data/Makefile \
          applications/texture/examplesTexture/paper1/reconstruction/Makefile \
          applications/texture/examplesTexture/paper1/texture/Makefile \
          applications/texture/examplesTexture/performance/Makefile \
          applications/texture/examplesTexture/regularity/Makefile \
          applications/texture/examplesTexture/regularity2/Makefile \
          applications/texture/examplesTexture/regularity3/Makefile \
          applications/texture/examplesTexture/texture_util/Makefile
)
fi
