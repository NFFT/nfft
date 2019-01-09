# Copyright (c) 2002, 2019 Jens Keiner, Stefan Kunis, Daniel Potts
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
################################################################################
# Process this file with autoconf to produce a configure script.
################################################################################

m4_define([nfft_version_major], [3])
m4_define([nfft_version_minor], [5])
m4_define([nfft_version_patch], [0])
m4_define([nfft_version_type], [alpha])
m4_append([NFFT_VERSION], m4_expand([nfft_version_major.nfft_version_minor.nfft_version_patch]))
m4_append([NFFT_VERSION], m4_expand([nfft_version_type]))
...
# autoconf initialization
AC_INIT([NFFT],[NFFT_VERSION],[mail@nfft.org])

# copyright notice
AC_COPYRIGHT([2003, 2019, Jens Keiner, Stefan Kunis, Daniel Potts])

AC_DEFINE_UNQUOTED([NFFT_VERSION_MAJOR], [nfft_version_major], [Major version number.])
AC_DEFINE_UNQUOTED([NFFT_VERSION_MINOR], [nfft_version_minor], [Minor version number.])
AC_DEFINE_UNQUOTED([NFFT_VERSION_PATCH], [nfft_version_patch], [Patch version number.])
AC_DEFINE_UNQUOTED([NFFT_VERSION_TYPE], [nfft_version_type], [Release type.])

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

# number of CPUs
#AX_COUNT_CPUS
#AC_DEFINE_UNQUOTED(NFFT_NUM_CORES,$cpu_count,[Define to number of cores.])

# substitute abs_srcdir in generated Makefiles.
AC_SUBST([abs_srcdir])
if test "x${host_os}" = "xmingw32" -o "x${host_os}" = "xmingw64"; then
  ABS_SRCDIR=`cd "$srcdir" && pwd -W`
else
  ABS_SRCDIR=`cd "$srcdir" && pwd`
fi
AC_DEFINE_UNQUOTED([ABS_SRCDIR], ["$ABS_SRCDIR"], [Absolute path to source directory.])
AC_SUBST([ABS_SRCDIR])

AC_SUBST([abs_top_srcdir])
AC_SUBST([top_builddir])
AC_SUBST([abs_top_builddir])

# minimum required libtool version
LT_PREREQ([2.4.2])

# libtool initialization
LT_INIT([win32-dll])

# substitute LIBTOOL_DEPS variable in generated Makefiles
AC_SUBST([LIBTOOL_DEPS])

# version information for shared library
SHARED_VERSION_INFO="3:0:1"

# substitute SHARED_VERSION_INFO in generated Makefiles
AC_SUBST(SHARED_VERSION_INFO)

# check for make
AC_PROG_MAKE_SET

# automake initialization
AM_INIT_AUTOMAKE(1.10)

# doxygen ${PACKAGE_STRING}
DX_DOT_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([NFFT], [doxygen/doxygen.Doxyfile], [doc])

################################################################################
# options for customizing the build process
################################################################################

# maintainer mode option
AM_MAINTAINER_MODE

# Precision.
AC_ARG_ENABLE(single, [AS_HELP_STRING([--enable-single],[compile NFFT in single precision])], ok=$enableval, ok=no)
AC_ARG_ENABLE(float, [AS_HELP_STRING([--enable-float],[synonym for --enable-single])], ok=$enableval)
if test "$ok" = "yes"; then
  AC_DEFINE(NFFT_SINGLE,1,[Define to compile in single precision.])
  PRECISION=s
else
  PRECISION=d
fi

AC_ARG_ENABLE(long-double, [AS_HELP_STRING([--enable-long-double],[compile NFFT in long-double precision])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  if test "$PRECISION" = "s"; then
    AC_MSG_ERROR([--enable-single/--enable-long-double conflict])
  fi
  AC_DEFINE(NFFT_LDOUBLE,1,[Define to compile in long-double precision.])
  PRECISION=l
fi
AC_SUBST(PRECISION)

# Library suffix
case "$PRECISION" in
     s) PREC_SUFFIX=f;;
     d) PREC_SUFFIX=;;
     l) PREC_SUFFIX=l;;
esac
AC_SUBST(PREC_SUFFIX)

# Library suffix
case "$PRECISION" in
     s) NFFT_PRECISION_MACRO=NFFT_PRECISION_SINGLE;;
     d) NFFT_PRECISION_MACRO=NFFT_PRECISION_DOUBLE;;
     l) NFFT_PRECISION_MACRO=NFFT_PRECISION_LONG_DOUBLE;;
esac
AC_SUBST(NFFT_PRECISION_MACRO)

AM_CONDITIONAL(HAVE_NON_DOUBLE_PRECISION, test "x$PRECISION" != "xd" )

# enable or disable parts of NFFT

# whether we need the fpt module
need_fpt="no"

# build all modules by default in maintainer mode or if option is given
AC_ARG_ENABLE(all, [AC_HELP_STRING([--enable-all],[build all modules])],
  ok=$enableval, ok=$USE_MAINTAINER_MODE)
nfft_module_default=$ok

# option for example programs
nfft_examples_default="yes"
AC_ARG_ENABLE(examples, [AC_HELP_STRING([--enable-examples],
  [enable example programs])], enable_examples=$enableval, enable_examples=$nfft_examples_default)
AM_CONDITIONAL(HAVE_EXAMPLES, test "x$enable_examples" = "xyes" )

# option for application programs
nfft_applications_default="yes"
AC_ARG_ENABLE(applications, [AC_HELP_STRING([--enable-applications],
  [enable application programs])], enable_applications=$enableval, enable_applications=$nfft_applications_default)
AM_CONDITIONAL(HAVE_APPLICATIONS, test "x$enable_applications" = "xyes" )

# option for julia interface
nfft_julia_default="no"
AC_ARG_ENABLE(julia, [AC_HELP_STRING([--enable-julia],
  [enable julia interface])], enable_julia=$enableval, enable_julia=$nfft_julia_default)
AM_CONDITIONAL(HAVE_JULIA, test "x$enable_julia" = "xyes" )

if test "x$enable_julia" = "xyes" -a "x$PRECISION" != "xd"; then
  AC_MSG_ERROR([Building the julia interface requires double precision.])
fi


# options for modules
AX_NFFT_MODULE([nfct],[NFCT],[nonequispaced fast cosine transform],["yes"])
AX_NFFT_MODULE([nfst],[NFST],[nonequispaced fast sine transform],["yes"])
AX_NFFT_MODULE([nfsft],[NFSFT],[nonequispaced fast spherical Fourier transform],["no"],
  [need_fpt="yes"])
AX_NFFT_MODULE([nfsoft],[NFSOFT],[nonequispaced fast SO(3) Fourier transform],["no"],
  [need_fpt="yes"])
AX_NFFT_MODULE([nnfft],[NNFFT],[nonequispaced fast Fourier transform -- ]#
  [nonequispaced in both domains],["no"])
AX_NFFT_MODULE([nsfft],[NSFFT],[nonequispaced sparse fast Fourier transform],["no"])
AX_NFFT_MODULE([mri],[MRI],[magnet resonance imaging],["no"])
AX_NFFT_MODULE([fpt],[FPT],[fast polynomial transform],["no"],[],[],
  [test "x$ok" = "xyes" -o "x$need_fpt" = "xyes"])

# multithreaded code
#AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],
#  [enable OpenMP multithreaded code])], [enable_threads=$enableval; AC_DEFINE(ENABLE_OPENMP, 1, ["Define to enable OpenMP code."])], enable_threads=no)
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],
  [enable OpenMP multithreaded code])], enable_threads=$enableval, enable_threads=no)
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

AC_ARG_ENABLE(mips_zbus_timer, [AC_HELP_STRING([--enable-mips-zbus-timer],
  [use MIPS ZBus cycle-counter])], have_mips_zbus_timer=$enableval,
  have_mips_zbus_timer=no)
if test "$have_mips_zbus_timer" = "yes"; then
  AC_DEFINE(HAVE_MIPS_ZBUS_TIMER,1,
    [Define to enable use of MIPS ZBus cycle-counter.])
fi

# select window function
AC_ARG_WITH(window, [AC_HELP_STRING([--with-window=ARG],[choose window function
  (ARG can be one of: kaiserbessel (default), gaussian, bspline, sinc,
  dirac)])],
  window=$withval, window="kaiserbessel")

AC_MSG_CHECKING([window function])
case "$window" in
  kaiserbessel)
    AC_DEFINE(KAISER_BESSEL,1,[Define to enable Kaiser-Bessel window function.]);;
  gaussian)
    AC_DEFINE(GAUSSIAN,1,[Define to enable Gaussian window function.]);;
  bspline)
    AC_DEFINE(B_SPLINE,1,[Define to enable B-spline window function.]);;
  sinc)
    AC_DEFINE(SINC_POWER,1,[Define to enable sinc power window function.]);;
  delta)
    AC_DEFINE(DIRAC_DELTA,1,[Define to enable Dirac delta window function.]);;
  *)
    AC_MSG_ERROR([Unknown window function "$window".]);;
esac
AC_MSG_RESULT([$window])
AC_DEFINE_UNQUOTED(WINDOW_NAME,$window,[The name of the window function.])

AC_ARG_ENABLE(exhaustive-unit-tests, [AS_HELP_STRING([--enable-exhaustive-unit-tests],[enable exhaustive unit tests])], ok=$enableval, ok=no)
matlab_exhaustive_unit_tests_flag=0
if test "$ok" = "yes"; then
  AC_DEFINE(NFFT_EXHAUSTIVE_UNIT_TESTS,1,[Define to enable exhaustive unit tests.])
  matlab_exhaustive_unit_tests_flag=1
fi
AC_SUBST(matlab_exhaustive_unit_tests_flag)

################################################################################
# compiler characteristis
################################################################################

# select programming language
AC_LANG(C)

# compiler vendor
AX_COMPILER_VENDOR

# check for C99 compliant mode (possibly with GNU extensions)
AC_PROG_CC_C99

# per-target flags
AM_PROG_CC_C_O

# enable "const" keyword
AC_C_CONST

# enable "restrict" keyword
AC_C_RESTRICT

# enable "inline" keyword
AC_C_INLINE

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

# mkdir
AC_PROG_MKDIR_P

# whether ln -s works
AC_PROG_LN_S

################################################################################
# 3rd party libraries
################################################################################

# Check for fftw3.
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

AX_OPENMP
AM_CONDITIONAL(HAVE_OPENMP, test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xunknown" )
AC_SUBST(OPENMP_CFLAGS)

if  test "x$enable_threads" = "xyes" -a "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" = "xunknown"; then
  AC_MSG_ERROR([You don't seem to have a C compiler with OpenMP support installed which is required for threaded NFFT.])
fi

# Check for MATLAB.
AX_PROG_MATLAB

if test "x$ax_prog_matlab" = "xyes" -a "x$PRECISION" != "xd"; then
  AC_MSG_ERROR([Building the Matlab interfaces requires double precision.])
fi

if test "x$matlab_threads" = "xyes" -a "x$enable_threads" != "xyes"; then
  AC_MSG_ERROR([The NFFT Matlab interface with thread support requires the threaded NFFT to be built. Please re-run configure with \"--enable-openmp\".])
fi


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
#  if test "x$ax_cv_c_compiler_vendor" = "xapple"; then
#    CFLAGS="$CFLAGS $ax_cv_apple_gcc_archflag"
#  fi
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

# use MinGW implementation of printf
if test "x${host_os}" = "xmingw32" -o "x${host_os}" = "xmingw64"; then
  CFLAGS="$CFLAGS -D__USE_MINGW_ANSI_STDIO=1"
fi  

CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"

# add Matlab CFLAGS
if test "x$ax_prog_matlab" = "xyes"; then
  CFLAGS="$CFLAGS $matlab_CFLAGS"
fi

################################################################################
# check if mex.h can be included without undefinition of
################################################################################
ax_prog_matlab_gcc_require_undef_stdc_utf_16=no
if test "x$ax_prog_matlab" = "xyes"; then
  AC_MSG_CHECKING([if inclusion of mex.h succeeds])
  saved_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $matlab_CPPFLAGS"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include "mex.h"]],[[]])],[ax_prog_matlab_check_compile_mex_h=yes], [ax_prog_matlab_check_compile_mex_h=no])
  AC_MSG_RESULT([$ax_prog_matlab_check_compile_mex_h])
  if test "x$ax_prog_matlab_check_compile_mex_h" = "xno"; then
    AC_MSG_CHECKING([if undefining __STDC_UTF_16__ solves this issue])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#undef __STDC_UTF_16__
    #include "mex.h"]],[[]])],[ax_prog_matlab_gcc_require_undef_stdc_utf_16=yes], [])
    AC_MSG_RESULT([$ax_prog_matlab_gcc_require_undef_stdc_utf_16])
  fi
  CPPFLAGS="$saved_CPPFLAGS"
fi
if test "x$ax_prog_matlab_gcc_require_undef_stdc_utf_16" =  "xyes"; then
  AC_DEFINE([HAVE_MATLAB_GCC_REQUIRE_UNDEF_STDC_UTF_16],[1],[Define to 1 if compilation of Matlab interface only succeeds when macro __STDC_UTF_16__ not defined.])
fi

################################################################################
# header files/data types/compiler characteristics
################################################################################

AC_CHECK_HEADERS([math.h stdio.h stdlib.h time.h  sys/time.h \
  complex.h string.h float.h limits.h stdarg.h stddef.h sys/types.h stdint.h \
  inttypes.h stdbool.h malloc.h c_asm.h intrinsics.h mach/mach_time.h])

AC_HEADER_TIME

AC_TYPE_SIZE_T
AC_CHECK_TYPE([long double],
  [AC_DEFINE(HAVE_LONG_DOUBLE, 1, [Define to 1 if the compiler supports]
  ['long double'])],[])
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
AC_CHECK_LIB(m, sin, [])
LIBS=$saved_LIBS

AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])

AC_CHECK_FUNCS([memset posix_memalign memalign])
AC_CHECK_FUNCS([_mm_malloc _mm_free sysctl])
AC_CHECK_FUNCS([abort snprintf sqrt])
AC_CHECK_FUNCS([sleep usleep nanosleep drand48 srand48])
AC_CHECK_FUNCS([gethostname])

AC_CHECK_DECLS([memalign, posix_memalign])
AC_CHECK_DECLS([sleep],[],[],[#include <unistd.h>])
AC_CHECK_DECLS([nanosleep],[],[],[#include <time.h>])
AC_CHECK_DECLS([drand48],[],[],[#include <stdlib.h>])
AC_CHECK_DECLS([srand48],[],[],[#include <stdlib.h>])

# Cray UNICOS _rtc() (real-time clock) intrinsic
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_TRY_LINK([#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif], [_rtc()], [AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])], [rtc_ok=no])
AC_MSG_RESULT($rtc_ok)

AC_MSG_CHECKING([whether a cycle counter is available])
save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$srcdir/include"
AC_TRY_CPP([#include "cycle.h"
#ifndef HAVE_TICK_COUNTER
#  error No cycle counter
#endif], [ok=yes], [ok=no])
AC_MSG_RESULT($ok)
TICKS_PER_SECOND=1
if test $ok = no; then
  echo "***************************************************************"
  echo "WARNING: No cycle counter found. Time measurements in NFFT will"
  echo "         show incorrect results.                               "
  echo "***************************************************************"
else
  AC_TRY_RUN([#include "cycle.h"
#include <stdio.h>
#if defined(HAVE_NANOSLEEP)
#include <stdlib.h>
static struct timespec ts = {1,0};
#define SLEEP nanosleep(&ts, 0)
#elif defined(HAVE_SLEEP)
#include <unistd.h>
#define SLEEP sleep(1)
#elif defined(HAVE_USLEEP)
#include <unistd.h>
#define SLEEP usleep(1000)
#else
#error wtf!
#endif
int main(){
ticks t0 = getticks();
SLEEP;{
  double tps = elapsed(getticks(),t0)/1.0;
  FILE *f = fopen("ticks.tmp","w");
  fprintf(f,"%.1f\n",tps);
  fclose(f);}
return 0;}],[read TICKS_PER_SECOND < ticks.tmp],[],[TICKS_PER_SECOND=1])
    rm -f ticks.tmp
  if test "$TICKS_PER_SECOND" = "1"; then
    echo "***************************************************************"
    echo "WARNING: Number of ticks per second could not be determined.   "
    echo "         Time measurements in NFFT will be in arbitrary units  "
    echo "         instead of seconds.                                   "
    echo "***************************************************************"
  fi
fi

CPPFLAGS=$save_CPPFLAGS
AC_SUBST(TICKS_PER_SECOND)

if test "$PRECISION" = "l"; then
#AC_CHECK_FUNCS([copysignl nextafterl nanl ceill floorl nearbyintl rintl roundl lrintl lroundl llrintl llroundl truncl fmodl remainderl remquol fdiml fmaxl fminl fmal fabsl sqrtl cbrtl hypotl expl exp2l expm1l logl log2l log10l log1pl logbl ilogbl modfl frexpl ldexpl scalbnl scalblnl powl cosl sinl tanl coshl sinhl tanhl acosl asinl atanl atan2l acoshl asinhl atanhl tgammal lgammal j0l j1l jnl jnl y0l y1l ynl erfl erfcl creall cimagl cabsl cargl conjl cprojl csqrtl cexpl clogl cpowl csinl ccosl ctanl casinl cacosl catanl csinhl ccoshl ctanhl casinhl cacoshl catanhl])
AC_CHECK_DECLS([copysignl, nextafterl, nanl, ceill, floorl, nearbyintl, rintl, roundl, lrintl, lroundl, llrintl, llroundl, truncl, fmodl, remainderl, remquol, fdiml, fmaxl, fminl, fmal, fabsl, sqrtl, cbrtl, hypotl, expl, exp2l, expm1l, logl, log2l, log10l, log1pl, logbl, ilogbl, modfl, frexpl, ldexpl, scalbnl, scalblnl, powl, cosl, sinl, tanl, coshl, sinhl, tanhl, acosl, asinl, atanl, atan2l, acoshl, asinhl, atanhl, tgammal, lgammal, j0l, j1l, jnl, y0l, y1l, ynl, erfl, erfcl, creall, cimagl, cabsl, cargl, conjl, cprojl, csqrtl, cexpl, clogl, cpowl, csinl, ccosl, ctanl, casinl, cacosl, catanl, csinhl, ccoshl, ctanhl, casinhl, cacoshl, catanhl],[],[],
  [#include <math.h>
#include <complex.h>])
fi

if test "$PRECISION" = "d"; then
#AC_CHECK_FUNCS([copysign nextafter nan ceil floor nearbyint rint round lrint lround llrint llround trunc fmod remainder remquo fdim fmax fmin fma fabs sqrt cbrt hypot exp exp2 expm1 log log2 log10 log1p logb ilogb modf frexp ldexp scalbn scalbln pow cos sin tan cosh sinh tanh acos asin atan atan2 acosh asinh atanh tgamma lgamma j0 j1 jn y0 y1 yn erf erfc creal cimag cabs carg conj cproj csqrt cexp clog cpow csin ccos ctan casin cacos catan csinh ccosh ctanh casinh cacosh catanh])
AC_CHECK_DECLS([copysign, nextafter, nan, ceil, floor, nearbyint, rint, round, lrint, lround, llrint, llround, trunc, fmod, remainder, remquo, fdim, fmax, fmin, fma, fabs, sqrt, cbrt, hypot, exp, exp2, expm1, log, log2, log10, log1p, logb, ilogb, modf, frexp, ldexp, scalbn, scalbln, pow, cos, sin, tan, cosh, sinh, tanh, acos, asin, atan, atan2, acosh, asinh, atanh, tgamma, lgamma, j0, j1, jn, y0, y1, yn, erf, erfc, creal, cimag, cabs, carg, conj, cproj, csqrt, cexp, clog, cpow, csin, ccos, ctan, casin, cacos, catan, csinh, ccosh, ctanh, casinh, cacosh, catanh],[],[],
  [#include <math.h>
#include <complex.h>])
fi

if test "$PRECISION" = "s"; then
#AC_CHECK_FUNCS([copysignf nextafterf nanf ceilf floorf nearbyintf rintf roundf lrintf lroundf llrintf llroundf truncf fmodf remainderf remquof fdimf fmaxf fminf fmaf fabsf sqrtf cbrtf hypotf expf exp2f expm1f logf log2f log10f log1pf logbf ilogbf modff frexpf ldexpf scalbnf scalblnf powf cosf sinf tanf coshf sinhf tanhf acosf asinf atanf atan2f acoshf asinhf atanhf tgammaf lgammaf j0f j1f jnf y0f y1f ynf erff erfcf crealf cimagf cabsf cargf conjf cprojf csqrtf cexpf clogf cpowf csinf ccosf ctanf casinf cacosf catanf csinhf ccoshf ctanhf casinhf cacoshf catanhf])
AC_CHECK_DECLS([copysignf, nextafterf, nanf, ceilf, floorf, nearbyintf, rintf, roundf, lrintf, lroundf, llrintf, llroundf, truncf, fmodf, remainderf, remquof, fdimf, fmaxf, fminf, fmaf, fabsf, sqrtf, cbrtf, hypotf, expf, exp2f, expm1f, logf, log2f, log10f, log1pf, logbf, ilogbf, modff, frexpf, ldexpf, scalbnf, scalblnf, powf, cosf, sinf, tanf, coshf, sinhf, tanhf, acosf, asinf, atanf, atan2f, acoshf, asinhf, atanhf, tgammaf, lgammaf, j0f, j1f, jnf, y0f, y1f, ynf, erff, erfcf, crealf, cimagf, cabsf, cargf, conjf, cprojf, csqrtf, cexpf, clogf, cpowf, csinf, ccosf, ctanf, casinf, cacosf, catanf, csinhf, ccoshf, ctanhf, casinhf, cacoshf, catanhf],[],[],
  [#include <math.h>
#include <complex.h>])
fi

# Get machine epsilon by 2^(1-MAND_DIG) for selected precision per float.h.
AX_EPS_DEF([$PRECISION])

if test "$ax_cv_eps_def" = "unknown"; then  
  AC_MSG_ERROR([Unable to determine floating-point epsilon.])
else
  AC_DEFINE_UNQUOTED([NFFT_EPSILON], [$ax_cv_eps_def], [Floating-point epsilon.])
fi

# CUnit
AX_CUNIT
AM_CONDITIONAL(HAVE_CUNIT, test "x$ax_cv_cunit" = "xyes" )

# The output files to be generated
AC_CONFIG_FILES(Makefile \
                nfft3.pc \
                doxygen/doxygen.Doxyfile \
                include/Makefile \
                include/ticks.h \
                3rdparty/Makefile \
                3rdparty/cstripack/Makefile \
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
                kernel/util/Makefile \
                tests/Makefile \
                tests/data/Makefile \
                examples/Makefile \
                examples/fpt/Makefile \
                examples/mri/Makefile \
                examples/nfct/Makefile \
                examples/nfct/simple_test.c \
                examples/nfft/Makefile \
                examples/nfft/simple_test.c \
                examples/nfft/simple_test_threads.c \
                examples/nfsft/Makefile \
                examples/nfsoft/Makefile \
                examples/nfst/Makefile \
                examples/nfst/simple_test.c \
                examples/nnfft/Makefile \
                examples/nsfft/Makefile \
                examples/solver/Makefile \
		examples/solver/glacier.c \
		examples/solver/simple_test.c \
                applications/Makefile \
                applications/fastgauss/Makefile \
		applications/fastgauss/fastgauss.c \
                applications/fastsum/Makefile \
                applications/fastsumS2/Makefile \
                applications/quadratureS2/Makefile \
                applications/mri/Makefile \
                applications/mri/mri2d/Makefile \
                applications/mri/mri3d/Makefile \
                applications/polarFFT/Makefile \
                applications/polarFFT/linogram_fft_test.c \
                applications/polarFFT/mpolar_fft_test.c \
                applications/polarFFT/polar_fft_test.c \
                applications/radon/Makefile \
		applications/radon/inverse_radon.c \
		applications/radon/radon.c \
                matlab/Makefile \
                matlab/fpt/Makefile \
                matlab/nfsft/Makefile \
                matlab/nfft/Makefile \
                matlab/nfct/Makefile \
                matlab/nfst/Makefile \
                matlab/nnfft/Makefile \
                matlab/fastsum/Makefile \
                matlab/infft1d/Makefile \
                matlab/nfsft/@f_hat/Makefile \
                matlab/nfsoft/Makefile \
                matlab/tests/Makefile \
		matlab/tests/check_nfft_matlab.sh \
		matlab/tests/check_nfft_octave.sh \
		matlab/tests/check_nfsft_matlab.sh \
		matlab/tests/check_nfsft_octave.sh \
		matlab/tests/check_nfsoft_matlab.sh \
		matlab/tests/check_nfsoft_octave.sh \
		matlab/tests/nfftUnitTests.m \
		matlab/tests/nfsftUnitTests.m \
		julia/Makefile \
		julia/nfft/Makefile \
                doxygen/Makefile \
                support/Makefile)
# temproarily removed:
#  applications/iterS2/Makefile \

AC_OUTPUT
