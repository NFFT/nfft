dnl $Id: nfsftmex.c 1548 2007-04-18 07:34:23Z keiner $
dnl
dnl Copyright (c) 2007 Jens Keiner
dnl
dnl This program is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by the Free
dnl Software Foundation; either version 2 of the License, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful, but WITHOUT
dnl ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
dnl FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
dnl more details.
dnl
dnl You should have received a copy of the GNU General Public License along with
dnl this program; if not, write to the Free Software Foundation, Inc., 51
dnl Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
dnl
dnl @synopsis AX_PROG_MATLAB
dnl @summary set up compiler and linker flags for Matlab mex file compilation.
dnl @category C
dnl
dnl @version 2007-04-24
dnl @license GPLWithACException
dnl @author Jens Keiner
dnl ############################################################################
dnl Process this file with autoconf to produce a configure script.
dnl ############################################################################

dnl Initialize Autoconf.
AC_INIT(NFFT, 3.0.3, potts@mathematik.tu-chemnitz.de)

dnl Add Copyright notice.
AC_COPYRIGHT([2003, 2007, Jens Keiner, Stefan Kunis, Daniel Potts])

dnl Add revision number of configure.ac.
AC_REVISION([$Id$])

dnl Set configure source directory.
AC_CONFIG_SRCDIR(include/nfft3.h)

dnl Set directory for auxilliary files.
AC_CONFIG_AUX_DIR(config)

dnl Define config headers.
AC_CONFIG_HEADERS([include/config.h:include/config.h.in])

dnl Check if make needs MAKE environment variable to be set.
AC_PROG_MAKE_SET

dnl Initialise automake.
AM_INIT_AUTOMAKE

dnl Determine host system type.
AC_CANONICAL_HOST

dnl "To hell with shared libraries"
AM_ENABLE_SHARED(no)

dnl Make abs_srcdir available in Makefiles.
AC_SUBST([abs_srcdir])

dnl ############################################################################
dnl Configure options
dnl ############################################################################

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile nfft with extra runtime checks for debugging])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(NFFT_DEBUG,1,[Define to enable extra NFFT debugging code.])
fi

dnl Make maintainer mode available.
AM_MAINTAINER_MODE

dnl ############################################################################
dnl Checks for programs.
dnl ############################################################################

dnl Set programming language to C.
AC_LANG(C)

dnl Search for C compiler. Automatically selects options to support C89.
AC_PROG_CC

dnl Determine compiler vendor (thanks fftw)
AX_COMPILER_VENDOR

dnl Check if the compiler accepts C99.
AC_PROG_CC_STDC

if test $ac_cv_prog_cc_stdc = "no"; then
  AC_MSG_ERROR([Your compiler does not accept C99!])
fi

dnl Check for assembler (thanks fftw).
AMX_PROG_AS dnl AM_PROG_AS breaks at least one version of libtool.

dnl Search for the C preprocessor.
AC_PROG_CPP_WERROR

dnl Check for BSD-compatible install program. Enforces inclusion of own install
dnl script.
AC_PROG_INSTALL

dnl Check if ln -s works.
AC_PROG_LN_S

dnl Make libtool ddl capabilities available.
AC_LIBTOOL_WIN32_DLL

dnl Check for libtool.
AC_PROG_LIBTOOL

dnl Check for MATLAB.
AX_PROG_MATLAB

dnl Determine (hopefully) optimal CFLAGS (thanks fftw)
AX_CC_MAXOPT

dnl Make CFLAGS available.
AC_SUBST(CFLAGS)

AX_LIB_GSL

# Checks for libraries
AC_CHECK_LIB(fftw3,fftw_execute,,,)
LIBS=

# Checks for library functions.

# Checks for typedefs, structures, and compiler characteristics



# The output files to be generated
AC_OUTPUT(Makefile \
          include/Makefile \
          3rdparty/Makefile \
          3rdparty/gsl/Makefile \
          3rdparty/gsl/err/Makefile \
          3rdparty/gsl/sys/Makefile \
          3rdparty/gsl/specfunc/Makefile \
          3rdparty/gsl/ieee-utils/Makefile \
          3rdparty/gsl/utils/Makefile \
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
	  applications/mri/Makefile \
          applications/mri/mri2d/Makefile \
          applications/mri/mri3d/Makefile \
	  applications/polarFFT/Makefile \
	  applications/radon/Makefile \
	  matlab/Makefile \
	  matlab/nfsft/Makefile
)
#	  applications/quadratureS2/Makefile \

if test $USE_MAINTAINER_MODE = yes; then
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
					applications/texture/examplesTexture/performance/Makefile \
					applications/texture/examplesTexture/regularity/Makefile \
					applications/texture/examplesTexture/regularity2/Makefile \
					applications/texture/examplesTexture/regularity3/Makefile \
					applications/texture/examplesTexture/texture_util/Makefile
)
fi

