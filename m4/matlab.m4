# $Id$
#
# Copyright (c) 2007 Jens Keiner
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
# @synopsis AX_PROG_MATLAB
# @summary set up compiler and linker flags for Matlab mex file compilation.
# @category C
#
# @version 2008-11-25
# @license GPLWithACException
# @author Jens Keiner

AC_DEFUN([AX_PROG_MATLAB],
[
  AC_REQUIRE([AC_CANONICAL_HOST])

  # option to enable mex file compilation
  AC_ARG_WITH(matlab,
  [AC_HELP_STRING([--with-matlab=DIR],[the directory where Matlab is installed])],
  matlab_dir=${withval},matlab_dir=)

  AC_MSG_CHECKING([whether to check for Matlab])

  # Check whether Matlab directory exists.
  if test -d "${matlab_dir}"; then
    # directory existent
    AC_MSG_RESULT([yes])

    matlab_include_dir="${matlab_dir}/extern/include"
    matlab_src_dir="${matlab_dir}/extern/src"
    matlab_bin_dir="${matlab_dir}/bin"
    matlab_CPPFLAGS="-I${matlab_include_dir} -I${matlab_src_dir} -DMATLAB_MEX_FILE"

    # Save environment
    saved_CPPFLAGS=$CPPFLAGS

    # augmented environment
    CPPFLAGS="$CPPFLAGS $matlab_CPPFLAGS"

    # header files
    AC_CHECK_FILE([${matlab_include_dir}/mat.h],[],AC_MSG_ERROR([Required file ]${matlab_include_dir}[/mat.h not found]))
    AC_CHECK_FILE([${matlab_include_dir}/matrix.h],[],AC_MSG_ERROR([Required file ]${matlab_include_dir}[/matrix.h not found]))
    AC_CHECK_FILE([${matlab_include_dir}/mex.h],[],AC_MSG_ERROR([Required file ]${matlab_include_dir}[/mex.h not found]))
    AC_CHECK_HEADER(mat.h,[],AC_MSG_ERROR([Required header mat.h not usable]))
    AC_CHECK_HEADER(matrix.h,[],AC_MSG_ERROR([Required header matrix.h not usable]))
    AC_CHECK_HEADER(mex.h,[],AC_MSG_ERROR([Required header mex.h not usable]))

    # host specific stuff
    case $host in
      powerpc-*darwin*) # Mac (PowerPC)
        matlab_check_mexversion_c="yes"
        matlab_dir_prefix="mac"
        matlab_mexext=".mexmac"
        matlab_libext=".dylib";;
      i686-*darwin*) # Mac (Intel)
        matlab_check_mexversion_c="yes"
        matlab_dir_prefix="maci"
        matlab_mexext=".mexmaci"
        matlab_libext=".dylib";;
      *86_64*linux*) # Linux (x86, 64 bit)
        matlab_check_mexversion_c="yes"
        matlab_dir_prefix="glna64"
        matlab_mexext=".mexa64"
        matlab_libext=".so";;
      *86*linux*) # Linux (x86, 32 bit)
        matlab_check_mexversion_c="yes"
        matlab_dir_prefix="glnx86"
        matlab_mexext=".mexglx"
        matlab_libext=".so";;
      *irix*) # IRIX (mips)
        matlab_check_mexversion_c="yes"
        matlab_dir_prefix="sgi"
        matlab_mexext=".mexsg"
        matlab_libext=".so";;
      *solaris*) # FIXME Matlab on Solaris only supported on 64bit systems
        matlab_check_mexversion_c="yes"
        matlab_dir_prefix="sol"
        matlab_mexext=".mexsg"
        matlab_libext=".so";;
      *cygwin*|*mingw*)
        # Windows (Cygwin or MinGW)
        matlab_check_mexversion_c="no"
        matlab_dir_prefix="win32"
        if test ! -e "${matlab_dir}/bin/win32/libmx.a"; then
          cd ${matlab_bin_dir}/win32
          libmx=`dlltool -llibmx.a -d${matlab_include_dir}/libmx.def -Dlibmx.dll`
          cd -
        fi
        if test ! -e "${matlab_dir}/bin/win32/libmex.a"; then
          cd ${matlab_bin_dir}/win32
          libmex=`dlltool -llibmex.a -d${matlab_include_dir}/libmex.def -Dlibmex.dll`
          cd -
        fi
        if test ! -e "${matlab_dir}/bin/win32/libmat.a"; then
          cd ${matlab_bin_dir}/win32
          libmat=`dlltool -llibmat.a -d${matlab_include_dir}/libmat.def -Dlibmat.dll`
          cd -
        fi
        matlab_mexext=".dll"
        matlab_libext=".a";;
    esac

    matlab_host_bin_dir="${matlab_bin_dir}/${matlab_dir_prefix}"

    if test "$matlab_check_mexversion_c" = "yes"; then
      AC_CHECK_FILE([${matlab_src_dir}/mexversion.c],AC_DEFINE([HAVE_MEXVERSION_C],[1],[Define to have the file mexversion.c]),AC_MSG_ERROR([Required file ]${matlab_src_dir}[/mexversion.c not found]))
    fi

    AC_CHECK_FILE([${matlab_host_bin_dir}/libmat${matlab_libext}],[],AC_MSG_ERROR([Required library ]${matlab_host_bin_dir}[/libmat]${matlab_libext}[ not found]))
    AC_CHECK_FILE([${matlab_host_bin_dir}/libmex${matlab_libext}],[],AC_MSG_ERROR([Required library ]${matlab_host_bin_dir}[/libmex]${matlab_libext}[ not found]))
    AC_CHECK_FILE([${matlab_host_bin_dir}/libmx${matlab_libext}],[],AC_MSG_ERROR([Required library ]${matlab_host_bin_dir}[/libmx]${matlab_libext}[ not found]))

    matlab_LIBADD="-lmx -lmex -lmat"
    matlab_LDFLAGS="-L${matlab_host_bin_dir}"

    AC_SUBST([matlab_dir])
    AC_SUBST([matlab_CPPFLAGS])
    AC_SUBST([matlab_LDFLAGS])
    AC_SUBST([matlab_LIBADD])
    AC_SUBST([matlab_mexext])
    AM_CONDITIONAL(HAVE_MATLAB, test "xyes" = "xyes" )

    # Restore environment.
    CPPFLAGS=$saved_CPPFLAGS
  else
    AC_MSG_RESULT([no])
    AM_CONDITIONAL(HAVE_MATLAB, test "xno" = "xyes" )
  fi
])
