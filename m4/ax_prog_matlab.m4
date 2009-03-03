# $Id$
# 
# Copyright (c) 2007, 2008 Jens Keiner
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
# @version 2008-12-07
# @license GPLWithACException
# @author Jens Keiner <keiner@math.uni-luebeck.de>

AC_DEFUN([AX_PROG_MATLAB],
[
  AC_REQUIRE([AC_CANONICAL_HOST])

  # option to enable mex file compilation
  AC_ARG_WITH(matlab,
    [AC_HELP_STRING([--with-matlab=DIR],
      [the directory where Matlab is installed])],
    matlab_dir=${withval},matlab_dir="no")

  AC_ARG_WITH(matlab-arch,
    [AC_HELP_STRING([--with-matlab-arch=DIR],
      [Matlab architecture acronym])],
    matlab_arch=${withval},matlab_arch="yes")

  AC_ARG_ENABLE(matlab-argchecks,
    [AC_HELP_STRING([--enable-matlab-argchecks],
      [Compile Matlab interface with argument checks [default=yes]])],
      [ok="$enableval"],
      [ok="yes"])

  if test "x$ok" = "xyes"; then
    AC_DEFINE(MATLAB_ARGCHECKS,1,[Define to enable Matlab argument checks.])
  fi

  AC_MSG_CHECKING([whether to check for Matlab])

  if test "x${matlab_dir}" = "xno"; then
    AC_MSG_RESULT([no])
    ax_prog_matlab="no"
  else
    AC_MSG_RESULT([yes])

    # Matlab root
    AX_CHECK_DIR([${matlab_dir}],[],
      [AC_MSG_ERROR([Please supply a valid path to a Matlab root directory or run configure without the option --with-matlab.])])

    # subdirectories
    matlab_bin_dir="${matlab_dir}/bin"
    matlab_src_dir="${matlab_dir}/extern/src"
    matlab_include_dir="${matlab_dir}/extern/include"
    AX_CHECK_DIR([${matlab_bin_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])
    AX_CHECK_DIR([${matlab_src_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])
    AX_CHECK_DIR([${matlab_include_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])

    # architecture and mex file extension
    if test ! "x${matlab_arch}" = "xyes"; then
      AC_MSG_CHECKING([for architecture])
      AC_MSG_RESULT([${matlab_arch}])
      # mex file extension for architecture
      AC_MSG_CHECKING([for mex file extension])
      case $matlab_arch in
        glnx86) matlab_mexext="mexglx";;
        glnxa64) matlab_mexext="mexa64";;
        mac) matlab_mexext="mexmac";;
        maci) matlab_mexext="mexmaci";;
        maci64) matlab_mexext="mexmaci64";;
        sol64) matlab_mexext="mexs64";;
        win32) matlab_mexext="mexw32";;
        win64) matlab_mexext="mexw64";;
        *) AC_MSG_ERROR([Unsupported or invalid architecture ${matlab_arch}.]);;
      esac
      AC_MSG_RESULT([${matlab_mexext}])
    else 
      matlab_mexext="unknown"
      matlab_arch="unknown"

      # mex file extension, maybe
      for matlab_check_prog_mexext in "mexext mexext.sh mexext.bat"; do
        AC_PATH_PROG([matlab_prog_mexext],[$matlab_check_prog_mexext],[no],
          [$PATH$PATH_SEPARATOR$matlab_bin_dir])
        if test ! "x${matlab_prog_mexext}" = "xno"; then
          AC_MSG_CHECKING([for mex file extension])
          matlab_mexext=`${matlab_prog_mexext}`
          AC_MSG_RESULT([${matlab_mexext}])
          break
        fi
      done

      # architecture, maybe
      if test "x${matlab_mexext}" = "xunknown"; then
        # Try guessing the architecture based on host CPU
        case $host in
          *86_64*linux*) matlab_arch_test="glnxa64";;
          *86*linux*) matlab_arch_test="glnx86";;
          *powerpc*darwin*) matlab_arch_test="mac mac64";;
          *86*darwin*) matlab_arch_test="maci maci64";; 
          *solaris*) matlab_arch_test="sol sol64";;
          *cygwin*) matlab_arch_test="win32 win64";;
          *mingw*) matlab_arch_test="win32 win64";;
          *) AC_MSG_ERROR([Cannot guess Matlab architecture based on host type.]);;
        esac
        AC_MSG_CHECKING([for architecture])
        for matlab_arch in "$matlab_arch_test"; do
          if test -d "${matlab_bin_dir}/${matlab_arch}" -a -f "${matlab_bin_dir}/${matlab_arch}/MATLAB"; then
            AC_MSG_RESULT([${matlab_arch}])
            break
          fi
          matlab_arch="unkown"
          AC_MSG_RESULT([unknown])
        done
      fi

      # mex file extension or architecture found
      if test "x${matlab_mexext}" = "xunknown" -a "x${matlab_arch}" = "xunknown"; then
        AC_MSG_ERROR([Could not determine mex file extension nor Matlab architecture. Please supply a valid architecture flag using the option --with-matlab-arch])
      fi

      if test "x${matlab_arch}" = "xunknown"; then
        # architecture for mex file extension 
        AC_MSG_CHECKING([for architecture])
        case $matlab_mexext in
          mexglx) matlab_arch="glnx86";;
          mexa64) matlab_arch="glnxa64";;
          mexmac) matlab_arch="mac";;
          mexmaci) matlab_arch="maci";;
          mexmaci64) matlab_arch="maci64";;
          mexs64) matlab_arch="sol64";;
          mexw32) matlab_arch="win32";;
          mexw64) matlab_arch="win64";;
          *) AC_MSG_ERROR([Unsupported mex file extension ${matlab_mexext}.]);;
        esac
        AC_MSG_RESULT([${matlab_arch}])
      elif test "x${matlab_mexext}" = "xunknown"; then
        AC_MSG_CHECKING([for mex file extension])
        case $matlab_arch in
          glnx86) matlab_mexext="mexglx";;
          glnxa64) matlab_mexext="mexa64";;
          mac) matlab_mexext="mexmac";;
          maci) matlab_mexext="mexmaci";;
          maci64) matlab_mexext="mexmaci64";;
          sol64) matlab_mexext="mexs64";;
          win32) matlab_mexext="mexw32";;
          win64) matlab_mexext="mexw64";;
          *) AC_MSG_ERROR([Unsupported or invalid architecture ${matlab_arch}.]);;
        esac
        AC_MSG_RESULT([${matlab_mexext}])
      fi
    fi

    # add "." to mex file extension
    matlab_mexext=".$matlab_mexext"

    # arch bin dir
    matlab_arch_bin_dir="${matlab_bin_dir}/${matlab_arch}"
    AX_CHECK_DIR([${matlab_arch_bin_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])

    # dynamic library extension for architecture
    case $matlab_arch in
      glnx86|glnxa64|sol|sol64) matlab_libext=".so";;
      mac|mac64|maci|maci64) matlab_libext=".dylib";;
      win32|win64) matlab_libext=".dll";;
      *) AC_MSG_ERROR([Unsupported or invalid architecture ${matlab_arch}.]);;
    esac

    # libraries
    matlab_LDFLAGS="-L${matlab_arch_bin_dir}"
    saved_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS ${matlab_LDFLAGS}"
    matlab_LIBS=""
    AC_CHECK_LIB([mx],[mxMalloc],[matlab_LIBS="$matlab_LIBS -lmx"],
      [AC_CHECK_FILE([${matlab_arch_bin_dir}/libmx${matlab_libext}],
      [matlab_LIBS="$matlab_LIBS -lmx"],
      [AC_MSG_ERROR([Needed Matlab library mx not usable.])])],[])
    AC_CHECK_LIB([mex],[mexCallMATLAB],[matlab_LIBS="$matlab_LIBS -lmex"],
      [AC_CHECK_FILE([${matlab_arch_bin_dir}/libmex${matlab_libext}],
      [matlab_LIBS="$matlab_LIBS -lmex"],
      [AC_MSG_ERROR([Needed Matlab library mex not usable.])])],[])
    AC_CHECK_LIB([mat],[matGetVariable],[matlab_LIBS="$matlab_LIBS -lmat"],
      [AC_CHECK_FILE([${matlab_arch_bin_dir}/libmat${matlab_libext}],
      [matlab_LIBS="$matlab_LIBS -lmat"],
      [AC_MSG_ERROR([Needed Matlab library mat not usable.])])],[])
    LDFLAGS="$saved_LDFLAGS"

    matlab_CPPFLAGS="-I${matlab_include_dir}"

    # mexversion.c
    AC_CHECK_FILE([${matlab_src_dir}/mexversion.c],[matlab_CPPFLAGS="${matlab_CPPFLAGS} -I${matlab_src_dir}"; AC_DEFINE([HAVE_MEXVERSION_C],[1],[Define to have the file mexversion.c])],[AC_MSG_WARN([Required file ]${matlab_src_dir}[/mexversion.c not found])])

    ax_prog_matlab="yes"
  fi
])
