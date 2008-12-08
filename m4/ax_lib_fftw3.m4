# $Id: matlab.m4 2652 2008-12-04 13:19:40Z keiner $
#
# Copyright (c) 2008 Jens Keiner
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
# @synopsis AX_LIB_FFTW3
# @summary Check configure options and assign variables related to the fftw3 library.
# @category C
#
# @version 2008-12-07
# @license GPLWithACException
# @author Jens Keiner <keiner@math.uni-luebeck.de>
#
#  If we find the library, set the shell variable `ax_have_fftw3' to `yes'.
#  Otherwise, set `ax_have_fftw3' to `no'.

AC_DEFUN([AX_LIB_FFTW3],
[
  AC_ARG_WITH(fftw3, [AC_HELP_STRING([--with-fftw3=DIR],
  [compile with fftw3 in DIR])], with_fftw3=$withval, with_fftw3="yes")

  AC_ARG_WITH(fftw3-libdir, [AC_HELP_STRING([--with-fftw3-libdir=DIR],
  [compile with fftw3 library directory DIR])], fftw3_lib_dir=$withval, 
    fftw3_lib_dir="yes")

  AC_ARG_WITH(fftw3-includedir, [AC_HELP_STRING([--with-fftw3-includedir=DIR],
  [compile with fftw3 include directory DIR])], fftw3_include_dir=$withval, 
    fftw3_include_dir="yes")

  if test "x$with_fftw3" = "xyes"; then
    if test "x${fftw3_include_dir}" = "xyes"; then 
      fftw3_include_dir="/usr/local/include"
    fi
    if test "x${fftw3_lib_dir}" = "xyes"; then 
      fftw3_lib_dir="/usr/local/lib"
    fi
  else
    fftw3_include_dir="$with_fftw3/include"
    fftw3_lib_dir="$with_fftw3/lib"
  fi

  AX_CHECK_DIR([${fftw3_include_dir}],[],
      [AC_MSG_ERROR([The directory ${fftw3_include_dir} does not exist.])])

  AX_CHECK_DIR([${fftw3_lib_dir}],[],
      [AC_MSG_ERROR([The directory ${fftw3_lib_dir} does not exist.])])

  # Save compiler and linker flags.
  fftw3_CPPFLAGS="-I$fftw3_include_dir"
  fftw3_LDFLAGS="-L$fftw3_lib_dir"
  fftw3_LIBS="-lfftw3"
  saved_LDFLAGS="$LDFLAGS"
  saved_CPPFLAGS="$CPPFLAGS"
  saved_LIBS="$LIBS"
  CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"
  LDFLAGS="$LDFLAGS $fftw3_LDFLAGS"
  LIBS="$LIBS $fftw3_LIBS"

  # Check if library is present and usable.
  AC_CHECK_HEADER([fftw3.h],
    [AC_CHECK_LIB(fftw3, fftw_execute,
      ax_lib_fftw3=yes,
      ax_lib_fftw3=no)], ax_lib_fftw3=no)

  # Check if library was found.
  ax_lib_fftw3_libtool="no"
  if test "x$ax_lib_fftw3" = "xyes"; then
    if test -f "$fftw3_lib_dir/libfftw3.la"; then
        ax_lib_fftw3_libtool="yes"
    elif test "x$fftw3_lib_dir" = "x/usr/local/lib" -a -f "/usr/lib/libfftw3.la"; then
      fftw3_include_dir="/usr/include"
      fftw3_CFLAGS="-I$fftw3_include_dir"
      fftw3_lib_dir="/usr/lib"
      fftw3_LDFLAGS="-L$fftw3_lib_dir"
      ax_lib_fftw3_libtool="yes"
    fi
  fi

  # Restore saved flags.
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"
])
