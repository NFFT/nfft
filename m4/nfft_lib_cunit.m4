# SYNOPSIS
#
#   NFFT_LIB_CUNIT
#
# DESCRIPTION
#
#   This macro tries to find out if the CUnit library is present and usable.
#
#   On success, it sets the nfft_cunit_have_lib output variable to "yes",
#   the nfft_cunit_CPPFLAGS output variable to the preprocessor flags necessary
#   to find the CUnit header files, and the nfft_cunit_LDFLAGS and
#   nfft_cunit_LIBS output variables to the linker flags necessary to link
#   against the CUnit library.
#
#   On failure, it sets the nfft_cunit_have_lib output variable to "no".
#   The content of the other output variables is undefined.
#
# LICENSE
#
#   Copyright (c) 2025 Jens Keiner <jens.keiner@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([NFFT_LIB_CUNIT],
[
  # Save current flags.
  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBSS="$LIBS"

  cunit_have_lib="no"
  cunit_CPPFLAGS=""
  cunit_LDFLAGS=""
  cunit_LIBS=""

  AC_ARG_WITH(cunit-includedir, [AS_HELP_STRING([--with-cunit-includedir=DIR],
  [compile with CUnit include directory DIR])], cunit_include_dir=$withval, 
    cunit_include_dir="")

  AC_ARG_WITH(cunit-libdir, [AS_HELP_STRING([--with-cunit-libdir=DIR],
  [compile with CUnit library directory DIR])], cunit_lib_dir=$withval, 
    cunit_lib_dir="")

  if test "x${cunit_include_dir}" != "x"; then
    AX_CHECK_DIR([${cunit_include_dir}],[],
      [AC_MSG_ERROR([The given CUnit include directory ${cunit_include_dir} does not exist.])])
    cunit_CPPFLAGS="-I${cunit_include_dir}"
  fi

  if test "x${cunit_lib_dir}" != "x"; then 
    AX_CHECK_DIR([${cunit_lib_dir}],[],
      [AC_MSG_ERROR([The given CUnit library directory ${cunit_lib_dir} does not exist.])])
    cunit_LDFLAGS="-L${cunit_lib_dir}"
  fi

  CPPFLAGS="${cunit_CPPFLAGS} ${saved_CPPFLAGS}"
  AC_CHECK_HEADER([CUnit/CUnit.h], [cunit_have_header=yes], [cunit_have_header=no])

  if test "x$cunit_have_header" = "xyes"; then
    LDFLAGS="${cunit_LDFLAGS} ${saved_LDFLAGS}"

    # Check if CUnit library is present and usable.
    AS_UNSET([ac_cv_search_cunit])
    AC_SEARCH_LIBS([CU_initialize_registry], [cunit], [], [], [])

    if test "x$ac_res" != "xno"; then
      cunit_have_lib="yes"
      if test "x$ac_res" != "xnone required"; then
        cunit_LIBS="$ac_res"
      fi
    fi
  fi

  # Restore saved flags.
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"
])dnl NFFT_CUNIT
