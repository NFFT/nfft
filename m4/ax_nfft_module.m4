# Copyright (c) 2009-2010, 2014, 2021 The NFFT Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/>.

AC_DEFUN([AX_NFFT_MODULE],
[
  nfft_module_default_local=$nfft_module_default
  if test "x$PRECISION" = "xs" -o "x$PRECISION" = "xl"; then
    if test "x$4" = "xno"; then
      nfft_module_default_local="no"
    fi
  fi
  AC_ARG_ENABLE($1, [AS_HELP_STRING([--enable-]$1,[build ]$2[ module (]$3[)])], 
    ok=$enableval, ok=$nfft_module_default_local)
  AC_MSG_CHECKING([Whether to compile $2 module])
  if m4_default($7,test "x$ok" = "xyes"); then
    AC_MSG_RESULT([yes])
    if test "x$PRECISION" = "xs" -o "x$PRECISION" = "xl"; then
      if test "x$4" = "xno"; then
        AC_MSG_ERROR([The $2 module cannot be used with the selected floating point precision.])
      fi
    fi
    AC_DEFINE(HAVE_$2, 1, [Define to enable ]$2[module.])
    HAVE_$2="#define HAVE_$2 1"
    $5
  else
    AC_MSG_RESULT([no])
    HAVE_$2="#undef HAVE_$2"
    $6
  fi
  AM_CONDITIONAL(HAVE_$2, m4_default($7,test "x$ok" = "xyes"))
  AC_SUBST(HAVE_$2)
])