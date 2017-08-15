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
# @synopsis AX_EPS_DEF
# @summary Dtermine machine epsilon for given floating-point type (d = double, 
#     f = float, l = long double) from corresponding macro defined in float.h.
# @category C
#
# @version 2017-08-13
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>
AC_DEFUN([AX_EPS_DEF],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for floating-point epsilon as per macro, ax_cv_eps_def,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>
#include <float.h>
#define d 1
#define l 2
#define s 3
#define PRECISION $1
 
#if PRECISION == 1
#define __FE__ "%.16lE"
#define EPSILON DBL_EPSILON
#elif PRECISION == 2
#define __FE__ "%.32LE"
#define EPSILON LDBL_EPSILON
#elif PRECISION == 3
#define __FE__ "%.8E"
#define EPSILON FLT_EPSILON
#else
#error "Unknown floating-point precision."
#endif], [  FILE *f;
  f = fopen("conftest_eps_def", "w"); if (!f) return 1;
  fprintf(f, __FE__ "\n", EPSILON);
  fclose(f);
  return 0;])], 
  [ax_cv_eps_def=`cat conftest_eps_def`; rm -f conftest_eps_def],
  [ax_cv_eps_def=unknown; rm -f conftest_eps_def],
  [ax_cv_eps_def=unknown])])
AC_LANG_POP([C])
])
