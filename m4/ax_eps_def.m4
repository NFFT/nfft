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
# @summary Determine machine epsilon for given floating-point type (d = double, 
#     f = float, l = long double) from 2^(1-MAND_DIG) defined in float.h.
# @category C
#
# @version 2017-08-13
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>
AC_DEFUN([AX_EPS_DEF],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for floating-point epsilon as per length of mantissa, ax_cv_eps_def,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>
#include <float.h>
#define d 1
#define l 2
#define s 3
#define PRECISION $1

#if PRECISION == 1
typedef double R;
#define __FE__ "%.16lE"
#define K(x) ((R) x)
#define MANT_DIG DBL_MANT_DIG
#elif PRECISION == 2
typedef long double R;
#define __FE__ "%.32LE"
#define K(x) ((R) x##L)
#define MANT_DIG LDBL_MANT_DIG
#elif PRECISION == 3
typedef float R;
#define __FE__ "%.8E"
#define K(x) ((R) x)
#define MANT_DIG FLT_MANT_DIG
#else
#error "Unknown floating-point precision."
#endif], [  FILE *f;
  f = fopen("conftest_eps_def", "w"); if (!f) return 1;
  R epsilon = K(1.0);
  for (int i=0; i<MANT_DIG-1; i++)
    epsilon /= K(2.0);
  fprintf(f, __FE__ "\n", epsilon);
  fclose(f);
  return 0;])], 
  [ax_cv_eps_def=`cat conftest_eps_def`; rm -f conftest_eps_def],
  [ax_cv_eps_def=unknown; rm -f conftest_eps_def],
  [ax_cv_eps_def=unknown])])
AC_LANG_POP([C])
])
