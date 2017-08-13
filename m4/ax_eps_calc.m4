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
# @synopsis AX_EPS_CALC
# @summary Calculate machine epsilon for given floating-point type (d = double, 
#     f = float, l = long double).
# @category C
#
# @version 2017-08-13
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>
AC_DEFUN([AX_EPS_CALC],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for floating-point epsilon as per calculation, ax_cv_eps_calc,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>
#define d 1
#define l 2
#define s 3
#define PRECISION $1
 
#if PRECISION == 1
typedef double R;
#define __FE__ "%.16lE"
#define K(x) ((R) x)
#elif PRECISION == 2
typedef long double R;
#define __FE__ "%.32LE"
#define K(x) ((R) x##L)
#elif PRECISION == 3
typedef float R;
#define __FE__ "%.8E"
#define K(x) ((R) x)
#else
#error "Unknown floating-point precision."
#endif], [  FILE *f;
  R eps = K(1.0);
  int i = 0;
  for (i = 0; i < 100000; i++) {
    if (K(1.0) + eps == K(1.0))
      break;
    eps /= K(2.0);
  }
  if (eps == K(0.0))
    return 1;
  if (K(1.0) + eps != K(1.0))
    return 1;
  eps *= K(2.0);
  f = fopen("conftest_eps_calc", "w"); if (!f) return 1;
  fprintf(f, __FE__ "\n", eps);
  fclose(f);
  return 0;])], 
  [ax_cv_eps_calc=`cat conftest_eps_calc`; rm -f conftest_eps_calc],
  [ax_cv_eps_calc=unknown; rm -f conftest_eps_calc],
  [ax_cv_eps_calc=unknown])])
AC_LANG_POP([C])
])
