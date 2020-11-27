/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*! \file kernels.c
 *  \brief File with predefined kernels for the fast summation algorithm.
 */
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "kernels.h"

/**
 * \addtogroup applications_fastsum
 * \{
 */

C gaussian(R x, int der, const R *param)    /* K(x)=EXP(-x^2/c^2) */
{
  R c = param[0];
  R value = K(0.0);

  switch (der)
  {
    case  0 : value = EXP(-x*x/(c*c)); break;
    case  1 : value = -K(2.0)*x/(c*c)*EXP(-x*x/(c*c)); break;
    case  2 : value = K(2.0)*EXP(-x*x/(c*c))*(-c*c+K(2.0)*x*x)/(c*c*c*c); break;
    case  3 : value = -K(4.0)*x*EXP(-x*x/(c*c))*(-K(3.0)*c*c+K(2.0)*x*x)/(c*c*c*c*c*c); break;
    case  4 : value = K(4.0)*EXP(-x*x/(c*c))*(K(3.0)*c*c*c*c-K(12.0)*c*c*x*x+K(4.0)*x*x*x*x)/(c*c*c*c*c*c*c*c); break;
    case  5 : value = -K(8.0)*x*EXP(-x*x/(c*c))*(K(15.0)*c*c*c*c-K(20.0)*c*c*x*x+K(4.0)*x*x*x*x)/POW(c,K(10.0)); break;
    case  6 : value = K(8.0)*EXP(-x*x/(c*c))*(-K(15.0)*c*c*c*c*c*c+K(90.0)*x*x*c*c*c*c-K(60.0)*x*x*x*x*c*c+K(8.0)*x*x*x*x*x*x)/POW(c,K(12.0)); break;
    case  7 : value = -K(16.0)*x*EXP(-x*x/(c*c))*(-K(105.0)*c*c*c*c*c*c+K(210.0)*x*x*c*c*c*c-K(84.0)*x*x*x*x*c*c+K(8.0)*x*x*x*x*x*x)/POW(c,K(14.0)); break;
    case  8 : value = K(16.0)*EXP(-x*x/(c*c))*(K(105.0)*c*c*c*c*c*c*c*c-K(840.0)*x*x*c*c*c*c*c*c+K(840.0)*x*x*x*x*c*c*c*c-K(224.0)*x*x*x*x*x*x*c*c+K(16.0)*x*x*x*x*x*x*x*x)/POW(c,K(16.0)); break;
    case  9 : value = -K(32.0)*x*EXP(-x*x/(c*c))*(K(945.0)*c*c*c*c*c*c*c*c-K(2520.0)*x*x*c*c*c*c*c*c+K(1512.0)*x*x*x*x*c*c*c*c-K(288.0)*x*x*x*x*x*x*c*c+K(16.0)*x*x*x*x*x*x*x*x)/POW(c,K(18.0)); break;
    case 10 : value = K(32.0)*EXP(-x*x/(c*c))*(-K(945.0)*POW(c,K(10.0))+K(9450.0)*x*x*c*c*c*c*c*c*c*c-K(12600.0)*x*x*x*x*c*c*c*c*c*c+K(5040.0)*x*x*x*x*x*x*c*c*c*c-K(720.0)*x*x*x*x*x*x*x*x*c*c+K(32.0)*POW(x,K(10.0)))/POW(c,K(20.0)); break;
    case 11 : value = -K(64.0)*x*EXP(-x*x/(c*c))*(-K(10395.0)*POW(c,K(10.0))+K(34650.0)*x*x*c*c*c*c*c*c*c*c-K(27720.0)*x*x*x*x*c*c*c*c*c*c+K(7920.0)*x*x*x*x*x*x*c*c*c*c-K(880.0)*x*x*x*x*x*x*x*x*c*c+K(32.0)*POW(x,K(10.0)))/POW(c,K(22.0)); break;
    case 12 : value = K(64.0)*EXP(-x*x/(c*c))*(K(10395.0)*POW(c,K(12.0))-K(124740.0)*x*x*POW(c,K(10.0))+K(207900.0)*x*x*x*x*c*c*c*c*c*c*c*c-K(110880.0)*x*x*x*x*x*x*c*c*c*c*c*c+K(23760.0)*x*x*x*x*x*x*x*x*c*c*c*c-K(2112.0)*POW(x,K(10.0))*c*c+K(64.0)*POW(x,K(12.0)))/POW(c,K(24.0)); break;
    default : value = K(0.0);
  }

  return value;
}

C multiquadric(R x, int der, const R *param)    /* K(x)=SQRT(x^2+c^2) */
{
  R c=param[0];
  R value=K(0.0);

  switch (der)
  {
    case  0 : value=SQRT(x*x+c*c); break;
    case  1 : value=K(1.0)/(SQRT(x*x+c*c))*x; break;
    case  2 : value=c*c/SQRT(POW(x*x+c*c,K(3.0))); break;
    case  3 : value=-K(3.0)*x*c*c/SQRT(POW(x*x+c*c,K(5.0))); break;
    case  4 : value=K(3.0)*c*c*(K(4.0)*x*x-c*c)/SQRT(POW(x*x+c*c,K(7.0))); break;
    case  5 : value=-K(15.0)*x*c*c*(K(4.0)*x*x-K(3.0)*c*c)/SQRT(POW(x*x+c*c,K(9.0))); break;
    case  6 : value=K(45.0)*c*c*(K(8.0)*x*x*x*x-K(12.0)*x*x*c*c+c*c*c*c)/SQRT(POW(x*x+c*c,K(11.0))); break;
    case  7 : value=-K(315.0)*x*c*c*(K(8.0)*x*x*x*x-K(20.0)*x*x*c*c+K(5.0)*c*c*c*c)/SQRT(POW(x*x+c*c,K(13.0))); break;
    case  8 : value=K(315.0)*c*c*(K(64.0)*x*x*x*x*x*x-K(240.0)*x*x*x*x*c*c+K(120.0)*x*x*c*c*c*c-K(5.0)*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(15.0))); break;
    case  9 : value=-K(2835.0)*x*c*c*(K(64.0)*x*x*x*x*x*x-K(336.0)*x*x*x*x*c*c+K(280.0)*x*x*c*c*c*c-K(35.0)*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(17.0))); break;
    case 10 : value=K(14175.0)*c*c*(K(128.0)*x*x*x*x*x*x*x*x-K(896.0)*x*x*x*x*x*x*c*c+K(1120.0)*x*x*x*x*c*c*c*c-K(280.0)*x*x*c*c*c*c*c*c+K(7.0)*c*c*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(19.0))); break;
    case 11 : value=-K(155925.0)*x*c*c*(K(128.0)*x*x*x*x*x*x*x*x-K(1152.0)*x*x*x*x*x*x*c*c+K(2016.0)*x*x*x*x*c*c*c*c-K(840.0)*x*x*c*c*c*c*c*c+K(63.0)*c*c*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(21.0))); break;
    case 12 : value=K(467775.0)*c*c*(K(1260.0)*x*x*c*c*c*c*c*c*c*c-K(21.0)*POW(c,K(10.0))+K(512.0)*POW(x,K(10.0))-K(5760.0)*x*x*x*x*x*x*x*x*c*c+K(13440.0)*x*x*x*x*x*x*c*c*c*c-K(8400.0)*x*x*x*x*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(23.0))); break;
    default : value=K(0.0);
  }

  return value;
}

C inverse_multiquadric(R x, int der, const R *param)    /* K(x)=1/SQRT(x^2+c^2) */
{
  R c=param[0];
  R value=K(0.0);

  switch (der)
  {
    case  0 : value=K(1.0)/SQRT(x*x+c*c); break;
    case  1 : value=-K(1.0)/(SQRT(POW(x*x+c*c,K(3.0))))*x; break;
    case  2 : value=(K(2.0)*x*x-c*c)/SQRT(POW(x*x+c*c,K(5.0))); break;
    case  3 : value=-K(3.0)*x*(K(2.0)*x*x-K(3.0)*c*c)/SQRT(POW(x*x+c*c,K(7.0))); break;
    case  4 : value=K(3.0)*(K(8.0)*x*x*x*x-K(24.0)*x*x*c*c+K(3.0)*c*c*c*c)/SQRT(POW(x*x+c*c,K(9.0))); break;
    case  5 : value=-K(15.0)*x*(K(8.0)*x*x*x*x-K(40.0)*x*x*c*c+K(15.0)*c*c*c*c)/SQRT(POW(x*x+c*c,K(11.0))); break;
    case  6 : value=K(45.0)*(K(16.0)*x*x*x*x*x*x-K(120.0)*x*x*x*x*c*c+K(90.0)*x*x*c*c*c*c-K(5.0)*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(13.0))); break;
    case  7 : value=-K(315.0)*x*(K(16.0)*x*x*x*x*x*x-K(168.0)*x*x*x*x*c*c+K(210.0)*x*x*c*c*c*c-K(35.0)*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(15.0))); break;
    case  8 : value=K(315.0)*(K(128.0)*x*x*x*x*x*x*x*x-K(1792.0)*x*x*x*x*x*x*c*c+K(3360.0)*x*x*x*x*c*c*c*c-K(1120.0)*x*x*c*c*c*c*c*c+K(35.0)*c*c*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(17.0))); break;
    case  9 : value=-K(2835.0)*x*(K(128.0)*x*x*x*x*x*x*x*x-K(2304.0)*x*x*x*x*x*x*c*c+K(6048.0)*x*x*x*x*c*c*c*c-K(3360.0)*x*x*c*c*c*c*c*c+K(315.0)*c*c*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(19.0))); break;
    case 10 : value=K(14175.0)*(K(256.0)*POW(x,K(10.0))-K(5760.0)*x*x*x*x*x*x*x*x*c*c+K(20160.0)*x*x*x*x*x*x*c*c*c*c-K(16800.0)*x*x*x*x*c*c*c*c*c*c+K(3150.0)*x*x*c*c*c*c*c*c*c*c-K(63.0)*POW(c,K(10.0)))/SQRT(POW(x*x+c*c,K(21.0))); break;
    case 11 : value=-K(155925.0)*x*(K(256.0)*POW(x,K(10.0))-K(7040.0)*x*x*x*x*x*x*x*x*c*c+K(31680.0)*x*x*x*x*x*x*c*c*c*c-K(36960.0)*x*x*x*x*c*c*c*c*c*c+K(11550.0)*x*x*c*c*c*c*c*c*c*c-K(693.0)*POW(c,K(10.0)))/SQRT(POW(x*x+c*c,K(23.0))); break;
    case 12 : value=K(467775.0)*(K(231.0)*POW(c,K(12.0))+K(190080.0)*x*x*x*x*x*x*x*x*c*c*c*c-K(16632.0)*x*x*POW(c,K(10.0))-K(295680.0)*x*x*x*x*x*x*c*c*c*c*c*c+K(138600.0)*x*x*x*x*c*c*c*c*c*c*c*c+K(1024.0)*POW(x,K(12.0))-K(33792.0)*POW(x,K(10.0))*c*c)/SQRT(POW(x*x+c*c,K(25.0))); break;
    default : value=K(0.0);
  }

  return value;
}

C logarithm(R x, int der, const R *param)    /* K(x)=LOG |x| */
{
  R value=K(0.0);

  (void)param;

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value=LOG(FABS(x)); break;
    case  1 : value=(x<0 ? -1 : 1)/FABS(x); break;
    case  2 : value=-1/(x*x); break;
    case  3 : value=K(2.0)*(x<0 ? -1 : 1)/POW(FABS(x),K(3.0)); break;
    case  4 : value=-K(6.0)/(x*x*x*x); break;
    case  5 : value=K(24.0)*(x<0 ? -1 : 1)/POW(FABS(x),K(5.0)); break;
    case  6 : value=-K(120.0)/(x*x*x*x*x*x); break;
    case  7 : value=K(720.0)*(x<0 ? -1 : 1)/POW(FABS(x),K(7.0)); break;
    case  8 : value=-K(5040.0)/(x*x*x*x*x*x*x*x); break;
    case  9 : value=K(40320.0)*(x<0 ? -1 : 1)/POW(FABS(x),K(9.0)); break;
    case 10 : value=-K(362880.0)/POW(x,K(10.0)); break;
    case 11 : value=K(3628800.0)*(x<0 ? -1 : 1)/POW(FABS(x),K(11.0)); break;
    case 12 : value=-K(39916800.0)/POW(x,K(12.0)); break;
    case 13 : value=K(479001600.0)/POW(x,K(13.0)); break;
    case 14 : value=-K(6227020800.0)/POW(x,K(14.0)); break;
    case 15 : value=K(87178291200.0)/POW(x,K(15.0)); break;
    case 16 : value=-K(1307674368000.0)/POW(x,K(16.0)); break;
    case 17 : value=K(20922789888000.0)/POW(x,K(17.0)); break;
    default : value=K(0.0);
  }

  return value;
}

C thinplate_spline(R x, int der, const R *param)    /* K(x) = x^2 LOG |x| */
{
  R value=K(0.0);

  (void)param;

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value=x*x*LOG(FABS(x)); break;
    case  1 : value=K(2.0)*x*LOG(FABS(x))+x; break;
    case  2 : value=K(2.0)*LOG(FABS(x))+K(3.0); break;
    case  3 : value=K(2.0)/x; break;
    case  4 : value=-K(2.0)/(x*x); break;
    case  5 : value=K(4.0)/(x*x*x); break;
    case  6 : value=-K(12.0)/(x*x*x*x); break;
    case  7 : value=K(48.0)/(x*x*x*x*x); break;
    case  8 : value=-K(240.0)/(x*x*x*x*x*x); break;
    case  9 : value=K(1440.0)/(x*x*x*x*x*x*x); break;
    case 10 : value=-K(10080.0)/(x*x*x*x*x*x*x*x); break;
    case 11 : value=K(80640.0)/(x*x*x*x*x*x*x*x*x); break;
    case 12 : value=-K(725760.0)/POW(x,K(10.0)); break;
    default : value=K(0.0);
  }

  return value;
}

C one_over_square(R x, int der, const R *param)    /* K(x) = 1/x^2 */
{
  R value=K(0.0);

  (void)param;

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value=K(1.0)/(x*x); break;
    case  1 : value=-K(2.0)/(x*x*x); break;
    case  2 : value=K(6.0)/(x*x*x*x); break;
    case  3 : value=-K(24.0)/(x*x*x*x*x); break;
    case  4 : value=K(120.0)/(x*x*x*x*x*x); break;
    case  5 : value=-K(720.0)/(x*x*x*x*x*x*x); break;
    case  6 : value=K(5040.0)/(x*x*x*x*x*x*x*x); break;
    case  7 : value=-K(40320.0)/(x*x*x*x*x*x*x*x*x); break;
    case  8 : value=K(362880.0)/POW(x,K(10.0)); break;
    case  9 : value=-K(3628800.0)/POW(x,K(11.0)); break;
    case 10 : value=K(39916800.0)/POW(x,K(12.0)); break;
    case 11 : value=-K(479001600.0)/POW(x,K(13.0)); break;
    case 12 : value=K(6227020800.0)/POW(x,K(14.0)); break;
    default : value=K(0.0);
  }

  return value;
}

C one_over_modulus(R x, int der, const R *param)    /* K(x) = 1/|x| */
{
  R value=K(0.0);

  (void)param;

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value=K(1.0)/FABS(x); break;
    case  1 : value=-1/x/FABS(x); break;
    case  2 : value=K(2.0)/POW(FABS(x),K(3.0)); break;
    case  3 : value=-K(6.0)/(x*x*x)/FABS(x); break;
    case  4 : value=K(24.0)/POW(FABS(x),K(5.0)); break;
    case  5 : value=-K(120.0)/(x*x*x*x*x)/FABS(x); break;
    case  6 : value=K(720.0)/POW(FABS(x),K(7.0)); break;
    case  7 : value=-K(5040.0)/(x*x*x*x*x*x*x)/FABS(x); break;
    case  8 : value=K(40320.0)/POW(FABS(x),K(9.0)); break;
    case  9 : value=-K(362880.0)/(x*x*x*x*x*x*x*x*x)/FABS(x); break;
    case 10 : value=K(3628800.0)/POW(FABS(x),K(11.0)); break;
    case 11 : value=-K(39916800.0)/POW(x,K(11.0))/FABS(x); break;
    case 12 : value=K(479001600.0)/POW(FABS(x),K(13.0)); break;
    default : value=K(0.0);
  }

  return value;
}

C one_over_x(R x, int der, const R *param)    /* K(x) = 1/x */
{
  R value=K(0.0);

  (void)param;

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value=K(1.0)/x; break;
    case  1 : value=-K(1.0)/(x*x); break;
    case  2 : value=K(2.0)/(x*x*x); break;
    case  3 : value=-K(6.0)/(x*x*x*x); break;
    case  4 : value=K(24.0)/(x*x*x*x*x); break;
    case  5 : value=-K(120.0)/(x*x*x*x*x*x); break;
    case  6 : value=K(720.0)/(x*x*x*x*x*x*x); break;
    case  7 : value=-K(5040.0)/(x*x*x*x*x*x*x*x); break;
    case  8 : value=K(40320.0)/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-K(362880.0)/POW(x,K(10.0)); break;
    case 10 : value=K(3628800.0)/POW(x,K(11.0)); break;
    case 11 : value=-K(39916800.0)/POW(x,K(12.0)); break;
    case 12 : value=K(479001600.0)/POW(x,K(13.0)); break;
    default : value=K(0.0);
  }

  return value;
}

C inverse_multiquadric3(R x, int der, const R *param)    /* K(x) = 1/SQRT(x^2+c^2)^3 */
{
  R c=param[0];
  R value=K(0.0);

  switch (der)
  {
    case  0 : value=K(1.0)/(SQRT(POW(x*x+c*c,K(3.0)))); break;
    case  1 : value=-K(3.0)/SQRT(POW(x*x+c*c,K(5.0)))*x; break;
    case  2 : value=K(3.0)*(K(4.0)*x*x-c*c)/SQRT(POW(x*x+c*c,K(7.0))); break;
    case  3 : value=-K(15.0)*x*(K(4.0)*x*x-K(3.0)*c*c)/SQRT(POW(x*x+c*c,K(9.0))); break;
    case  4 : value=K(45.0)*(K(8.0)*x*x*x*x-K(12.0)*x*x*c*c+c*c*c*c)/SQRT(POW(x*x+c*c,K(11.0))); break;
    case  5 : value=-K(315.0)*x*(K(8.0)*x*x*x*x-K(20.0)*x*x*c*c+K(5.0)*c*c*c*c)/SQRT(POW(x*x+c*c,K(13.0))); break;
    case  6 : value=K(315.0)*(K(64.0)*x*x*x*x*x*x-K(240.0)*x*x*x*x*c*c+K(120.0)*x*x*c*c*c*c-K(5.0)*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(15.0))); break;
    case  7 : value=-K(2835.0)*x*(K(64.0)*x*x*x*x*x*x-K(336.0)*x*x*x*x*c*c+K(280.0)*x*x*c*c*c*c-K(35.0)*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(17.0))); break;
    case  8 : value=K(14175.0)*(K(128.0)*x*x*x*x*x*x*x*x-K(896.0)*x*x*x*x*x*x*c*c+K(1120.0)*x*x*x*x*c*c*c*c-K(280.0)*x*x*c*c*c*c*c*c+K(7.0)*c*c*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(19.0))); break;
    case  9 : value=-K(155925.0)*x*(K(128.0)*x*x*x*x*x*x*x*x-K(1152.0)*x*x*x*x*x*x*c*c+K(2016.0)*x*x*x*x*c*c*c*c-K(840.0)*x*x*c*c*c*c*c*c+K(63.0)*c*c*c*c*c*c*c*c)/SQRT(POW(x*x+c*c,K(21.0))); break;
    case 10 : value=K(467775.0)*(K(512.0)*POW(x,K(10.0))-K(5760.0)*x*x*x*x*x*x*x*x*c*c+K(13440.0)*x*x*x*x*x*x*c*c*c*c-K(8400.0)*x*x*x*x*c*c*c*c*c*c+K(1260.0)*x*x*c*c*c*c*c*c*c*c-K(21.0)*POW(c,K(10.0)))/SQRT(POW(x*x+c*c,K(23.0))); break;
    case 11 : value=-K(6081075.0)*x*(K(512.0)*POW(x,K(10.0))-K(7040.0)*x*x*x*x*x*x*x*x*c*c+K(21120.0)*x*x*x*x*x*x*c*c*c*c-K(18480.0)*x*x*x*x*c*c*c*c*c*c+K(4620.0)*x*x*c*c*c*c*c*c*c*c-K(231.0)*POW(c,K(10.0)))/SQRT(POW(x*x+c*c,K(25.0))); break;
    case 12 : value=K(42567525.0)*(K(1024.0)*POW(x,K(12.0))+K(27720.0)*x*x*x*x*c*c*c*c*c*c*c*c+K(33.0)*POW(c,K(12.0))-K(2772.0)*x*x*POW(c,K(10.0))-K(73920.0)*x*x*x*x*x*x*c*c*c*c*c*c+K(63360.0)*x*x*x*x*x*x*x*x*c*c*c*c-K(16896.0)*POW(x,K(10.0))*c*c)/SQRT(POW(x*x+c*c,K(27.0))); break;
    default : value=K(0.0);
  }

  return value;
}

C sinc_kernel(R x, int der, const R *param)    /* K(x) = SIN(cx)/x */
{
  R c=param[0];
  R value=K(0.0);

  if (FABS(x)<EPSILON) value=c;
  else switch (der)
  {
    case  0 : value=SIN(c*x)/x; break;
    case  1 : value=(COS(c*x)*c*x-SIN(c*x))/(x*x); break;
    case  2 : value=-(SIN(c*x)*c*c*x*x+K(2.0)*COS(c*x)*c*x-K(2.0)*SIN(c*x))/(x*x*x); break;
    case  3 : value=-(COS(c*x)*c*c*c*x*x*x-K(3.0)*SIN(c*x)*c*c*x*x-K(6.0)*COS(c*x)*c*x+K(6.0)*SIN(c*x))/(x*x*x*x); break;
    case  4 : value=(SIN(c*x)*c*c*c*c*x*x*x*x+K(4.0)*COS(c*x)*c*c*c*x*x*x-K(12.0)*SIN(c*x)*c*c*x*x-K(24.0)*COS(c*x)*c*x+K(24.0)*SIN(c*x))/(x*x*x*x*x); break;
    case  5 : value=(COS(c*x)*c*c*c*c*c*x*x*x*x*x-K(5.0)*SIN(c*x)*c*c*c*c*x*x*x*x-K(20.0)*COS(c*x)*c*c*c*x*x*x+K(60.0)*SIN(c*x)*c*c*x*x+K(120.0)*COS(c*x)*c*x-K(120.0)*SIN(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+K(6.0)*COS(c*x)*c*c*c*c*c*x*x*x*x*x-K(30.0)*SIN(c*x)*c*c*c*c*x*x*x*x-K(120.0)*COS(c*x)*c*c*c*x*x*x+K(360.0)*SIN(c*x)*c*c*x*x+K(720.0)*COS(c*x)*c*x-K(720.0)*SIN(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=-(COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-K(7.0)*SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(42.0)*COS(c*x)*c*c*c*c*c*x*x*x*x*x+K(210.0)*SIN(c*x)*c*c*c*c*x*x*x*x+K(840.0)*COS(c*x)*c*c*c*x*x*x-K(2520.0)*SIN(c*x)*c*c*x*x-K(5040.0)*COS(c*x)*c*x+K(5040.0)*SIN(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+K(8.0)*COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-K(56.0)*SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(336.0)*COS(c*x)*c*c*c*c*c*x*x*x*x*x+K(1680.0)*SIN(c*x)*c*c*c*c*x*x*x*x+K(6720.0)*COS(c*x)*c*c*c*x*x*x-K(20160.0)*SIN(c*x)*c*c*x*x-K(40320.0)*COS(c*x)*c*x+K(40320.0)*SIN(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=(COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-K(9.0)*SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-K(72.0)*COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+K(504.0)*SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+K(3024.0)*COS(c*x)*c*c*c*c*c*x*x*x*x*x-K(15120.0)*SIN(c*x)*c*c*c*c*x*x*x*x-K(60480.0)*COS(c*x)*c*c*c*x*x*x+K(181440.0)*SIN(c*x)*c*c*x*x+K(362880.0)*COS(c*x)*c*x-K(362880.0)*SIN(c*x))/POW(x,K(10.0)); break;
    case 10 : value=-(SIN(c*x)*POW(c,K(10.0))*POW(x,K(10.0))+K(10.0)*COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-K(90.0)*SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-K(720.0)*COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+K(5040.0)*SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+K(30240.0)*COS(c*x)*c*c*c*c*c*x*x*x*x*x-K(151200.0)*SIN(c*x)*c*c*c*c*x*x*x*x-K(604800.0)*COS(c*x)*c*c*c*x*x*x+K(1814400.0)*SIN(c*x)*c*c*x*x+K(3628800.0)*COS(c*x)*c*x-K(3628800.0)*SIN(c*x))/POW(x,K(11.0)); break;
    case 11 : value=-(COS(c*x)*POW(c,K(11.0))*POW(x,K(11.0))-K(11.0)*SIN(c*x)*POW(c,K(10.0))*POW(x,K(10.0))-K(110.0)*COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+K(990.0)*SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+K(7920.0)*COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-K(55440.0)*SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(332640.0)*COS(c*x)*c*c*c*c*c*x*x*x*x*x+K(1663200.0)*SIN(c*x)*c*c*c*c*x*x*x*x+K(6652800.0)*COS(c*x)*c*c*c*x*x*x-K(19958400.0)*SIN(c*x)*c*c*x*x-K(39916800.0)*COS(c*x)*c*x+K(39916800.0)*SIN(c*x))/POW(x,K(12.0)); break;
    case 12 : value=(SIN(c*x)*POW(c,K(12.0))*POW(x,K(12.0))+K(12.0)*COS(c*x)*POW(c,K(11.0))*POW(x,K(11.0))-K(132.0)*SIN(c*x)*POW(c,K(10.0))*POW(x,K(10.0))-K(1320.0)*COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+K(11880.0)*SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+K(95040.0)*COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-K(665280.0)*SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(3991680.0)*COS(c*x)*c*c*c*c*c*x*x*x*x*x+K(19958400.0)*SIN(c*x)*c*c*c*c*x*x*x*x+K(79833600.0)*COS(c*x)*c*c*c*x*x*x-K(239500800.0)*SIN(c*x)*c*c*x*x-K(479001600.0)*COS(c*x)*c*x+K(479001600.0)*SIN(c*x))/POW(x,K(13.0)); break;
    default : value=K(0.0);
  }

  return value;
}

C cosc(R x, int der, const R *param)    /* K(x) = COS(cx)/x */
{
  R c=param[0];
  R value=K(0.0);
  R sign;

  if (x<0) sign=-K(1.0); else sign=K(1.0);
  x=FABS(x);

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value=COS(c*x)/x; break;
    case  1 : value=-(SIN(c*x)*c*x+COS(c*x))/(x*x); break;
    case  2 : value=(-COS(c*x)*c*c*x*x+K(2.0)*SIN(c*x)*c*x+K(2.0)*COS(c*x))/(x*x*x); break;
    case  3 : value=(SIN(c*x)*c*c*c*x*x*x+K(3.0)*COS(c*x)*c*c*x*x-K(6.0)*SIN(c*x)*c*x-K(6.0)*COS(c*x))/(x*x*x*x); break;
    case  4 : value=(COS(c*x)*c*c*c*c*x*x*x*x-K(4.0)*SIN(c*x)*c*c*c*x*x*x-K(12.0)*COS(c*x)*c*c*x*x+K(24.0)*SIN(c*x)*c*x+K(24.0)*COS(c*x))/(x*x*x*x*x); break;
    case  5 : value=-(SIN(c*x)*c*c*c*c*c*x*x*x*x*x+K(5.0)*COS(c*x)*c*c*c*c*x*x*x*x-K(20.0)*SIN(c*x)*c*c*c*x*x*x-K(60.0)*COS(c*x)*c*c*x*x+K(120.0)*SIN(c*x)*c*x+K(120.0)*COS(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(6.0)*SIN(c*x)*c*c*c*c*c*x*x*x*x*x-K(30.0)*COS(c*x)*c*c*c*c*x*x*x*x+K(120.0)*SIN(c*x)*c*c*c*x*x*x+K(360.0)*COS(c*x)*c*c*x*x-K(720.0)*SIN(c*x)*c*x-K(720.0)*COS(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=(SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+K(7.0)*COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(42.0)*SIN(c*x)*c*c*c*c*c*x*x*x*x*x-K(210.0)*COS(c*x)*c*c*c*c*x*x*x*x+K(840.0)*SIN(c*x)*c*c*c*x*x*x+K(2520.0)*COS(c*x)*c*c*x*x-K(5040.0)*SIN(c*x)*c*x-K(5040.0)*COS(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-K(8.0)*SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-K(56.0)*COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+K(336.0)*SIN(c*x)*c*c*c*c*c*x*x*x*x*x+K(1680.0)*COS(c*x)*c*c*c*c*x*x*x*x-K(6720.0)*SIN(c*x)*c*c*c*x*x*x-K(20160.0)*COS(c*x)*c*c*x*x+K(40320.0)*SIN(c*x)*c*x+K(40320.0)*COS(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-(SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+K(9.0)*COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-K(72.0)*SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-K(504.0)*COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+K(3024.0)*SIN(c*x)*c*c*c*c*c*x*x*x*x*x+K(15120.0)*COS(c*x)*c*c*c*c*x*x*x*x-K(60480.0)*SIN(c*x)*c*c*c*x*x*x-K(181440.0)*COS(c*x)*c*c*x*x+K(362880.0)*SIN(c*x)*c*x+K(362880.0)*COS(c*x))/POW(x,K(10.0)); break;
    case 10 : value=-(COS(c*x)*POW(c,K(10.0))*POW(x,K(10.0))-K(10.0)*SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-K(90.0)*COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+K(720.0)*SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+K(5040.0)*COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(30240.0)*SIN(c*x)*c*c*c*c*c*x*x*x*x*x-K(151200.0)*COS(c*x)*c*c*c*c*x*x*x*x+K(604800.0)*SIN(c*x)*c*c*c*x*x*x+K(1814400.0)*COS(c*x)*c*c*x*x-K(3628800.0)*SIN(c*x)*c*x-K(3628800.0)*COS(c*x))/POW(x,K(11.0)); break;
    case 11 : value=(SIN(c*x)*POW(c,K(11.0))*POW(x,K(11.0))+K(11.0)*COS(c*x)*POW(c,K(10.0))*POW(x,K(10.0))-K(110.0)*SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-K(990.0)*COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+K(7920.0)*SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+K(55440.0)*COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-K(332640.0)*SIN(c*x)*c*c*c*c*c*x*x*x*x*x-K(1663200.0)*COS(c*x)*c*c*c*c*x*x*x*x+K(6652800.0)*SIN(c*x)*c*c*c*x*x*x+K(19958400.0)*COS(c*x)*c*c*x*x-K(39916800.0)*SIN(c*x)*c*x-K(39916800.0)*COS(c*x))/POW(x,K(12.0)); break;
    case 12 : value=(COS(c*x)*POW(c,K(12.0))*POW(x,K(12.0))-K(12.0)*SIN(c*x)*POW(c,K(11.0))*POW(x,K(11.0))-K(132.0)*COS(c*x)*POW(c,K(10.0))*POW(x,K(10.0))+K(1320.0)*SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+K(11880.0)*COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-K(95040.0)*SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-K(665280.0)*COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+K(3991680.0)*SIN(c*x)*c*c*c*c*c*x*x*x*x*x+K(19958400.0)*COS(c*x)*c*c*c*c*x*x*x*x-K(79833600.0)*SIN(c*x)*c*c*c*x*x*x-K(239500800.0)*COS(c*x)*c*c*x*x+K(479001600.0)*SIN(c*x)*c*x+K(479001600.0)*COS(c*x))/POW(x,K(13.0)); break;
    default : value=K(0.0);
  }
  value *= POW(sign, (R)(der));

  return value;
}

C kcot(R x, int der, const R *param)   /* K(x) = cot(cx) */
{
  R c=param[0];
  R value=K(0.0);

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value = K(1.0)/TAN(c * x); break;
    case  1 : value = -(K(1.0) + POW(K(1.0)/TAN(c * x), K(2.0))) * c; break;
    case  2 : value = K(2.0) * K(1.0)/TAN(c * x) * (K(1.0) + POW(K(1.0)/TAN(c * x), K(2.0))) * c * c; break;
    case  3 : value = -K(2.0) * (K(1.0) + POW(K(1.0)/TAN(c * x), K(2.0))) * POW(c, K(3.0)) * (K(1.0) + K(3.0) * POW(K(1.0)/TAN(c * x), K(2.0))); break;
    case  4 : value = K(8.0) * (K(1.0) + POW(K(1.0)/TAN(c * x), K(2.0))) * POW(c, K(4.0)) * K(1.0)/TAN(c * x) * (K(2.0) + K(3.0) * POW(K(1.0)/TAN(c * x), K(2.0))); break;
    case  5 : value = -K(0.8e1) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.5e1)) * (K(0.15e2) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.15e2) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.2e1)); break;
    case  6 : value = K(0.16e2) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.6e1)) * K(1.0)/TAN(c * x) * (K(0.60e2) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.45e2) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.17e2)); break;
    case  7 : value = -K(0.16e2) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.7e1)) * (K(0.525e3) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.315e3) * POW(K(1.0)/TAN(c * x), K(0.6e1)) + K(0.231e3) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.17e2)); break;
    case  8 : value = K(0.128e3) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.8e1)) * K(1.0)/TAN(c * x) * (K(0.630e3) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.315e3) * POW(K(1.0)/TAN(c * x), K(0.6e1)) + K(0.378e3) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.62e2)); break;
    case  9 : value = -K(0.128e3) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.9e1)) * (K(0.6615e4) * POW(K(1.0)/TAN(c * x), K(0.6e1)) + K(0.2835e4) * POW(K(1.0)/TAN(c * x), K(0.8e1)) + K(0.5040e4) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.1320e4) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.62e2)); break;
    case 10 : value = K(0.256e3) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.10e2)) * K(1.0)/TAN(c * x) * (K(0.37800e5) * POW(K(1.0)/TAN(c * x), K(0.6e1)) + K(0.14175e5) * POW(K(1.0)/TAN(c * x), K(0.8e1)) + K(0.34965e5) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.12720e5) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.1382e4)); break;
    case 11 : value = -K(0.256e3) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.11e2)) * (K(0.467775e6) * POW(K(1.0)/TAN(c * x), K(0.8e1)) + K(0.155925e6) * POW(K(1.0)/TAN(c * x), K(0.10e2)) + K(0.509355e6) * POW(K(1.0)/TAN(c * x), K(0.6e1)) + K(0.238425e6) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.42306e5) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.1382e4)); break;
    case 12 : value = K(0.1024e4) * (K(0.1e1) + POW(K(1.0)/TAN(c * x), K(0.2e1))) * POW(c, K(0.12e2)) * K(1.0)/TAN(c * x) * (K(0.1559250e7) * POW(K(1.0)/TAN(c * x), K(0.8e1)) + K(0.467775e6) * POW(K(1.0)/TAN(c * x), K(0.10e2)) + K(0.1954260e7) * POW(K(1.0)/TAN(c * x), K(0.6e1)) + K(0.1121670e7) * POW(K(1.0)/TAN(c * x), K(0.4e1)) + K(0.280731e6) * POW(K(1.0)/TAN(c * x), K(0.2e1)) + K(0.21844e5)); break;
    default : value=K(0.0);
  }

  return value;
}


C one_over_cube(R x, int der, const R *param)
{
  R value=K(0.0);
  UNUSED(param);

  if (FABS(x)<EPSILON) value=K(0.0);
  else switch (der)
  {
    case  0 : value = K(1.0)/(x*x*x); break;
    case  1 : value = -K(3.0)/(x*x*x*x); break;
    case  2 : value = K(12.0)/(x*x*x*x*x); break;
    case  3 : value = -K(60.0)/(x*x*x*x*x*x); break;
    case  4 : value = K(360.0)/(x*x*x*x*x*x*x); break;
    case  5 : value = -K(2520.0)/(x*x*x*x*x*x*x*x); break;
    case  6 : value = K(20160.0)/POW(x, K(9.0)); break;
    case  7 : value = -K(181440.0)/POW(x, K(10.0)); break;
    case  8 : value = K(1814400.0)/POW(x, K(11.0)); break;
    case  9 : value = -K(19958400.0)/POW(x, K(12.0)); break;
    case  10 : value = K(239500800.0)/POW(x, K(13.0)); break;
    case  11 : value = -K(3113510400.0)/POW(x, K(14.0)); break;
    case  12 : value = K(43589145600.0)/POW(x, K(15.0)); break;
    default : value=K(0.0);
  }

  return value;
}


C log_sin(R x, int der, const R *param)   /* K(x) = log(|sin(cx)|) */
{
  R c=param[0];
  R value=K(0.0);

  if (FABS(x)<EPSILON) value=K(0.0);
  else
  {
      if (der == 0) value = LOG(FABS(SIN(c * x)));
      else value = c * kcot(x, der-1, param);
  }
  
  return value;
}

C laplacian_rbf(R x, int der, const R *param)    /* K(x)=EXP(-|x|/c) */
{
  R c = param[0];
  R value = K(0.0);

  switch (der)
  {
    case  0: value = EXP(-FABS(x)/c); break;
    default:
      value = EXP(-FABS(x)/c)/POW(-c,(R)der);
      if (x < K(0.0))
        value *= POW(K(-1.0),(R)der);
  }

  return value;
}

/* \} */

/* kernels.c */
