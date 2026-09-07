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

NFFT_C gaussian(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=EXP(-x^2/c^2) */
{
  NFFT_R c = param[0];
  NFFT_R value = NFFT_K(0.0);

  switch (der)
  {
    case  0 : value = NFFT_EXP(-x*x/(c*c)); break;
    case  1 : value = -NFFT_K(2.0)*x/(c*c)*NFFT_EXP(-x*x/(c*c)); break;
    case  2 : value = NFFT_K(2.0)*NFFT_EXP(-x*x/(c*c))*(-c*c+NFFT_K(2.0)*x*x)/(c*c*c*c); break;
    case  3 : value = -NFFT_K(4.0)*x*NFFT_EXP(-x*x/(c*c))*(-NFFT_K(3.0)*c*c+NFFT_K(2.0)*x*x)/(c*c*c*c*c*c); break;
    case  4 : value = NFFT_K(4.0)*NFFT_EXP(-x*x/(c*c))*(NFFT_K(3.0)*c*c*c*c-NFFT_K(12.0)*c*c*x*x+NFFT_K(4.0)*x*x*x*x)/(c*c*c*c*c*c*c*c); break;
    case  5 : value = -NFFT_K(8.0)*x*NFFT_EXP(-x*x/(c*c))*(NFFT_K(15.0)*c*c*c*c-NFFT_K(20.0)*c*c*x*x+NFFT_K(4.0)*x*x*x*x)/NFFT_POW(c,NFFT_K(10.0)); break;
    case  6 : value = NFFT_K(8.0)*NFFT_EXP(-x*x/(c*c))*(-NFFT_K(15.0)*c*c*c*c*c*c+NFFT_K(90.0)*x*x*c*c*c*c-NFFT_K(60.0)*x*x*x*x*c*c+NFFT_K(8.0)*x*x*x*x*x*x)/NFFT_POW(c,NFFT_K(12.0)); break;
    case  7 : value = -NFFT_K(16.0)*x*NFFT_EXP(-x*x/(c*c))*(-NFFT_K(105.0)*c*c*c*c*c*c+NFFT_K(210.0)*x*x*c*c*c*c-NFFT_K(84.0)*x*x*x*x*c*c+NFFT_K(8.0)*x*x*x*x*x*x)/NFFT_POW(c,NFFT_K(14.0)); break;
    case  8 : value = NFFT_K(16.0)*NFFT_EXP(-x*x/(c*c))*(NFFT_K(105.0)*c*c*c*c*c*c*c*c-NFFT_K(840.0)*x*x*c*c*c*c*c*c+NFFT_K(840.0)*x*x*x*x*c*c*c*c-NFFT_K(224.0)*x*x*x*x*x*x*c*c+NFFT_K(16.0)*x*x*x*x*x*x*x*x)/NFFT_POW(c,NFFT_K(16.0)); break;
    case  9 : value = -NFFT_K(32.0)*x*NFFT_EXP(-x*x/(c*c))*(NFFT_K(945.0)*c*c*c*c*c*c*c*c-NFFT_K(2520.0)*x*x*c*c*c*c*c*c+NFFT_K(1512.0)*x*x*x*x*c*c*c*c-NFFT_K(288.0)*x*x*x*x*x*x*c*c+NFFT_K(16.0)*x*x*x*x*x*x*x*x)/NFFT_POW(c,NFFT_K(18.0)); break;
    case 10 : value = NFFT_K(32.0)*NFFT_EXP(-x*x/(c*c))*(-NFFT_K(945.0)*NFFT_POW(c,NFFT_K(10.0))+NFFT_K(9450.0)*x*x*c*c*c*c*c*c*c*c-NFFT_K(12600.0)*x*x*x*x*c*c*c*c*c*c+NFFT_K(5040.0)*x*x*x*x*x*x*c*c*c*c-NFFT_K(720.0)*x*x*x*x*x*x*x*x*c*c+NFFT_K(32.0)*NFFT_POW(x,NFFT_K(10.0)))/NFFT_POW(c,NFFT_K(20.0)); break;
    case 11 : value = -NFFT_K(64.0)*x*NFFT_EXP(-x*x/(c*c))*(-NFFT_K(10395.0)*NFFT_POW(c,NFFT_K(10.0))+NFFT_K(34650.0)*x*x*c*c*c*c*c*c*c*c-NFFT_K(27720.0)*x*x*x*x*c*c*c*c*c*c+NFFT_K(7920.0)*x*x*x*x*x*x*c*c*c*c-NFFT_K(880.0)*x*x*x*x*x*x*x*x*c*c+NFFT_K(32.0)*NFFT_POW(x,NFFT_K(10.0)))/NFFT_POW(c,NFFT_K(22.0)); break;
    case 12 : value = NFFT_K(64.0)*NFFT_EXP(-x*x/(c*c))*(NFFT_K(10395.0)*NFFT_POW(c,NFFT_K(12.0))-NFFT_K(124740.0)*x*x*NFFT_POW(c,NFFT_K(10.0))+NFFT_K(207900.0)*x*x*x*x*c*c*c*c*c*c*c*c-NFFT_K(110880.0)*x*x*x*x*x*x*c*c*c*c*c*c+NFFT_K(23760.0)*x*x*x*x*x*x*x*x*c*c*c*c-NFFT_K(2112.0)*NFFT_POW(x,NFFT_K(10.0))*c*c+NFFT_K(64.0)*NFFT_POW(x,NFFT_K(12.0)))/NFFT_POW(c,NFFT_K(24.0)); break;
    default : value = NFFT_K(0.0);
  }

  return value;
}

NFFT_C multiquadric(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=SQRT(x^2+c^2) */
{
  NFFT_R c=param[0];
  NFFT_R value=NFFT_K(0.0);

  switch (der)
  {
    case  0 : value=NFFT_SQRT(x*x+c*c); break;
    case  1 : value=NFFT_K(1.0)/(NFFT_SQRT(x*x+c*c))*x; break;
    case  2 : value=c*c/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(3.0))); break;
    case  3 : value=-NFFT_K(3.0)*x*c*c/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(5.0))); break;
    case  4 : value=NFFT_K(3.0)*c*c*(NFFT_K(4.0)*x*x-c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(7.0))); break;
    case  5 : value=-NFFT_K(15.0)*x*c*c*(NFFT_K(4.0)*x*x-NFFT_K(3.0)*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(9.0))); break;
    case  6 : value=NFFT_K(45.0)*c*c*(NFFT_K(8.0)*x*x*x*x-NFFT_K(12.0)*x*x*c*c+c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(11.0))); break;
    case  7 : value=-NFFT_K(315.0)*x*c*c*(NFFT_K(8.0)*x*x*x*x-NFFT_K(20.0)*x*x*c*c+NFFT_K(5.0)*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(13.0))); break;
    case  8 : value=NFFT_K(315.0)*c*c*(NFFT_K(64.0)*x*x*x*x*x*x-NFFT_K(240.0)*x*x*x*x*c*c+NFFT_K(120.0)*x*x*c*c*c*c-NFFT_K(5.0)*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(15.0))); break;
    case  9 : value=-NFFT_K(2835.0)*x*c*c*(NFFT_K(64.0)*x*x*x*x*x*x-NFFT_K(336.0)*x*x*x*x*c*c+NFFT_K(280.0)*x*x*c*c*c*c-NFFT_K(35.0)*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(17.0))); break;
    case 10 : value=NFFT_K(14175.0)*c*c*(NFFT_K(128.0)*x*x*x*x*x*x*x*x-NFFT_K(896.0)*x*x*x*x*x*x*c*c+NFFT_K(1120.0)*x*x*x*x*c*c*c*c-NFFT_K(280.0)*x*x*c*c*c*c*c*c+NFFT_K(7.0)*c*c*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(19.0))); break;
    case 11 : value=-NFFT_K(155925.0)*x*c*c*(NFFT_K(128.0)*x*x*x*x*x*x*x*x-NFFT_K(1152.0)*x*x*x*x*x*x*c*c+NFFT_K(2016.0)*x*x*x*x*c*c*c*c-NFFT_K(840.0)*x*x*c*c*c*c*c*c+NFFT_K(63.0)*c*c*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(21.0))); break;
    case 12 : value=NFFT_K(467775.0)*c*c*(NFFT_K(1260.0)*x*x*c*c*c*c*c*c*c*c-NFFT_K(21.0)*NFFT_POW(c,NFFT_K(10.0))+NFFT_K(512.0)*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(5760.0)*x*x*x*x*x*x*x*x*c*c+NFFT_K(13440.0)*x*x*x*x*x*x*c*c*c*c-NFFT_K(8400.0)*x*x*x*x*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(23.0))); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C inverse_multiquadric(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=1/SQRT(x^2+c^2) */
{
  NFFT_R c=param[0];
  NFFT_R value=NFFT_K(0.0);

  switch (der)
  {
    case  0 : value=NFFT_K(1.0)/NFFT_SQRT(x*x+c*c); break;
    case  1 : value=-NFFT_K(1.0)/(NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(3.0))))*x; break;
    case  2 : value=(NFFT_K(2.0)*x*x-c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(5.0))); break;
    case  3 : value=-NFFT_K(3.0)*x*(NFFT_K(2.0)*x*x-NFFT_K(3.0)*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(7.0))); break;
    case  4 : value=NFFT_K(3.0)*(NFFT_K(8.0)*x*x*x*x-NFFT_K(24.0)*x*x*c*c+NFFT_K(3.0)*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(9.0))); break;
    case  5 : value=-NFFT_K(15.0)*x*(NFFT_K(8.0)*x*x*x*x-NFFT_K(40.0)*x*x*c*c+NFFT_K(15.0)*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(11.0))); break;
    case  6 : value=NFFT_K(45.0)*(NFFT_K(16.0)*x*x*x*x*x*x-NFFT_K(120.0)*x*x*x*x*c*c+NFFT_K(90.0)*x*x*c*c*c*c-NFFT_K(5.0)*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(13.0))); break;
    case  7 : value=-NFFT_K(315.0)*x*(NFFT_K(16.0)*x*x*x*x*x*x-NFFT_K(168.0)*x*x*x*x*c*c+NFFT_K(210.0)*x*x*c*c*c*c-NFFT_K(35.0)*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(15.0))); break;
    case  8 : value=NFFT_K(315.0)*(NFFT_K(128.0)*x*x*x*x*x*x*x*x-NFFT_K(1792.0)*x*x*x*x*x*x*c*c+NFFT_K(3360.0)*x*x*x*x*c*c*c*c-NFFT_K(1120.0)*x*x*c*c*c*c*c*c+NFFT_K(35.0)*c*c*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(17.0))); break;
    case  9 : value=-NFFT_K(2835.0)*x*(NFFT_K(128.0)*x*x*x*x*x*x*x*x-NFFT_K(2304.0)*x*x*x*x*x*x*c*c+NFFT_K(6048.0)*x*x*x*x*c*c*c*c-NFFT_K(3360.0)*x*x*c*c*c*c*c*c+NFFT_K(315.0)*c*c*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(19.0))); break;
    case 10 : value=NFFT_K(14175.0)*(NFFT_K(256.0)*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(5760.0)*x*x*x*x*x*x*x*x*c*c+NFFT_K(20160.0)*x*x*x*x*x*x*c*c*c*c-NFFT_K(16800.0)*x*x*x*x*c*c*c*c*c*c+NFFT_K(3150.0)*x*x*c*c*c*c*c*c*c*c-NFFT_K(63.0)*NFFT_POW(c,NFFT_K(10.0)))/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(21.0))); break;
    case 11 : value=-NFFT_K(155925.0)*x*(NFFT_K(256.0)*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(7040.0)*x*x*x*x*x*x*x*x*c*c+NFFT_K(31680.0)*x*x*x*x*x*x*c*c*c*c-NFFT_K(36960.0)*x*x*x*x*c*c*c*c*c*c+NFFT_K(11550.0)*x*x*c*c*c*c*c*c*c*c-NFFT_K(693.0)*NFFT_POW(c,NFFT_K(10.0)))/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(23.0))); break;
    case 12 : value=NFFT_K(467775.0)*(NFFT_K(231.0)*NFFT_POW(c,NFFT_K(12.0))+NFFT_K(190080.0)*x*x*x*x*x*x*x*x*c*c*c*c-NFFT_K(16632.0)*x*x*NFFT_POW(c,NFFT_K(10.0))-NFFT_K(295680.0)*x*x*x*x*x*x*c*c*c*c*c*c+NFFT_K(138600.0)*x*x*x*x*c*c*c*c*c*c*c*c+NFFT_K(1024.0)*NFFT_POW(x,NFFT_K(12.0))-NFFT_K(33792.0)*NFFT_POW(x,NFFT_K(10.0))*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(25.0))); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C logarithm(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=LOG |x| */
{
  NFFT_R value=NFFT_K(0.0);

  (void)param;

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value=NFFT_LOG(NFFT_FABS(x)); break;
    case  1 : value=(x<0 ? -1 : 1)/NFFT_FABS(x); break;
    case  2 : value=-1/(x*x); break;
    case  3 : value=NFFT_K(2.0)*(x<0 ? -1 : 1)/NFFT_POW(NFFT_FABS(x),NFFT_K(3.0)); break;
    case  4 : value=-NFFT_K(6.0)/(x*x*x*x); break;
    case  5 : value=NFFT_K(24.0)*(x<0 ? -1 : 1)/NFFT_POW(NFFT_FABS(x),NFFT_K(5.0)); break;
    case  6 : value=-NFFT_K(120.0)/(x*x*x*x*x*x); break;
    case  7 : value=NFFT_K(720.0)*(x<0 ? -1 : 1)/NFFT_POW(NFFT_FABS(x),NFFT_K(7.0)); break;
    case  8 : value=-NFFT_K(5040.0)/(x*x*x*x*x*x*x*x); break;
    case  9 : value=NFFT_K(40320.0)*(x<0 ? -1 : 1)/NFFT_POW(NFFT_FABS(x),NFFT_K(9.0)); break;
    case 10 : value=-NFFT_K(362880.0)/NFFT_POW(x,NFFT_K(10.0)); break;
    case 11 : value=NFFT_K(3628800.0)*(x<0 ? -1 : 1)/NFFT_POW(NFFT_FABS(x),NFFT_K(11.0)); break;
    case 12 : value=-NFFT_K(39916800.0)/NFFT_POW(x,NFFT_K(12.0)); break;
    case 13 : value=NFFT_K(479001600.0)/NFFT_POW(x,NFFT_K(13.0)); break;
    case 14 : value=-NFFT_K(6227020800.0)/NFFT_POW(x,NFFT_K(14.0)); break;
    case 15 : value=NFFT_K(87178291200.0)/NFFT_POW(x,NFFT_K(15.0)); break;
    case 16 : value=-NFFT_K(1307674368000.0)/NFFT_POW(x,NFFT_K(16.0)); break;
    case 17 : value=NFFT_K(20922789888000.0)/NFFT_POW(x,NFFT_K(17.0)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C thinplate_spline(NFFT_R x, int der, const NFFT_R *param)    /* K(x) = x^2 LOG |x| */
{
  NFFT_R value=NFFT_K(0.0);

  (void)param;

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value=x*x*NFFT_LOG(NFFT_FABS(x)); break;
    case  1 : value=NFFT_K(2.0)*x*NFFT_LOG(NFFT_FABS(x))+x; break;
    case  2 : value=NFFT_K(2.0)*NFFT_LOG(NFFT_FABS(x))+NFFT_K(3.0); break;
    case  3 : value=NFFT_K(2.0)/x; break;
    case  4 : value=-NFFT_K(2.0)/(x*x); break;
    case  5 : value=NFFT_K(4.0)/(x*x*x); break;
    case  6 : value=-NFFT_K(12.0)/(x*x*x*x); break;
    case  7 : value=NFFT_K(48.0)/(x*x*x*x*x); break;
    case  8 : value=-NFFT_K(240.0)/(x*x*x*x*x*x); break;
    case  9 : value=NFFT_K(1440.0)/(x*x*x*x*x*x*x); break;
    case 10 : value=-NFFT_K(10080.0)/(x*x*x*x*x*x*x*x); break;
    case 11 : value=NFFT_K(80640.0)/(x*x*x*x*x*x*x*x*x); break;
    case 12 : value=-NFFT_K(725760.0)/NFFT_POW(x,NFFT_K(10.0)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C one_over_square(NFFT_R x, int der, const NFFT_R *param)    /* K(x) = 1/x^2 */
{
  NFFT_R value=NFFT_K(0.0);

  (void)param;

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value=NFFT_K(1.0)/(x*x); break;
    case  1 : value=-NFFT_K(2.0)/(x*x*x); break;
    case  2 : value=NFFT_K(6.0)/(x*x*x*x); break;
    case  3 : value=-NFFT_K(24.0)/(x*x*x*x*x); break;
    case  4 : value=NFFT_K(120.0)/(x*x*x*x*x*x); break;
    case  5 : value=-NFFT_K(720.0)/(x*x*x*x*x*x*x); break;
    case  6 : value=NFFT_K(5040.0)/(x*x*x*x*x*x*x*x); break;
    case  7 : value=-NFFT_K(40320.0)/(x*x*x*x*x*x*x*x*x); break;
    case  8 : value=NFFT_K(362880.0)/NFFT_POW(x,NFFT_K(10.0)); break;
    case  9 : value=-NFFT_K(3628800.0)/NFFT_POW(x,NFFT_K(11.0)); break;
    case 10 : value=NFFT_K(39916800.0)/NFFT_POW(x,NFFT_K(12.0)); break;
    case 11 : value=-NFFT_K(479001600.0)/NFFT_POW(x,NFFT_K(13.0)); break;
    case 12 : value=NFFT_K(6227020800.0)/NFFT_POW(x,NFFT_K(14.0)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C one_over_modulus(NFFT_R x, int der, const NFFT_R *param)    /* K(x) = 1/|x| */
{
  NFFT_R value=NFFT_K(0.0);

  (void)param;

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value=NFFT_K(1.0)/NFFT_FABS(x); break;
    case  1 : value=-1/x/NFFT_FABS(x); break;
    case  2 : value=NFFT_K(2.0)/NFFT_POW(NFFT_FABS(x),NFFT_K(3.0)); break;
    case  3 : value=-NFFT_K(6.0)/(x*x*x)/NFFT_FABS(x); break;
    case  4 : value=NFFT_K(24.0)/NFFT_POW(NFFT_FABS(x),NFFT_K(5.0)); break;
    case  5 : value=-NFFT_K(120.0)/(x*x*x*x*x)/NFFT_FABS(x); break;
    case  6 : value=NFFT_K(720.0)/NFFT_POW(NFFT_FABS(x),NFFT_K(7.0)); break;
    case  7 : value=-NFFT_K(5040.0)/(x*x*x*x*x*x*x)/NFFT_FABS(x); break;
    case  8 : value=NFFT_K(40320.0)/NFFT_POW(NFFT_FABS(x),NFFT_K(9.0)); break;
    case  9 : value=-NFFT_K(362880.0)/(x*x*x*x*x*x*x*x*x)/NFFT_FABS(x); break;
    case 10 : value=NFFT_K(3628800.0)/NFFT_POW(NFFT_FABS(x),NFFT_K(11.0)); break;
    case 11 : value=-NFFT_K(39916800.0)/NFFT_POW(x,NFFT_K(11.0))/NFFT_FABS(x); break;
    case 12 : value=NFFT_K(479001600.0)/NFFT_POW(NFFT_FABS(x),NFFT_K(13.0)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C one_over_x(NFFT_R x, int der, const NFFT_R *param)    /* K(x) = 1/x */
{
  NFFT_R value=NFFT_K(0.0);

  (void)param;

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value=NFFT_K(1.0)/x; break;
    case  1 : value=-NFFT_K(1.0)/(x*x); break;
    case  2 : value=NFFT_K(2.0)/(x*x*x); break;
    case  3 : value=-NFFT_K(6.0)/(x*x*x*x); break;
    case  4 : value=NFFT_K(24.0)/(x*x*x*x*x); break;
    case  5 : value=-NFFT_K(120.0)/(x*x*x*x*x*x); break;
    case  6 : value=NFFT_K(720.0)/(x*x*x*x*x*x*x); break;
    case  7 : value=-NFFT_K(5040.0)/(x*x*x*x*x*x*x*x); break;
    case  8 : value=NFFT_K(40320.0)/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-NFFT_K(362880.0)/NFFT_POW(x,NFFT_K(10.0)); break;
    case 10 : value=NFFT_K(3628800.0)/NFFT_POW(x,NFFT_K(11.0)); break;
    case 11 : value=-NFFT_K(39916800.0)/NFFT_POW(x,NFFT_K(12.0)); break;
    case 12 : value=NFFT_K(479001600.0)/NFFT_POW(x,NFFT_K(13.0)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C inverse_multiquadric3(NFFT_R x, int der, const NFFT_R *param)    /* K(x) = 1/SQRT(x^2+c^2)^3 */
{
  NFFT_R c=param[0];
  NFFT_R value=NFFT_K(0.0);

  switch (der)
  {
    case  0 : value=NFFT_K(1.0)/(NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(3.0)))); break;
    case  1 : value=-NFFT_K(3.0)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(5.0)))*x; break;
    case  2 : value=NFFT_K(3.0)*(NFFT_K(4.0)*x*x-c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(7.0))); break;
    case  3 : value=-NFFT_K(15.0)*x*(NFFT_K(4.0)*x*x-NFFT_K(3.0)*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(9.0))); break;
    case  4 : value=NFFT_K(45.0)*(NFFT_K(8.0)*x*x*x*x-NFFT_K(12.0)*x*x*c*c+c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(11.0))); break;
    case  5 : value=-NFFT_K(315.0)*x*(NFFT_K(8.0)*x*x*x*x-NFFT_K(20.0)*x*x*c*c+NFFT_K(5.0)*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(13.0))); break;
    case  6 : value=NFFT_K(315.0)*(NFFT_K(64.0)*x*x*x*x*x*x-NFFT_K(240.0)*x*x*x*x*c*c+NFFT_K(120.0)*x*x*c*c*c*c-NFFT_K(5.0)*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(15.0))); break;
    case  7 : value=-NFFT_K(2835.0)*x*(NFFT_K(64.0)*x*x*x*x*x*x-NFFT_K(336.0)*x*x*x*x*c*c+NFFT_K(280.0)*x*x*c*c*c*c-NFFT_K(35.0)*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(17.0))); break;
    case  8 : value=NFFT_K(14175.0)*(NFFT_K(128.0)*x*x*x*x*x*x*x*x-NFFT_K(896.0)*x*x*x*x*x*x*c*c+NFFT_K(1120.0)*x*x*x*x*c*c*c*c-NFFT_K(280.0)*x*x*c*c*c*c*c*c+NFFT_K(7.0)*c*c*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(19.0))); break;
    case  9 : value=-NFFT_K(155925.0)*x*(NFFT_K(128.0)*x*x*x*x*x*x*x*x-NFFT_K(1152.0)*x*x*x*x*x*x*c*c+NFFT_K(2016.0)*x*x*x*x*c*c*c*c-NFFT_K(840.0)*x*x*c*c*c*c*c*c+NFFT_K(63.0)*c*c*c*c*c*c*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(21.0))); break;
    case 10 : value=NFFT_K(467775.0)*(NFFT_K(512.0)*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(5760.0)*x*x*x*x*x*x*x*x*c*c+NFFT_K(13440.0)*x*x*x*x*x*x*c*c*c*c-NFFT_K(8400.0)*x*x*x*x*c*c*c*c*c*c+NFFT_K(1260.0)*x*x*c*c*c*c*c*c*c*c-NFFT_K(21.0)*NFFT_POW(c,NFFT_K(10.0)))/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(23.0))); break;
    case 11 : value=-NFFT_K(6081075.0)*x*(NFFT_K(512.0)*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(7040.0)*x*x*x*x*x*x*x*x*c*c+NFFT_K(21120.0)*x*x*x*x*x*x*c*c*c*c-NFFT_K(18480.0)*x*x*x*x*c*c*c*c*c*c+NFFT_K(4620.0)*x*x*c*c*c*c*c*c*c*c-NFFT_K(231.0)*NFFT_POW(c,NFFT_K(10.0)))/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(25.0))); break;
    case 12 : value=NFFT_K(42567525.0)*(NFFT_K(1024.0)*NFFT_POW(x,NFFT_K(12.0))+NFFT_K(27720.0)*x*x*x*x*c*c*c*c*c*c*c*c+NFFT_K(33.0)*NFFT_POW(c,NFFT_K(12.0))-NFFT_K(2772.0)*x*x*NFFT_POW(c,NFFT_K(10.0))-NFFT_K(73920.0)*x*x*x*x*x*x*c*c*c*c*c*c+NFFT_K(63360.0)*x*x*x*x*x*x*x*x*c*c*c*c-NFFT_K(16896.0)*NFFT_POW(x,NFFT_K(10.0))*c*c)/NFFT_SQRT(NFFT_POW(x*x+c*c,NFFT_K(27.0))); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C sinc_kernel(NFFT_R x, int der, const NFFT_R *param)    /* K(x) = SIN(cx)/x */
{
  NFFT_R c=param[0];
  NFFT_R value=NFFT_K(0.0);

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=c;
  else switch (der)
  {
    case  0 : value=NFFT_SIN(c*x)/x; break;
    case  1 : value=(NFFT_COS(c*x)*c*x-NFFT_SIN(c*x))/(x*x); break;
    case  2 : value=-(NFFT_SIN(c*x)*c*c*x*x+NFFT_K(2.0)*NFFT_COS(c*x)*c*x-NFFT_K(2.0)*NFFT_SIN(c*x))/(x*x*x); break;
    case  3 : value=-(NFFT_COS(c*x)*c*c*c*x*x*x-NFFT_K(3.0)*NFFT_SIN(c*x)*c*c*x*x-NFFT_K(6.0)*NFFT_COS(c*x)*c*x+NFFT_K(6.0)*NFFT_SIN(c*x))/(x*x*x*x); break;
    case  4 : value=(NFFT_SIN(c*x)*c*c*c*c*x*x*x*x+NFFT_K(4.0)*NFFT_COS(c*x)*c*c*c*x*x*x-NFFT_K(12.0)*NFFT_SIN(c*x)*c*c*x*x-NFFT_K(24.0)*NFFT_COS(c*x)*c*x+NFFT_K(24.0)*NFFT_SIN(c*x))/(x*x*x*x*x); break;
    case  5 : value=(NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(5.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x-NFFT_K(20.0)*NFFT_COS(c*x)*c*c*c*x*x*x+NFFT_K(60.0)*NFFT_SIN(c*x)*c*c*x*x+NFFT_K(120.0)*NFFT_COS(c*x)*c*x-NFFT_K(120.0)*NFFT_SIN(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(NFFT_SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+NFFT_K(6.0)*NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(30.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x-NFFT_K(120.0)*NFFT_COS(c*x)*c*c*c*x*x*x+NFFT_K(360.0)*NFFT_SIN(c*x)*c*c*x*x+NFFT_K(720.0)*NFFT_COS(c*x)*c*x-NFFT_K(720.0)*NFFT_SIN(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=-(NFFT_COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-NFFT_K(7.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(42.0)*NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(210.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x+NFFT_K(840.0)*NFFT_COS(c*x)*c*c*c*x*x*x-NFFT_K(2520.0)*NFFT_SIN(c*x)*c*c*x*x-NFFT_K(5040.0)*NFFT_COS(c*x)*c*x+NFFT_K(5040.0)*NFFT_SIN(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+NFFT_K(8.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-NFFT_K(56.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(336.0)*NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(1680.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x+NFFT_K(6720.0)*NFFT_COS(c*x)*c*c*c*x*x*x-NFFT_K(20160.0)*NFFT_SIN(c*x)*c*c*x*x-NFFT_K(40320.0)*NFFT_COS(c*x)*c*x+NFFT_K(40320.0)*NFFT_SIN(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=(NFFT_COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-NFFT_K(9.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-NFFT_K(72.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+NFFT_K(504.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+NFFT_K(3024.0)*NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(15120.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x-NFFT_K(60480.0)*NFFT_COS(c*x)*c*c*c*x*x*x+NFFT_K(181440.0)*NFFT_SIN(c*x)*c*c*x*x+NFFT_K(362880.0)*NFFT_COS(c*x)*c*x-NFFT_K(362880.0)*NFFT_SIN(c*x))/NFFT_POW(x,NFFT_K(10.0)); break;
    case 10 : value=-(NFFT_SIN(c*x)*NFFT_POW(c,NFFT_K(10.0))*NFFT_POW(x,NFFT_K(10.0))+NFFT_K(10.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-NFFT_K(90.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-NFFT_K(720.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+NFFT_K(5040.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+NFFT_K(30240.0)*NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(151200.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x-NFFT_K(604800.0)*NFFT_COS(c*x)*c*c*c*x*x*x+NFFT_K(1814400.0)*NFFT_SIN(c*x)*c*c*x*x+NFFT_K(3628800.0)*NFFT_COS(c*x)*c*x-NFFT_K(3628800.0)*NFFT_SIN(c*x))/NFFT_POW(x,NFFT_K(11.0)); break;
    case 11 : value=-(NFFT_COS(c*x)*NFFT_POW(c,NFFT_K(11.0))*NFFT_POW(x,NFFT_K(11.0))-NFFT_K(11.0)*NFFT_SIN(c*x)*NFFT_POW(c,NFFT_K(10.0))*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(110.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+NFFT_K(990.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+NFFT_K(7920.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-NFFT_K(55440.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(332640.0)*NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(1663200.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x+NFFT_K(6652800.0)*NFFT_COS(c*x)*c*c*c*x*x*x-NFFT_K(19958400.0)*NFFT_SIN(c*x)*c*c*x*x-NFFT_K(39916800.0)*NFFT_COS(c*x)*c*x+NFFT_K(39916800.0)*NFFT_SIN(c*x))/NFFT_POW(x,NFFT_K(12.0)); break;
    case 12 : value=(NFFT_SIN(c*x)*NFFT_POW(c,NFFT_K(12.0))*NFFT_POW(x,NFFT_K(12.0))+NFFT_K(12.0)*NFFT_COS(c*x)*NFFT_POW(c,NFFT_K(11.0))*NFFT_POW(x,NFFT_K(11.0))-NFFT_K(132.0)*NFFT_SIN(c*x)*NFFT_POW(c,NFFT_K(10.0))*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(1320.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+NFFT_K(11880.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+NFFT_K(95040.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-NFFT_K(665280.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(3991680.0)*NFFT_COS(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(19958400.0)*NFFT_SIN(c*x)*c*c*c*c*x*x*x*x+NFFT_K(79833600.0)*NFFT_COS(c*x)*c*c*c*x*x*x-NFFT_K(239500800.0)*NFFT_SIN(c*x)*c*c*x*x-NFFT_K(479001600.0)*NFFT_COS(c*x)*c*x+NFFT_K(479001600.0)*NFFT_SIN(c*x))/NFFT_POW(x,NFFT_K(13.0)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}

NFFT_C cosc(NFFT_R x, int der, const NFFT_R *param)    /* K(x) = COS(cx)/x */
{
  NFFT_R c=param[0];
  NFFT_R value=NFFT_K(0.0);
  NFFT_R sign;

  if (x<0) sign=-NFFT_K(1.0); else sign=NFFT_K(1.0);
  x=NFFT_FABS(x);

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value=NFFT_COS(c*x)/x; break;
    case  1 : value=-(NFFT_SIN(c*x)*c*x+NFFT_COS(c*x))/(x*x); break;
    case  2 : value=(-NFFT_COS(c*x)*c*c*x*x+NFFT_K(2.0)*NFFT_SIN(c*x)*c*x+NFFT_K(2.0)*NFFT_COS(c*x))/(x*x*x); break;
    case  3 : value=(NFFT_SIN(c*x)*c*c*c*x*x*x+NFFT_K(3.0)*NFFT_COS(c*x)*c*c*x*x-NFFT_K(6.0)*NFFT_SIN(c*x)*c*x-NFFT_K(6.0)*NFFT_COS(c*x))/(x*x*x*x); break;
    case  4 : value=(NFFT_COS(c*x)*c*c*c*c*x*x*x*x-NFFT_K(4.0)*NFFT_SIN(c*x)*c*c*c*x*x*x-NFFT_K(12.0)*NFFT_COS(c*x)*c*c*x*x+NFFT_K(24.0)*NFFT_SIN(c*x)*c*x+NFFT_K(24.0)*NFFT_COS(c*x))/(x*x*x*x*x); break;
    case  5 : value=-(NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(5.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x-NFFT_K(20.0)*NFFT_SIN(c*x)*c*c*c*x*x*x-NFFT_K(60.0)*NFFT_COS(c*x)*c*c*x*x+NFFT_K(120.0)*NFFT_SIN(c*x)*c*x+NFFT_K(120.0)*NFFT_COS(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(NFFT_COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(6.0)*NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(30.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x+NFFT_K(120.0)*NFFT_SIN(c*x)*c*c*c*x*x*x+NFFT_K(360.0)*NFFT_COS(c*x)*c*c*x*x-NFFT_K(720.0)*NFFT_SIN(c*x)*c*x-NFFT_K(720.0)*NFFT_COS(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=(NFFT_SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+NFFT_K(7.0)*NFFT_COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(42.0)*NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(210.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x+NFFT_K(840.0)*NFFT_SIN(c*x)*c*c*c*x*x*x+NFFT_K(2520.0)*NFFT_COS(c*x)*c*c*x*x-NFFT_K(5040.0)*NFFT_SIN(c*x)*c*x-NFFT_K(5040.0)*NFFT_COS(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(NFFT_COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-NFFT_K(8.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-NFFT_K(56.0)*NFFT_COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+NFFT_K(336.0)*NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(1680.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x-NFFT_K(6720.0)*NFFT_SIN(c*x)*c*c*c*x*x*x-NFFT_K(20160.0)*NFFT_COS(c*x)*c*c*x*x+NFFT_K(40320.0)*NFFT_SIN(c*x)*c*x+NFFT_K(40320.0)*NFFT_COS(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-(NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+NFFT_K(9.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-NFFT_K(72.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-NFFT_K(504.0)*NFFT_COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+NFFT_K(3024.0)*NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(15120.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x-NFFT_K(60480.0)*NFFT_SIN(c*x)*c*c*c*x*x*x-NFFT_K(181440.0)*NFFT_COS(c*x)*c*c*x*x+NFFT_K(362880.0)*NFFT_SIN(c*x)*c*x+NFFT_K(362880.0)*NFFT_COS(c*x))/NFFT_POW(x,NFFT_K(10.0)); break;
    case 10 : value=-(NFFT_COS(c*x)*NFFT_POW(c,NFFT_K(10.0))*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(10.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-NFFT_K(90.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+NFFT_K(720.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+NFFT_K(5040.0)*NFFT_COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(30240.0)*NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(151200.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x+NFFT_K(604800.0)*NFFT_SIN(c*x)*c*c*c*x*x*x+NFFT_K(1814400.0)*NFFT_COS(c*x)*c*c*x*x-NFFT_K(3628800.0)*NFFT_SIN(c*x)*c*x-NFFT_K(3628800.0)*NFFT_COS(c*x))/NFFT_POW(x,NFFT_K(11.0)); break;
    case 11 : value=(NFFT_SIN(c*x)*NFFT_POW(c,NFFT_K(11.0))*NFFT_POW(x,NFFT_K(11.0))+NFFT_K(11.0)*NFFT_COS(c*x)*NFFT_POW(c,NFFT_K(10.0))*NFFT_POW(x,NFFT_K(10.0))-NFFT_K(110.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-NFFT_K(990.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+NFFT_K(7920.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+NFFT_K(55440.0)*NFFT_COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-NFFT_K(332640.0)*NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x-NFFT_K(1663200.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x+NFFT_K(6652800.0)*NFFT_SIN(c*x)*c*c*c*x*x*x+NFFT_K(19958400.0)*NFFT_COS(c*x)*c*c*x*x-NFFT_K(39916800.0)*NFFT_SIN(c*x)*c*x-NFFT_K(39916800.0)*NFFT_COS(c*x))/NFFT_POW(x,NFFT_K(12.0)); break;
    case 12 : value=(NFFT_COS(c*x)*NFFT_POW(c,NFFT_K(12.0))*NFFT_POW(x,NFFT_K(12.0))-NFFT_K(12.0)*NFFT_SIN(c*x)*NFFT_POW(c,NFFT_K(11.0))*NFFT_POW(x,NFFT_K(11.0))-NFFT_K(132.0)*NFFT_COS(c*x)*NFFT_POW(c,NFFT_K(10.0))*NFFT_POW(x,NFFT_K(10.0))+NFFT_K(1320.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+NFFT_K(11880.0)*NFFT_COS(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-NFFT_K(95040.0)*NFFT_SIN(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-NFFT_K(665280.0)*NFFT_COS(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+NFFT_K(3991680.0)*NFFT_SIN(c*x)*c*c*c*c*c*x*x*x*x*x+NFFT_K(19958400.0)*NFFT_COS(c*x)*c*c*c*c*x*x*x*x-NFFT_K(79833600.0)*NFFT_SIN(c*x)*c*c*c*x*x*x-NFFT_K(239500800.0)*NFFT_COS(c*x)*c*c*x*x+NFFT_K(479001600.0)*NFFT_SIN(c*x)*c*x+NFFT_K(479001600.0)*NFFT_COS(c*x))/NFFT_POW(x,NFFT_K(13.0)); break;
    default : value=NFFT_K(0.0);
  }
  value *= NFFT_POW(sign, (NFFT_R)(der));

  return value;
}

NFFT_C kcot(NFFT_R x, int der, const NFFT_R *param)   /* K(x) = cot(cx) */
{
  NFFT_R c=param[0];
  NFFT_R value=NFFT_K(0.0);

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value = NFFT_K(1.0)/NFFT_TAN(c * x); break;
    case  1 : value = -(NFFT_K(1.0) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(2.0))) * c; break;
    case  2 : value = NFFT_K(2.0) * NFFT_K(1.0)/NFFT_TAN(c * x) * (NFFT_K(1.0) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(2.0))) * c * c; break;
    case  3 : value = -NFFT_K(2.0) * (NFFT_K(1.0) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(2.0))) * NFFT_POW(c, NFFT_K(3.0)) * (NFFT_K(1.0) + NFFT_K(3.0) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(2.0))); break;
    case  4 : value = NFFT_K(8.0) * (NFFT_K(1.0) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(2.0))) * NFFT_POW(c, NFFT_K(4.0)) * NFFT_K(1.0)/NFFT_TAN(c * x) * (NFFT_K(2.0) + NFFT_K(3.0) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(2.0))); break;
    case  5 : value = -NFFT_K(0.8e1) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.5e1)) * (NFFT_K(0.15e2) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.15e2) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.2e1)); break;
    case  6 : value = NFFT_K(0.16e2) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.6e1)) * NFFT_K(1.0)/NFFT_TAN(c * x) * (NFFT_K(0.60e2) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.45e2) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.17e2)); break;
    case  7 : value = -NFFT_K(0.16e2) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.7e1)) * (NFFT_K(0.525e3) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.315e3) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.6e1)) + NFFT_K(0.231e3) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.17e2)); break;
    case  8 : value = NFFT_K(0.128e3) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.8e1)) * NFFT_K(1.0)/NFFT_TAN(c * x) * (NFFT_K(0.630e3) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.315e3) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.6e1)) + NFFT_K(0.378e3) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.62e2)); break;
    case  9 : value = -NFFT_K(0.128e3) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.9e1)) * (NFFT_K(0.6615e4) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.6e1)) + NFFT_K(0.2835e4) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.8e1)) + NFFT_K(0.5040e4) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.1320e4) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.62e2)); break;
    case 10 : value = NFFT_K(0.256e3) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.10e2)) * NFFT_K(1.0)/NFFT_TAN(c * x) * (NFFT_K(0.37800e5) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.6e1)) + NFFT_K(0.14175e5) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.8e1)) + NFFT_K(0.34965e5) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.12720e5) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.1382e4)); break;
    case 11 : value = -NFFT_K(0.256e3) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.11e2)) * (NFFT_K(0.467775e6) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.8e1)) + NFFT_K(0.155925e6) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.10e2)) + NFFT_K(0.509355e6) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.6e1)) + NFFT_K(0.238425e6) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.42306e5) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.1382e4)); break;
    case 12 : value = NFFT_K(0.1024e4) * (NFFT_K(0.1e1) + NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1))) * NFFT_POW(c, NFFT_K(0.12e2)) * NFFT_K(1.0)/NFFT_TAN(c * x) * (NFFT_K(0.1559250e7) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.8e1)) + NFFT_K(0.467775e6) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.10e2)) + NFFT_K(0.1954260e7) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.6e1)) + NFFT_K(0.1121670e7) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.4e1)) + NFFT_K(0.280731e6) * NFFT_POW(NFFT_K(1.0)/NFFT_TAN(c * x), NFFT_K(0.2e1)) + NFFT_K(0.21844e5)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}


NFFT_C one_over_cube(NFFT_R x, int der, const NFFT_R *param)
{
  NFFT_R value=NFFT_K(0.0);
  NFFT_UNUSED(param);

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else switch (der)
  {
    case  0 : value = NFFT_K(1.0)/(x*x*x); break;
    case  1 : value = -NFFT_K(3.0)/(x*x*x*x); break;
    case  2 : value = NFFT_K(12.0)/(x*x*x*x*x); break;
    case  3 : value = -NFFT_K(60.0)/(x*x*x*x*x*x); break;
    case  4 : value = NFFT_K(360.0)/(x*x*x*x*x*x*x); break;
    case  5 : value = -NFFT_K(2520.0)/(x*x*x*x*x*x*x*x); break;
    case  6 : value = NFFT_K(20160.0)/NFFT_POW(x, NFFT_K(9.0)); break;
    case  7 : value = -NFFT_K(181440.0)/NFFT_POW(x, NFFT_K(10.0)); break;
    case  8 : value = NFFT_K(1814400.0)/NFFT_POW(x, NFFT_K(11.0)); break;
    case  9 : value = -NFFT_K(19958400.0)/NFFT_POW(x, NFFT_K(12.0)); break;
    case  10 : value = NFFT_K(239500800.0)/NFFT_POW(x, NFFT_K(13.0)); break;
    case  11 : value = -NFFT_K(3113510400.0)/NFFT_POW(x, NFFT_K(14.0)); break;
    case  12 : value = NFFT_K(43589145600.0)/NFFT_POW(x, NFFT_K(15.0)); break;
    default : value=NFFT_K(0.0);
  }

  return value;
}


NFFT_C log_sin(NFFT_R x, int der, const NFFT_R *param)   /* K(x) = log(|sin(cx)|) */
{
  NFFT_R c=param[0];
  NFFT_R value=NFFT_K(0.0);

  if (NFFT_FABS(x)<NFFT_R_EPSILON) value=NFFT_K(0.0);
  else
  {
      if (der == 0) value = NFFT_LOG(NFFT_FABS(NFFT_SIN(c * x)));
      else value = c * kcot(x, der-1, param);
  }
  
  return value;
}

NFFT_C laplacian_rbf(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=EXP(-|x|/c) */
{
  NFFT_R c = param[0];
  NFFT_R value = NFFT_K(0.0);

  switch (der)
  {
    case  0: value = NFFT_EXP(-NFFT_FABS(x)/c); break;
    default:
      value = NFFT_EXP(-NFFT_FABS(x)/c)/NFFT_POW(-c,(NFFT_R)der);
      if (x < NFFT_K(0.0))
        value *= NFFT_POW(NFFT_K(-1.0),(NFFT_R)der);
  }

  return value;
}

NFFT_C der_laplacian_rbf(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=|x|/c EXP(-|x|/c) */
{
  NFFT_R c = param[0];
  NFFT_R value = NFFT_K(0.0);

  switch (der)
  {
    case  0 : value = (NFFT_FABS(x)/c)*NFFT_EXP(-NFFT_FABS(x)/c); break;
    default:
        value = (NFFT_POW(NFFT_K(-1.0),(NFFT_R)der))*((NFFT_FABS(x)-(NFFT_R)der*c)/NFFT_POW(c,(NFFT_R)der+1))*NFFT_EXP(-NFFT_FABS(x)/c);
        value *= 1 - 2 * ((x < NFFT_K(0.0)) && (der % 2));
  }

  return value;
}

NFFT_C xx_gaussian(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=x^2/c^2 EXP(-x^2/c^2) */
{
  NFFT_R c = param[0];
  NFFT_R value = NFFT_K(0.0);

  switch (der)
  {
    case  0 : value = x*x*NFFT_EXP(-x*x/(c*c)); break;
    case  1 : value = NFFT_K(2.0)*x*NFFT_EXP(-x*x/(c*c))-NFFT_K(2.0)*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c); break;
    case  2 : value = NFFT_K(4.0)*x*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c)-NFFT_K(10.0)*x*x*NFFT_EXP(-x*x/(c*c))/(c*c)+NFFT_K(2.0)*NFFT_EXP(-x*x/(c*c)); break;
    case  3 : value = -NFFT_K(8.0)*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c)+NFFT_K(36.0)*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c)-NFFT_K(24.0)*x*NFFT_EXP(-x*x/(c*c))/(c*c); break;
    case  4 : value = NFFT_K(16.0)*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c*c*c)-NFFT_K(112.0)*x*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c)+NFFT_K(156.0)*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c)-NFFT_K(24.0)*NFFT_EXP(-x*x/(c*c))/(c*c); break;
    case  5 : value = -NFFT_K(32.0)*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0))+NFFT_K(320.0)*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c*c*c)-NFFT_K(760.0)*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c)+NFFT_K(360.0)*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c); break;
    case  6 : value = NFFT_K(64.0)*x*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(12.0))-NFFT_K(864.0)*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0))+NFFT_K(3120.0)*x*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c*c*c)-NFFT_K(3000.0)*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c)+NFFT_K(360.0)*NFFT_EXP(-x*x/(c*c))/(c*c*c*c); break;
    case  7 : value = -NFFT_K(128.0)*x*x*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(14.0))+NFFT_K(2240.0)*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(12.0))-NFFT_K(11424.0)*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0))+NFFT_K(18480.0)*x*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c*c*c)-NFFT_K(6720.0)*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c); break;
    case  8 : value = NFFT_K(256.0)*NFFT_POW(x,NFFT_K(10.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(16.0))-NFFT_K(5632.0)*x*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(14.0))+NFFT_K(38528.0)*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(12.0))-NFFT_K(94080.0)*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0))+NFFT_K(68880.0)*x*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c*c*c)-NFFT_K(6720.0)*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c); break;
    case  9 : value = -NFFT_K(512.0)*NFFT_POW(x,NFFT_K(11.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(18.0))+NFFT_K(13824.0)*x*x*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(16.0))-NFFT_K(122112.0)*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(14.0))+NFFT_K(419328.0)*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(12.0))-NFFT_K(514080.0)*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0))+NFFT_K(151200.0)*x*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c*c*c); break;
    case 10 : value = NFFT_K(1024.0)*NFFT_POW(x,NFFT_K(12.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(20.0))-NFFT_K(33280.0)*NFFT_POW(x,NFFT_K(10.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(18.0))+NFFT_K(368640.0)*x*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(16.0))-NFFT_K(1693440.0)*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(14.0))+NFFT_K(3124800.0)*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(12.0))-NFFT_K(1844640.0)*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0))+NFFT_K(151200.0)*NFFT_EXP(-x*x/(c*c))/(c*c*c*c*c*c*c*c); break;
    case 11 : value = -NFFT_K(2048.0)*NFFT_POW(x,NFFT_K(13.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(22.0))+NFFT_K(78848.0)*NFFT_POW(x,NFFT_K(11.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(20.0))-NFFT_K(1070080.0)*x*x*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(18.0))+NFFT_K(6336000.0)*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(16.0))-NFFT_K(16410240.0)*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(14.0))+NFFT_K(16188480.0)*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(12.0))-NFFT_K(3991680.0)*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0)); break;
    case 12 : value = NFFT_K(4096.0)*NFFT_POW(x,NFFT_K(14.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(24.0))-NFFT_K(184320.0)*NFFT_POW(x,NFFT_K(12.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(22.0))+NFFT_K(3007488.0)*NFFT_POW(x,NFFT_K(10.0))*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(20.0))-NFFT_K(22302720.0)*x*x*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(18.0))+NFFT_K(77172480.0)*x*x*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(16.0))-NFFT_K(114428160.0)*x*x*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(14.0))+NFFT_K(56548800.0)*x*x*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(12.0))-NFFT_K(3991680.0)*NFFT_EXP(-x*x/(c*c))/NFFT_POW(c,NFFT_K(10.0)); break;
    default : value = NFFT_K(0.0);
  }

  return value / (c*c);
}

NFFT_C absx(NFFT_R x, int der, const NFFT_R *param)    /* K(x)=|x| */
{
  NFFT_R value=NFFT_K(0.0);

  (void)param;
  
  if (der == 0) value=NFFT_FABS(x);
  else if (der == 1){
    if (x<0) value=NFFT_K(-1.0);
    else value=NFFT_K(1.0);
  }
  else value=NFFT_K(0.0);
  
  return value;
}

/* \} */

/* kernels.c */
