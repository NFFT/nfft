/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: util.c 3483 2010-04-23 19:02:34Z keiner $ */

#include "infft.h"

static const R P1[15] = {
    K(-2.2335582639474375249e+15),
    K(-5.5050369673018427753e+14),
    K(-3.2940087627407749166e+13),
    K(-8.4925101247114157499e+11),
    K(-1.1912746104985237192e+10),
    K(-1.0313066708737980747e+08),
    K(-5.9545626019847898221e+05),
    K(-2.4125195876041896775e+03),
    K(-7.0935347449210549190e+00),
    K(-1.5453977791786851041e-02),
    K(-2.5172644670688975051e-05),
    K(-3.0517226450451067446e-08),
    K(-2.6843448573468483278e-11),
    K(-1.5982226675653184646e-14),
    K(-5.2487866627945699800e-18)
};

static const R Q1[6] = {
    K(-2.2335582639474375245e+15),
    K(7.8858692566751002988e+12),
    K(-1.2207067397808979846e+10),
    K(1.0377081058062166144e+07),
    K(-4.8527560179962773045e+03),
    K(1.0)
};

static const R P2[7] = {
    K(-2.2210262233306573296e-04),
    K(1.3067392038106924055e-02),
    K(-4.4700805721174453923e-01),
    K(5.5674518371240761397e+00),
    K(-2.3517945679239481621e+01),
    K(3.1611322818701131207e+01),
    K(-9.6090021968656180000e+00)
};

static const R Q2[8] = {
    K(-5.5194330231005480228e-04),
    K(3.2547697594819615062e-02),
    K(-1.1151759188741312645e+00),
    K(1.3982595353892851542e+01),
    K(-6.0228002066743340583e+01),
    K(8.5539563258012929600e+01),
    K(-3.1446690275135491500e+01),
    K(1.0)
};

static inline R evaluate_polynomial(const R* poly, const R z, const int count)
{
   //A(count > 0)
   R sum = poly[count - 1];
   int i;

   for (i = count - 2; i >= 0; --i)
   {
     sum *= z;
     sum += poly[i];
   }

   return sum;
}

R X(bessel_i0)(R x)
{
  if (x < 0)
  {
    /* even function */
    x = -x;
  }

  if (x == K(0.0))
    return K(1.0);

  if (x <= K(15.0))
  {
    /* x in (0, 15] */
    const R y = x * x;
    return (evaluate_polynomial(P1, y, 15) / evaluate_polynomial(Q1, y, 6));
  }
  else
  {
    /* x in (15, \infty) */
    const R y = K(1.0) / x - K(1.0) / K(15.0);
    return ((EXP(x) / SQRT(x)) * (evaluate_polynomial(P2, y, 7) /
      evaluate_polynomial(Q2, y, 8)));
  }
}
