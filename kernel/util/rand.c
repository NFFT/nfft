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

#include "infft.h"

R Y(drand48)(void)
{
#ifdef HAVE_DRAND48
  return (R)(drand48());
#else
  return ((R)rand())/((R)RAND_MAX);
#endif
}

void Y(srand48)(long int seed)
{
#ifdef HAVE_SRAND48
  srand48(seed);
#else
  srand((unsigned int)seed);
#endif
}

void Y(vrand_unit_complex)(C *x, const INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = Y(drand48)() + II * Y(drand48)();
}

void Y(vrand_shifted_unit_double)(R *x, const INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = Y(drand48)() - K(0.5);
}

void Y(vrand_real)(R *x, const INT n, const R a, const R b)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = a + Y(drand48)() * (b - a);
}
