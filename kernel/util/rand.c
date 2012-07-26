/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

R X(drand48)(void)
{
#ifdef HAVE_DRAND48
  return drand48();
#else
  return ((R)rand())/((R)RAND_MAX);
#endif
}

void X(srand48)(long int seed)
{
#ifdef HAVE_SRAND48
  srand48(seed);
#else
  srand((unsigned int)seed);
#endif
}

void X(vrand_unit_complex)(C *x, const INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = X(drand48)() + II*X(drand48)();
}

void X(vrand_shifted_unit_double)(R *x, const INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = X(drand48)() - K(0.5);
}
