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

/** Copies \f$x \leftarrow y\f$. */
void Y(cp_complex)(C *x, C *y, INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = y[k];
}

/** Copies \f$x \leftarrow y\f$. */
void Y(cp_double)(R *x, R *y, INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = y[k];
}

/** Copies \f$x \leftarrow a y\f$. */
void Y(cp_a_complex)(C *x, R a, C *y, INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = a * y[k];
}

/** Copies \f$x \leftarrow a y\f$. */
void Y(cp_a_double)(R *x, R a, R *y, INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = a * y[k];
}


/** Copies \f$x \leftarrow w\odot y\f$. */
void Y(cp_w_complex)(C *x, R *w, C *y, INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = w[k]*y[k];
}

/** Copies \f$x \leftarrow w\odot y\f$. */
void Y(cp_w_double)(R *x, R *w, R *y, INT n)
{
  INT k;

  for (k = 0; k < n; k++)
    x[k] = w[k] * y[k];
}
