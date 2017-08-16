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

/** Computes the inner/dot product \f$x^H x\f$. */
R Y(dot_complex)(C *x, INT n)
{
  INT k;
  R dot;

  for (k = 0, dot = K(0.0); k < n; k++)
    dot += CONJ(x[k])*x[k];

  return dot;
}

/** Computes the inner/dot product \f$x^H x\f$. */
R Y(dot_double)(R *x, INT n)
{
  INT k;
  R dot;

  for (k = 0, dot = K(0.0); k < n; k++)
    dot += x[k]*x[k];

  return dot;
}


/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$. */
R Y(dot_w_complex)(C *x, R *w, INT n)
{
  INT k;
  R dot;

  for (k = 0, dot = K(0.0); k < n; k++)
    dot += w[k]*CONJ(x[k])*x[k];

  return dot;
}

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$. */
R Y(dot_w_double)(R *x, R *w, INT n)
{
  INT k;
  R dot;

  for (k = 0, dot = K(0.0); k < n; k++)
    dot += w[k]*x[k]*x[k];

  return dot;
}


/** Computes the weighted inner/dot product \f$x^H (w\odot w2\odot w2 \odot x)\f$. */
R Y(dot_w_w2_complex)(C *x, R *w, R *w2, INT n)
{
  INT k;
  R dot;

  for (k = 0, dot = K(0.0); k < n; k++)
    dot += w[k]*w2[k]*w2[k]*CONJ(x[k])*x[k];

  return dot;
}

/** Computes the weighted inner/dot product \f$x^H (w2\odot w2 \odot x)\f$. */
R Y(dot_w2_complex)(C *x, R *w2, INT n)
{
  INT k;
  R dot;

  for (k = 0, dot = K(0.0); k < n; k++)
    dot+=w2[k]*w2[k]*CONJ(x[k])*x[k];

  return dot;
}
