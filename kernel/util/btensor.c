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

/* The tensor product kernels of B and B^T (infft.h) for any dimension d >= 2,
 * taken by Y(B_tensor_trafo) and Y(B_tensor_adjoint) for d >= 4: an odometer
 * over the window indices of dimensions 0, ..., d-2, with the prefix products
 * of the weights and the prefix row offsets updated from the changed dimension
 * on, and one row of the last dimension per step. */

#include "config.h"

#include "nfft3.h"
#include "infft.h"

void Y(B_tensor_trafo_nd)(const INT d, const INT L, const INT es,
    const INT *np, const INT *lo, const R *const *w, const R *wl,
    const R *g, R *acc)
{
  const INT len = es * L;
  INT lj[d - 1], off[d - 1], s, t = 0;
  R prod[d - 1];

  for (s = 0; s < d - 1; s++)
    lj[s] = 0;

  for (;;)
  {
    for (s = t; s < d - 1; s++)
    {
      prod[s] = (s == 0 ? K(1.0) : prod[s - 1]) * w[s][lj[s]];
      off[s] = (s == 0 ? 0 : off[s - 1] * np[s]) + lo[s] + lj[s];
    }

    Y(B_row_trafo)(acc, wl, g + es * (off[d - 2] * np[d - 1] + lo[d - 1]),
        prod[d - 2], len);

    for (t = d - 2; t >= 0 && ++lj[t] == L; t--)
      lj[t] = 0;
    if (t < 0)
      break;
  }
}

void Y(B_tensor_adjoint_nd)(const INT d, const INT L, const INT es,
    const INT *np, const INT *lo, const R *const *w, const R *wf, R *g)
{
  const INT len = es * L;
  INT lj[d - 1], off[d - 1], s, t = 0;
  R prod[d - 1];

  for (s = 0; s < d - 1; s++)
    lj[s] = 0;

  for (;;)
  {
    for (s = t; s < d - 1; s++)
    {
      prod[s] = (s == 0 ? K(1.0) : prod[s - 1]) * w[s][lj[s]];
      off[s] = (s == 0 ? 0 : off[s - 1] * np[s]) + lo[s] + lj[s];
    }

    Y(B_row_adjoint)(g + es * (off[d - 2] * np[d - 1] + lo[d - 1]), wf,
        prod[d - 2], len);

    for (t = d - 2; t >= 0 && ++lj[t] == L; t--)
      lj[t] = 0;
    if (t < 0)
      break;
  }
}
