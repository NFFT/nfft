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

#include "api.h"
#include "cstripack.h"

/** Compute non periodic voronoi weights for ordered nodes x_j */
void Y(voronoi_weights_1d)(R *w, R *x, const INT M)
{
  INT j;

  w[0] = (x[1]-x[0])/K(2.0);

  for(j = 1; j < M-1; j++)
    w[j] = (x[j+1]-x[j-1])/K(2.0);

  w[M-1] = (x[M-1]-x[M-2])/K(2.0);
}
