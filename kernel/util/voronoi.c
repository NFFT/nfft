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

/** Compute non periodic voronoi weights for ordered nodes x_j */
void X(voronoi_weights_1d)(R *w, R *x, const INT M)
{
  INT j;

  w[0] = (x[1]-x[0])/K(2.0);

  for(j = 1; j < M-1; j++)
    w[j] = (x[j+1]-x[j-1])/K(2.0);

  w[M-1] = (x[M-1]-x[M-2])/K(2.0);
}

void X(voronoi_weights_S2)(R *w, R *xi, INT M)
{
  R *x;
  R *y;
  R *z;
  int j;
  int k;
  int el;
  int Mlocal = M;
  int lnew;
  int ier;
  int *list;
  int *lptr;
  int *lend;
  int *near;
  int *next;
  R  *dist;
  int *ltri;
  int *listc;
  int nb;
  R *xc;
  R *yc;
  R *zc;
  R *rc;
  R *vr;
  int lp;
  int lpl;
  int kv;
  R a;

  /* Allocate memory for auxilliary arrays. */
  x = (R*)X(malloc)(M * sizeof(R));
  y = (R*)X(malloc)(M * sizeof(R));
  z = (R*)X(malloc)(M * sizeof(R));

  list = (int*)X(malloc)((6*M-12+1)*sizeof(int));
  lptr = (int*)X(malloc)((6*M-12+1)*sizeof(int));
  lend = (int*)X(malloc)((M+1)*sizeof(int));
  near = (int*)X(malloc)((M+1)*sizeof(int));
  next = (int*)X(malloc)((M+1)*sizeof(int));
  dist = (R*)X(malloc)((M+1)*sizeof(R));
  ltri = (int*)X(malloc)((6*M+1)*sizeof(int));
  listc = (int*)X(malloc)((6*M-12+1)*sizeof(int));
  xc = (R*)X(malloc)((2*M-4+1)*sizeof(R));
  yc = (R*)X(malloc)((2*M-4+1)*sizeof(R));
  zc = (R*)X(malloc)((2*M-4+1)*sizeof(R));
  rc = (R*)X(malloc)((2*M-4+1)*sizeof(R));
  vr = (R*)X(malloc)(3*(2*M-4+1)*sizeof(R));

  /* Convert from spherical Coordinates in [0,1/2]x[-1/2,1/2) to Cartesian
   * coordinates. */
  for (k = 0; k < M; k++)
  {
    x[k] = SIN(K2PI*xi[2*k+1])*COS(K2PI*xi[2*k]);
    y[k] = SIN(K2PI*xi[2*k+1])*SIN(K2PI*xi[2*k]);
    z[k] = COS(K2PI*xi[2*k+1]);
  }

  /* Generate Delaunay triangulation. */
  trmesh_(&Mlocal, x, y, z, list, lptr, lend, &lnew, near, next, dist, &ier);

  /* Check error flag. */
  if (ier == 0)
  {
    /* Get Voronoi vertices. */
    crlist_(&Mlocal, &Mlocal, x, y, z, list, lend, lptr, &lnew, ltri, listc, &nb, xc,
      yc, zc, rc, &ier);

    if (ier == 0)
    {
      /* Calcuate sizes of Voronoi regions. */
      for (k = 0; k < M; k++)
      {
        /* Get last neighbour index. */
        lpl = lend[k];
        lp = lpl;

        j = 0;
        vr[3*j] = x[k];
        vr[3*j+1] = y[k];
        vr[3*j+2] = z[k];

        do
        {
          j++;
          /* Get next neighbour. */
          lp = lptr[lp-1];
          kv = listc[lp-1];
          vr[3*j] = xc[kv-1];
          vr[3*j+1] = yc[kv-1];
          vr[3*j+2] = zc[kv-1];
          /* fprintf(stderr, "lp = %ld\t", lp); */
        } while (lp != lpl);

        a = 0;
        for (el = 0; el < j; el++)
        {
          a += areas_(vr, &vr[3*(el+1)],&vr[3*(((el+1)%j)+1)]);
        }

        w[k] = a;
      }
    }
  }

  /* Deallocate memory. */
  X(free)(x);
  X(free)(y);
  X(free)(z);

  X(free)(list);
  X(free)(lptr);
  X(free)(lend);
  X(free)(near);
  X(free)(next);
  X(free)(dist);
  X(free)(ltri);
  X(free)(listc);
  X(free)(xc);
  X(free)(yc);
  X(free)(zc);
  X(free)(rc);
  X(free)(vr);
}
