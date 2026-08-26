/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* 3D CONV solver: Step C of the fast NFFT decomposition -- the node convolution 
 * (matrix B). Sums the oversampled grid g against the window psi at each 
 * nonequispaced node x_j (forward), or scatter-adds f onto g with the same psi 
 * weights (adjoint). psi depends on x/window/n/N/m), precomputed once at awake
 * (sparse PRE_PSI strategy in legacy code). */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

typedef struct
{
  plan super;
  INT n0, n1, n2, N0, N1, N2, M; /* geometry captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi;     /* length M*3*(2m+2): psi[(j*3+t)*(2m+2)+lj], built at awake */
  int precomputed;
} conv_3d_plan;

/* uo2: neighbor window start/end on axis of size n, wrapped mod n. */
static void uo2(INT *u, INT *o, const R x, const INT n, const INT m) {
  INT c = LRINT(FLOOR(x * (R)n));
  *u = (c - m + n) % n;
  *o = (c + 1 + m + n) % n;
}

/* precompute the psi table from ego->x. */
static void conv_3d_awake(plan *ego_, int wakefulness) {
  conv_3d_plan *pln = (conv_3d_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO) {
    if (!pln->precomputed) {
      int t;
      INT nn[3];
      nn[0] = pln->n0;
      nn[1] = pln->n1;
      nn[2] = pln->n2;
      INT NN[3];
      NN[0] = pln->N0;
      NN[1] = pln->N1;
      NN[2] = pln->N2;
      for (t = 0; t < 3; t++)
        Y(window_phi_precompute)
      (pln->window, nn[t], NN[t], pln->m,
       pln->x + t, 3, pln->M,
       pln->psi + t * (2 * pln->m + 2), 3 * (2 * pln->m + 2));
      pln->precomputed = 1;
    }
  } else
    pln->precomputed = 0; /* -> SLEEPY: psi values now stale */
}

/* Forward B (g -> f[j]). */
static void conv_trafo_3d_compute(C *fj, const C *g, const R *psij_const0,
                                  const R *psij_const1, const R *psij_const2, const R *xj0, const R *xj1,
                                  const R *xj2, const INT n0, const INT n1, const INT n2, const int m) {
  INT u0, o0, l0, u1, o1, l1, u2, o2, l2;
  const C *gj;
  const R *psij0, *psij1, *psij2;

  psij0 = psij_const0;
  psij1 = psij_const1;
  psij2 = psij_const2;

  uo2(&u0, &o0, *xj0, n0, m);
  uo2(&u1, &o1, *xj1, n1, m);
  uo2(&u2, &o2, *xj2, n2, m);

  *fj = K(0.0);

  if (u0 < o0)
    if (u1 < o1)
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      else
        /* asserts (u2>o2)*/
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
    else /* asserts (u1>o1)*/
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      else /* asserts (u2>o2) */
      {
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      }
  else /* asserts (u0>o0) */
    if (u1 < o1)
      if (u2 < o2) {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      } else /* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      }
    else /* asserts (u1>o1) */
      if (u2 < o2) {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      } else /* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + (l0 * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      }
}

/* Adjoint B^H (f[j] -> g, scatter-add). */
static void conv_adjoint_3d_compute_serial(const C *fj, C *g,
                                           const R *psij_const0, const R *psij_const1, const R *psij_const2,
                                           const R *xj0, const R *xj1, const R *xj2, const INT n0, const INT n1,
                                           const INT n2, const int m) {
  INT u0, o0, l0, u1, o1, l1, u2, o2, l2;
  C *gj;
  const R *psij0, *psij1, *psij2;

  psij0 = psij_const0;
  psij1 = psij_const1;
  psij2 = psij_const2;

  uo2(&u0, &o0, *xj0, n0, m);
  uo2(&u1, &o1, *xj1, n1, m);
  uo2(&u2, &o2, *xj2, n2, m);

  if (u0 < o0)
    if (u1 < o1)
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      else
        /* asserts (u2>o2)*/
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
    else /* asserts (u1>o1)*/
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      else /* asserts (u2>o2) */
      {
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      }
  else /* asserts (u0>o0) */
    if (u1 < o1)
      if (u2 < o2) {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      } else /* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      }
    else /* asserts (u1>o1) */
      if (u2 < o2) {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      } else /* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++) {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++) {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + (l0 * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      }
}

static void conv_3d_apply(const plan *ego_, const problem *p) {
  const conv_3d_plan *pln = (const conv_3d_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n0 = pln->n0, n1 = pln->n1, n2 = pln->n2, M = pln->M;
  int m = pln->m;
  const C *g = pc->g;
  C *f = pc->f;
  INT j;
  for (j = 0; j < M; j++) {
    const R *psij0 = &pln->psi[(j * 3 + 0) * (2 * m + 2)];
    const R *psij1 = &pln->psi[(j * 3 + 1) * (2 * m + 2)];
    const R *psij2 = &pln->psi[(j * 3 + 2) * (2 * m + 2)];
    const R *xj0 = &pln->x[j * 3 + 0];
    const R *xj1 = &pln->x[j * 3 + 1];
    const R *xj2 = &pln->x[j * 3 + 2];
    conv_trafo_3d_compute(&f[j], g, psij0, psij1, psij2, xj0, xj1, xj2, n0,
                          n1, n2, m);
  }
}

static void conv_3d_apply_adjoint(const plan *ego_, const problem *p) {
  const conv_3d_plan *pln = (const conv_3d_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n0 = pln->n0, n1 = pln->n1, n2 = pln->n2, M = pln->M;
  INT ntot = n0 * n1 * n2;
  int m = pln->m;
  const C *f = pc->f;
  C *g = pc->g;
  INT j;
  /* The scatter accumulates (+=) into an overlapping, node-dependent set that
   * does not cover the grid, so the whole grid must start zeroed. */
  memset(g, 0, (size_t)ntot * sizeof(C)); /* zero the oversampled grid */
  for (j = 0; j < M; j++) {
    const R *psij0 = &pln->psi[(j * 3 + 0) * (2 * m + 2)];
    const R *psij1 = &pln->psi[(j * 3 + 1) * (2 * m + 2)];
    const R *psij2 = &pln->psi[(j * 3 + 2) * (2 * m + 2)];
    const R *xj0 = &pln->x[j * 3 + 0];
    const R *xj1 = &pln->x[j * 3 + 1];
    const R *xj2 = &pln->x[j * 3 + 2];
    conv_adjoint_3d_compute_serial(&f[j], g, psij0, psij1, psij2, xj0, xj1,
                                   xj2, n0, n1, n2, m);
  }
}

static void conv_3d_print(const plan *ego_, printer *pr) {
  const conv_3d_plan *pln = (const conv_3d_plan *)ego_;
  pr->print(pr, "(conv_solver_3d pcost=%D)", (INT)pln->super.pcost);
}
static void conv_3d_destroy(plan *ego_) {
  conv_3d_plan *pln = (conv_3d_plan *)ego_;
  Y(free)
  (pln->psi);
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_3d_plan_adt = {conv_3d_apply, conv_3d_awake,
                                          conv_3d_print, conv_3d_destroy,
                                          conv_3d_apply_adjoint};

/* d == 3 only */
static plan *mkplan_conv_3d(const solver *ego, const problem *p, planner *pl) {
  const problem_conv *pc = (const problem_conv *)p;
  conv_3d_plan *pln;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pc->sz->rnk != 3)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL ||
      pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */

  pln = (conv_3d_plan *)Y(plan_create)(sizeof(conv_3d_plan), &conv_3d_plan_adt);
  pln->n0 = Y(problem_conv_n)(p, 0);
  pln->n1 = Y(problem_conv_n)(p, 1);
  pln->n2 = Y(problem_conv_n)(p, 2);
  pln->N0 = Y(problem_conv_N)(p, 0);
  pln->N1 = Y(problem_conv_N)(p, 1);
  pln->N2 = Y(problem_conv_N)(p, 2);
  pln->M = pc->M;
  pln->m = pc->m;
  pln->window = pc->window;
  pln->x = pc->x; /* borrowed alias */
  pln->psi = (R *)Y(malloc)((size_t)pln->M * 3 * (size_t)(2 * pln->m + 2) * sizeof(R));
  pln->precomputed = 0;
  pln->super.pcost = Y(conv_b_pcost)(p);
  return &pln->super;
}

static const solver_adt conv_3d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_3d};
void Y(conv_solver_3d_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_3d_adt));
}
