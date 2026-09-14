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

/* Shared infrastructure for Step C (node convolution, matrix B) of the NFFT
 * decomposition. Include after nfft3.h, infft.h and iplanner.h. */

#ifndef NFFT_CONV_H
#define NFFT_CONV_H

/* FMA fusion under -ffast-math can vary by call site, rounding the same
 * accumulation order differently for band-clipped calls vs whole-axis ones.
 * GCC ignores STDC FP_CONTRACT under -ffast-math, hence this override. As a
 * #pragma it covers every translation unit that includes this header, not
 * just the threaded solver -- every serial CONV solver gets it too. */
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC optimize("fp-contract=off")
#endif

/* Split one axis' window support -- len taps starting at wrapped grid cell u --
 * into two contiguous runs (the second is empty when the support does not wrap
 * past n). tof[r] indexes psi, gof[r] indexes the grid, rl[r] is the length. */
static inline void Y(conv_runs)(INT u, INT n, INT len, INT *tof, INT *gof,
                             INT *rl)
{
  INT head = n - u;
  if (head > len)
    head = len;
  tof[0] = 0;
  gof[0] = u;
  rl[0] = head;
  tof[1] = head;
  gof[1] = 0;
  rl[1] = len - head;
}

/* Intersect the grid run [*gof, *gof + *rl) with [lo, hi), carrying the psi
 * offset with it. Leaves *rl <= 0 when the run falls outside, which every
 * caller's loop then skips. A serial caller passes the whole axis and this is
 * a no-op. */
static inline void Y(conv_clip)(INT lo, INT hi, INT *tof, INT *gof, INT *rl)
{
  INT end = *gof + *rl;
  if (*gof < lo) {
    *tof += lo - *gof;
    *gof = lo;
  }
  if (end > hi)
    end = hi;
  *rl = end - *gof;
}

/* One node's tap block, rank 1. Forward returns the gathered contribution and
 * leaves g alone; adjoint scatter-adds fj and returns 0. [lo, hi) bounds the
 * grid rows on axis 0 the call may touch: the serial solvers pass (0, n), the
 * threaded adjoint passes the band it owns. */
static inline C
Y(conv_node_1d)(C *g, const R *psi, INT u, INT n, INT len, INT lo,
                             INT hi, C fj, int forward)
{
  INT tof[2], gof[2], rl[2], k;
  C acc = K(0.0);
  int r;
  Y(conv_runs)(u, n, len, tof, gof, rl);
  for (r = 0; r < 2; r++) {
    Y(conv_clip)(lo, hi, &tof[r], &gof[r], &rl[r]);
    if (forward)
      for (k = 0; k < rl[r]; k++)
        acc += g[gof[r] + k] * psi[tof[r] + k];
    else
      for (k = 0; k < rl[r]; k++)
        g[gof[r] + k] += fj * psi[tof[r] + k];
  }
  return acc;
}

/* Rank 2. Axis 1 is contiguous, so its runs do not vary over axis 0 and are
 * hoisted; axis 0 carries the clip. */
static inline C Y(conv_node_2d)(C *g, const R *psi0, const R *psi1, INT u0,
                             INT u1, INT n0, INT n1, INT len, INT lo, INT hi,
                             C fj, int forward)
{
  INT tof0[2], gof0[2], rl0[2], tof1[2], gof1[2], rl1[2], i, k;
  C acc = K(0.0);
  int a;
  Y(conv_runs)(u0, n0, len, tof0, gof0, rl0);
  Y(conv_runs)(u1, n1, len, tof1, gof1, rl1);
  {
    const INT ka = rl1[0], kb = rl1[1], goa = gof1[0], gob = gof1[1];
    const R *p1a = psi1 + tof1[0], *p1b = psi1 + tof1[1];
    for (a = 0; a < 2; a++) {
      Y(conv_clip)(lo, hi, &tof0[a], &gof0[a], &rl0[a]);
      for (i = 0; i < rl0[a]; i++) {
        const R p0 = psi0[tof0[a] + i];
        C *grow = g + (gof0[a] + i) * n1;
        C *ga = grow + goa, *gb = grow + gob;
        if (forward) {
          C sub = K(0.0);
          for (k = 0; k < ka; k++)
            sub += ga[k] * p1a[k];
          for (k = 0; k < kb; k++)
            sub += gb[k] * p1b[k];
          acc += sub * p0;
        } else {
          const C fp = fj * p0;
          for (k = 0; k < ka; k++)
            ga[k] += fp * p1a[k];
          for (k = 0; k < kb; k++)
            gb[k] += fp * p1b[k];
        }
      }
    }
  }
  return acc;
}

/* Rank 3. Axis 2 is contiguous and hoisted; axis 0 carries the clip. */
static inline C
Y(conv_node_3d)(C *g, const R *psi0, const R *psi1, const R *psi2,
                             INT u0, INT u1, INT u2, INT n0, INT n1, INT n2,
                             INT len, INT lo, INT hi, C fj, int forward)
{
  INT tof0[2], gof0[2], rl0[2], tof1[2], gof1[2], rl1[2];
  INT tof2[2], gof2[2], rl2[2], i, k, l;
  C acc = K(0.0);
  int a, b;
  Y(conv_runs)(u0, n0, len, tof0, gof0, rl0);
  Y(conv_runs)(u1, n1, len, tof1, gof1, rl1);
  Y(conv_runs)(u2, n2, len, tof2, gof2, rl2);
  {
    const INT la = rl2[0], lb = rl2[1], goa = gof2[0], gob = gof2[1];
    const R *p2a = psi2 + tof2[0], *p2b = psi2 + tof2[1];
    for (a = 0; a < 2; a++) {
      Y(conv_clip)(lo, hi, &tof0[a], &gof0[a], &rl0[a]);
      for (i = 0; i < rl0[a]; i++) {
        const R p0 = psi0[tof0[a] + i];
        C *gplane = g + (gof0[a] + i) * n1 * n2;
        for (b = 0; b < 2; b++) {
          const INT nk = rl1[b], to1 = tof1[b], go1 = gof1[b];
          for (k = 0; k < nk; k++) {
            const R p01 = p0 * psi1[to1 + k];
            C *grow = gplane + (go1 + k) * n2;
            C *ga = grow + goa, *gb = grow + gob;
            if (forward) {
              C sub = K(0.0);
              for (l = 0; l < la; l++)
                sub += ga[l] * p2a[l];
              for (l = 0; l < lb; l++)
                sub += gb[l] * p2b[l];
              acc += sub * p01;
            } else {
              const C fp = fj * p01;
              for (l = 0; l < la; l++)
                ga[l] += fp * p2a[l];
              for (l = 0; l < lb; l++)
                gb[l] += fp * p2b[l];
            }
          }
        }
      }
    }
  }
  return acc;
}

/* Rank >= 4. Out of line: at (2m+2)^d taps per node the call costs nothing,
 * and the body is too large to inline into two translation units. */
C Y(conv_node_nd)(C *g, const R *psi, const INT *u, const INT *n, int d, INT len,
               INT lo, INT hi, C fj, int forward);

/* Awake-time tables, shared by every CONV solver so serial and threaded plans
 * build them through the same code. perm == NULL keeps caller order;
 * otherwise row k belongs to node perm[k]. */
void Y(conv_fill_u)(const R *x, const INT *perm, const INT *n, int d, int m,
                 INT M, INT *u);
void Y(conv_fill_psi)(int window, const R *x, const INT *perm, const INT *n,
                   const INT *N, int d, int m, INT M, R *psi);

/* CONV solvers. */
void Y(conv_solver_1d_register)(planner *pl);
void Y(conv_solver_2d_register)(planner *pl);
void Y(conv_solver_3d_register)(planner *pl);
void Y(conv_solver_nd_register)(planner *pl);

#endif /* NFFT_CONV_H */
