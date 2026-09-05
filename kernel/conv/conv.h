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

double Y(conv_b_pcost)(const problem *p);

/* Split one axis' window support -- len taps starting at wrapped grid cell u --
 * into two contiguous runs (the second is empty when the support does not wrap
 * past n). tof[r] indexes psi, gof[r] indexes the grid, rl[r] is the length. */
static inline void Y(conv_runs)(INT u, INT n, INT len, INT *tof, INT *gof,
                                INT *rl) {
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

/* Placeholder window starts for PLNR_AWAKE_ZERO, u[j*d+t]. All-zero starts
 * would make every node read the same grid cells, so a raced candidate would
 * time an access pattern far more cache-friendly than the real scattered one.
 * Walking each axis by a stride keeps the timing honest and still evaluates no
 * window. */
static inline void Y(conv_spread_u)(INT *u, INT M, int d, const INT *n) {
  int t;
  for (t = 0; t < d; t++) {
    INT nt = n[t], step = (nt > 1) ? (nt / 2 + 1) : 0;
    INT c = 0, j;
    for (j = 0; j < M; j++) {
      u[j * d + t] = c;
      c += step;
      if (c >= nt)
        c -= nt;
    }
  }
}

/* CONV solvers. */
void Y(conv_solver_1d_register)(planner *pl);
void Y(conv_solver_2d_register)(planner *pl);
void Y(conv_solver_3d_register)(planner *pl);
void Y(conv_solver_nd_register)(planner *pl);

#endif /* NFFT_CONV_H */
