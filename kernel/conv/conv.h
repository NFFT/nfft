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

/* CONV solvers. */
void Y(conv_solver_1d_register)(planner *pl);
void Y(conv_solver_2d_register)(planner *pl);
void Y(conv_solver_3d_register)(planner *pl);
void Y(conv_solver_nd_register)(planner *pl);

#endif /* NFFT_CONV_H */
