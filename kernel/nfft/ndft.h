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

/* The two planner-native direct NDFT solvers (ndft-1d.c, ndft-nd.c) and the
 * plan type they share. Include after nfft3.h, infft.h and iplanner.h. Both
 * solvers read x/f_hat/f straight from problem_nfft.
 *
 * Blocked phase recurrence: one reduced phase serves NDFT_RECURRENCE_BLOCK
 * terms, which advance by a complex multiply with the fixed step
 * exp(-+i 2pi x). Two transcendentals per block replace two per term, and the
 * small argument keeps the phase error independent of N. */

#ifndef NFFT_NDFT_H
#define NFFT_NDFT_H

#define NDFT_RECURRENCE_BLOCK 32

/* Phase 2pi(base + k x) reduced to ~[-1/2,1/2) with a single-rounding FMA.
 * base is what the outer axes contribute (0 in one dimension) and must itself
 * already be reduced. */
static inline R Y(nfft_ndft_reduced_omega)(const R k, const R x, const R base) {
  const R s = FFMA(k, x, base);
  return K2PI * FFMA(k, x, base - RINT(s));
}

plan *Y(nfft_ndft_make_plan)(double pcost,
                             void (*fwd)(const problem_nfft *),
                             void (*adj)(const problem_nfft *),
                             const char *reg_nam);

/* Direct NDFT cost estimate, both ranks. */
double Y(nfft_ndft_pcost)(const problem *p);

#endif /* NFFT_NDFT_H */
