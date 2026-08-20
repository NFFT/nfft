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

/* Shared infrastructure for the planner-native direct NDFT solvers
 * (kernel/nfft/ndft.c and ndft-nd.c).  Include after nfft3.h, infft.h and
 * iplanner.h.  These solvers are self-contained: apply reads x/f_hat/f straight
 * from problem_nfft and owns no shared legacy core. */

#ifndef NFFT_NDFT_H
#define NFFT_NDFT_H

/* Allocate a native NDFT plan. */
plan *Y(nfft_ndft_make_plan)(double pcost,
                             void (*apply_native_fwd)(const problem_nfft *),
                             void (*apply_native_adj)(const problem_nfft *),
                             const char *reg_nam);

/* Plain direct NDFT cost estimate. */
double Y(nfft_ndft_plain_pcost)(const problem *p);

#endif /* NFFT_NDFT_H */
