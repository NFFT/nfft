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

#ifndef WINDOW_TEST_H
#define WINDOW_TEST_H

#include "infft.h"

/* The runtime window vtable must reproduce the compile-time PHI_HUT/PHI macros
 * (same bessel_i0, same b = pi(2-1/sigma)). */
void Y(check_window_vtable)(void);
/* The compile-time window reported as an NFFT_WINDOW_* ordinal. */
void Y(check_window_id)(void);
/* The KB vtable is peak-normalized (phi_hut(0)=1); ratios between window values
 * are unchanged from the legacy (unnormalized) shape. */
void Y(check_window_normalized)(void);
/* Gaussian/B-spline/sinc (raw, unnormalized) alongside KB; Dirac and
 * out-of-range ordinals return 0 (Dirac is removed from the new API). */
void Y(check_window_all)(void);
/* The window is part of the wisdom key: two problems differing only in window
 * hash differently. */
void Y(check_window_key)(void);
/* Robust log(I0) and scaled I0: log I0 matches LOG(bessel_i0) in the small
 * branch and the asymptotic expansion for large x, and bessel_i0_scaled stays
 * in (0, 1] and finite without ever materializing I0(x). */
void Y(check_window_bessel_log_scaled)(void);
/* Near the peak the rationalized KB evaluator beats the naive exp(a - lg_peak)
 * fold's ~3e-6 float residual; also the identity log I0(x) == x + logtail(x). */
void Y(check_window_cancellation)(void);
/* The range-apply API hoists the per-axis window constant across a whole axis.
 * phi_hut_apply fills the band identically to the single-arg window_phi_hut;
 * phi_precompute builds the per-node psi identically to window_phi. */
void Y(check_window_apply)(void);

#endif
