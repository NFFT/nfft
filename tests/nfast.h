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

#ifndef NFAST_TEST_H
#define NFAST_TEST_H

#include "infft.h"

/* The fast-NFFT decomposition problem types (DECONV/CONV) carry
 * exactly the attributes their stage needs and omit the ones that must be
 * absent; hash is stable/discriminating and (for CONV) data-blind in x. */
void Y(check_nfast_window_id)(void);
void Y(check_nfast_deconv_problem)(void);
void Y(check_nfast_conv_problem)(void);
/* The runtime KB window vtable must reproduce the compile-time
 * PHI_HUT/PHI macros (same bessel_i0, same b = pi(2-1/sigma)). */
void Y(check_nfast_window_vtable)(void);
/* The KB vtable is peak-normalized (phi_hut(0)=1); ratios between
 * window values are unchanged from the legacy (unnormalized) shape. */
void Y(check_nfast_window_normalized)(void);
/* The vtable ports Gaussian/B-spline/sinc (raw, unnormalized) from
 * the infft.h macros alongside KB; Dirac and out-of-range ordinals return 0
 * (Dirac is removed from the new API). */
void Y(check_nfast_window_all)(void);
/* Window is in the wisdom key; wrapper fasts decline a non-
 * compile-time window. */
void Y(check_nfast_window_key)(void);
/* The DECONV solver -- plans a 1D KB deconvolution
 * problem via the process-global planner and asserts output values directly.
 * Forward: single spike -> zero-padded grid divided by phi_hut.  Adjoint (D^H):
 * multiplies by phi_hut_inv (1/phi_hut -- a real diagonal is self-adjoint) and
 * gathers; asserted at k=0 and k=1 from a clean grid (NOT a round-trip, which
 * the inverse would also pass).  Also covers the type-II +1 frequency shift. */
void Y(check_nfast_deconv_solver)(void);
/* rank-1 DECONV across odd N, type-II and the n < N decline */
void Y(check_nfast_deconv_1d_general)(void);
/* rank-2 DECONV across a type-II axis and an odd axis */
void Y(check_nfast_deconv_2d_general)(void);
/* rank-3 DECONV across type-II and odd axes */
void Y(check_nfast_deconv_3d_general)(void);
/* rank-4 DECONV across mixed type-II and odd axes */
void Y(check_nfast_deconv_nd_general)(void);
/* The CONV solver -- Step C, the node convolution
 * (matrix B).  Plans a 1D KB conv problem via the process-global planner and
 * asserts output values directly.  Forward: single grid spike -> f[j] equals
 * a direct recompute via Y(window_phi) + the wrap formula.  Adjoint (B^H):
 * single node spike -> scatters psi weights onto its wrapped neighbors,
 * asserted directly from a clean grid (not merely a round-trip). */
void Y(check_nfast_conv_solver)(void);
/* The composed planner-native fast NFFT
 * solver (DECONV child + FFTW + CONV child, 1D/KB).  Plans a 1D even-N case
 * via plan_ng_guru with NFFT_NO_DIRECT so the native fast is the sole
 * survivor; asserts forward accuracy against the reference file, and
 * adjoint accuracy against the reference file (NOT a
 * round-trip -- proves no accidental 1/n normalization crept in). */
void Y(check_nfast_native_fast_accuracy)(void);
/* Two-flag isolation tests -- prove NFFT_NO_FAST_NATIVE
 * selects between the native fast and the direct fallback. */
void Y(check_nfast_native_tree)(void);
void Y(check_nfast_native_declines_window)(void);
void Y(check_nfast_flag_selective)(void);
/* The d=2 vertical slice -- deconv-2d + conv-2d leaf
 * solvers let the composed native-fast solver serve rank-2
 * problems.  Plans each of the 8 reference 2D cases via plan_ng_guru with
 * NFFT_NO_DIRECT so the native fast is the sole survivor; asserts forward/
 * adjoint accuracy against both the
 * reference file and an in-test legacy X(trafo_2d)/X(adjoint_2d) plan, and
 * the printed plan tree naming the composed solver and both children. */
void Y(check_nfast_native_fast_2d)(void);
void Y(check_nfast_native_fast_2d_adjoint)(void);
/* The d=3 vertical slice -- deconv-3d + conv-3d leaf
 * solvers let the composed native-fast solver serve rank-3
 * problems.  Same structure as the d=2 pair above, single reference
 * geometry N={10,10,10}, M=10. */
void Y(check_nfast_native_fast_3d)(void);
void Y(check_nfast_native_fast_3d_adjoint)(void);
/* The generic rnk>=4 slice -- deconv-nd + conv-nd (the
 * legacy carry-odometer leaves) let the composed native-fast solver
 * serve rank-4 (and higher) problems.  In-test-generated d=4 geometry (no
 * reference file exists for d>=4); oracle is an in-test legacy X(trafo)/
 * X(adjoint) plan (the generic legacy path IS the odometer this ports),
 * plus a direct-NDFT native cross-check. */
void Y(check_nfast_native_fast_4d)(void);
void Y(check_nfast_native_fast_4d_adjoint)(void);
/* Robust log(I0) and scaled I0 helpers -- log I0 matches
 * LOG(bessel_i0) in the small branch, matches the asymptotic expansion for
 * large x, and bessel_i0_scaled(x, log I0(peak)) stays in (0, 1] and finite
 * without ever materializing I0(x) in the asymptotic region. */
void Y(check_nfast_bessel_log_scaled)(void);
/* Native fast NFFT vs the true/direct 3D reference
 * transform (forward AND adjoint), checked against a precision-aware bound
 * (nfast_err_bound) in ALL precisions -- the single-precision regression gate
 * for the self-normalizing KB window fix (~1.0 in float before the fix,
 * ~2-3e-6 after). */
void Y(check_nfast_float_accuracy)(void);
/* Near the peak the rationalized KB evaluator beats
 * the naive exp(a - lg_peak) fold's ~3e-6 float residual -- checks
 * window_phi_hut against a stable log-domain reference, plus the
 * bessel_i0_logtail identity log I0(x) == x + logtail(x). */
void Y(check_nfast_window_cancellation)(void);
/* The window-agnostic range-apply API hoists the per-axis window
 * constant (KB: lg_peak) across a whole axis's evaluations. phi_hut_apply
 * fills the band bit-identically to the single-arg window_phi_hut; phi_precompute
 * builds the per-node psi bit-identically to single-arg window_phi. */
void Y(check_nfast_window_apply)(void);
/* The composed native fast (and its DECONV/CONV children) accept
 * KB/Gaussian/B-spline/sinc at runtime; Dirac (removed from the new API) and
 * out-of-range ordinals are declined -> NULL plan. */
void Y(check_nfast_native_window_select)(void);
/* Native fast NFFT (runtime window) vs the planner-native DIRECT
 * NDFT (window-independent exact oracle), forward + adjoint, for every
 * window selected via NFAST_WINDOWS (default: all four implemented
 * ordinals) -- one build validates KB/Gaussian/B-spline/sinc accuracy. */
void Y(check_nfast_window_accuracy)(void);

#endif
