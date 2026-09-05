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

void Y(check_nfast_window_id)(void);

/* The DECONV and CONV problem types carry the attributes their stage needs
 * and omit the rest; the key discriminates geometry and stays data-blind. */
void Y(check_nfast_deconv_problem)(void);
void Y(check_nfast_conv_problem)(void);

/* The runtime window vtable: positivity and monotone decay, exact
 * normalization, and an independent recompute of the infft.h macro math for
 * every implemented window. */
void Y(check_nfast_window_vtable)(void);
void Y(check_nfast_window_normalized)(void);
void Y(check_nfast_window_all)(void);

/* A Kaiser-Bessel and a Gaussian NFFT problem hash to different wisdom
 * keys. */
void Y(check_nfast_window_key)(void);

/* The DECONV solver: forward, the adjoint of the real diagonal, and the
 * type-II frequency shift, each asserted as values from a clean input. */
void Y(check_nfast_deconv_solver)(void);

/* DECONV per rank, across odd N, type-II axes and the n < N decline. */
void Y(check_nfast_deconv_1d_general)(void);
void Y(check_nfast_deconv_2d_general)(void);
void Y(check_nfast_deconv_3d_general)(void);
void Y(check_nfast_deconv_nd_general)(void);

/* The CONV solver: forward and adjoint against a recompute from
 * Y(window_phi) and the wrap formula. */
void Y(check_nfast_conv_solver)(void);

/* The composed native fast NFFT on the 1D reference case, each direction
 * against its own reference file. */
void Y(check_nfast_native_fast_accuracy)(void);

/* NFFT_NO_FAST_NATIVE selects between the composed native fast and the
 * direct NDFT. */
void Y(check_nfast_native_tree)(void);
void Y(check_nfast_native_declines_window)(void);
void Y(check_nfast_flag_selective)(void);

/* The composed native fast at rank 2 and 3: the reference cases, checked
 * against the file and against the rank-matched in-test legacy X(trafo_2d) /
 * X(adjoint_3d) plan. The printed plan tree must name both children. */
void Y(check_nfast_native_fast_2d)(void);
void Y(check_nfast_native_fast_2d_adjoint)(void);
void Y(check_nfast_native_fast_3d)(void);
void Y(check_nfast_native_fast_3d_adjoint)(void);

/* The generic rank >= 4 slice. tests/refgen writes no d >= 4 reference file,
 * so the geometry is built in-test; the oracles are an in-test legacy
 * X(trafo)/X(adjoint) plan and the direct NDFT native. */
void Y(check_nfast_native_fast_4d)(void);
void Y(check_nfast_native_fast_4d_adjoint)(void);

/* The 3D reference case forward and adjoint, in whatever precision is
 * compiled, against Y(test_err_bound). */
void Y(check_nfast_float_accuracy)(void);

/* Near the peak the Kaiser-Bessel evaluator matches a stable log-domain
 * reference. */
void Y(check_nfast_window_cancellation)(void);

/* The range-apply window API matches the single-argument form to a relative
 * tolerance: phi_hut_apply over a whole band, phi_precompute over the
 * per-node psi taps. */
void Y(check_nfast_window_apply)(void);

/* The composed native fast and its children accept KB, Gaussian, B-spline
 * and sinc at runtime; Dirac and out-of-range ordinals give a NULL plan. */
void Y(check_nfast_native_window_select)(void);

/* Each window selected by NFAST_WINDOWS, native fast against the
 * window-independent direct NDFT, forward and adjoint. */
void Y(check_nfast_window_accuracy)(void);

#endif
