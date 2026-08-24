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

#ifndef FAST_NATIVE_TEST_H
#define FAST_NATIVE_TEST_H

#include "infft.h"

/* Gate-flag isolation: NFFT_NO_FAST_NATIVE / NFFT_NO_DIRECT must select
 * between the composed native fast and the direct fallback, and the printed
 * plan tree must name the composed solver and both children. */
void Y(check_fast_native_tree)(void);
void Y(check_fast_native_declines_window)(void);
void Y(check_fast_native_flag_selective)(void);
/* Per-rank cross-check against the *legacy* X(trafo)/X(adjoint) implementation
 * of the same algorithm -- the one oracle tests/nfft_ng.c cannot provide, since
 * it never touches legacy code.  d=2 and d=3 run the reference geometries;
 * d>=4 has no reference file, so the geometry is generated in-test and the
 * direct NDFT is cross-checked as a second oracle. */
void Y(check_fast_native_2d)(void);
void Y(check_fast_native_2d_adjoint)(void);
void Y(check_fast_native_3d)(void);
void Y(check_fast_native_3d_adjoint)(void);
void Y(check_fast_native_4d)(void);
void Y(check_fast_native_4d_adjoint)(void);
/* The composed solver and its children accept KB/Gaussian/B-spline/sinc at
 * runtime; Dirac and out-of-range ordinals are declined -> NULL plan. */
void Y(check_fast_native_window_select)(void);
/* Native fast vs the planner-native direct NDFT (a window-independent exact
 * oracle), forward and adjoint, for every window selected via NFAST_WINDOWS
 * (default: all four implemented ordinals). */
void Y(check_fast_native_window_accuracy)(void);

#endif
