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

#ifndef NPLAN_TEST_H
#define NPLAN_TEST_H

#include "infft.h"

void Y(check_nplan_problem)(void);

/* The problem carries the caller's x/f_hat/f arrays (FFTW-parity). */
void Y(check_nplan_problem_owns_data)(void);

/* The wisdom key stays data-blind (pointers excluded). */
void Y(check_nplan_key_is_data_blind)(void);

/* Unit axes are elided at construction on the top-level (copy_x=1) path; the
 * borrowed path keeps full rank by design, and elision is drop-only (no
 * canonical sort). */
void Y(check_nplan_elides_unit_axes_geometry)(void);

/* The per-axis type-II variant participates in the wisdom key
 * (even N only; odd N normalizes to type-I). */
void Y(check_nplan_variant_key)(void);

/* The fast native serves odd N and type-II axes. At sigma == 1 it declines
 * and the direct NDFT wins; with NFFT_NO_DIRECT and no applicable fast plan
 * the guru returns NULL. */
void Y(check_nplan_guard_declines)(void);
void Y(check_nplan_solvers)(void);

/* For d == 1 the planner-native 1D NDFT solver wins. */
void Y(check_nplan_ndft_dispatch)(void);

/* For d >= 2 the single generic native multivariate NDFT solver wins. */
void Y(check_nplan_ndft_multivariate_dispatch)(void);

/* Accuracy of the native NDFT against a long-double, argument-reduced
 * reference (large N, nodes near +-1/2). */
void Y(check_nplan_ndft_accuracy)(void);
void Y(check_nplan_correct)(void);
void Y(check_nplan_wisdom_memo)(void);
void Y(check_nplan_measured)(void);
void Y(check_nplan_destructive_default)(void);
void Y(check_nplan_execute_on)(void);
void Y(check_nplan_print_includes_registrar_names)(void);
void Y(check_nplan_measured_wisdom)(void);
void Y(check_nplan_measured_prunes_by_estimate)(void);
void Y(check_nplan_timelimit_tight_degrades_to_estimate)(void);
void Y(check_nplan_timelimit_unset_measures_and_blesses)(void);
void Y(check_nplan_timelimit_partial_race_does_not_bless)(void);

/* An internal FFTW plan the caller's fftw_flags cannot produce makes the fast
 * solver decline. */
void Y(check_nplan_fftw_wisdom_only_declines)(void);
void Y(check_nplan_set_timelimit_roundtrip)(void);
void Y(check_nplan_public_api)(void);

/* Forward-only race, dual apply/apply_adjoint: the bundle races only
 * the forward problem (prob[ADJ] stays NULL); the single winning plan carries
 * both apply and apply_adjoint, and execute_adjoint dispatches to the same
 * plan's apply_adjoint. */
void Y(check_nplan_forward_only_race)(void);
void Y(check_nplan_apply_adjoint)(void);

/* plan_ng_guru copies the caller's x at guru time (f_hat/f are still aliased);
 * a post-guru mutation of the caller's array must not be visible to the plan. */
void Y(check_nplan_x_copied_not_aliased)(void);

/* guru, precompute and execute on a 1D fast plan. */
void Y(check_nplan_precompute_execute)(void);
void Y(check_nplan_destroy_keeps_caller_arrays)(void);

/* A bundle whose winner is the direct NDFT precomputes and executes
 * correctly, checked against trafo_direct. */
void Y(check_nplan_direct_only_bundle)(void);

/* plan_ng_guru accepts a per-axis NDFT variant array (NULL = all
 * type-I) and builds an executable plan from it. */
void Y(check_nplan_variant_guru)(void);

/* Odd per-axis N for the nD native (kernel/nfft/ndft-nd.c carry loop),
 * including non-outermost odd axes and odd+unit mixes. */
void Y(check_nplan_odd_n)(void);

/* The type-II 1D native computes the ascending (+N/2-at-last-slot) range:
 * uniform +1 shift of the type-I range. */
void Y(check_nplan_type_ii_1d)(void);

/* per-axis type-II in the nD native (kernel/nfft/ndft-nd.c carry loop):
 * a mixed type-I/type-II axis problem, uniform +1 shift on type-II axes only. */
void Y(check_nplan_type_ii_nd)(void);

/* data-driven acceptance tests for the native NDFT
 * solvers against the tests/refgen-generated reference data (roster:
 * 66 reused type-I even cases + 34 new odd/type-II cases, both
 * trafo and adjoint). */
void Y(check_nplan_data)(void);

/* The 3-state wakefulness enum (SLEEPY < AWAKE_ZERO < AWAKE). */
void Y(check_nplan_awake_states)(void);

/* mkproblem_nfft copies x for the top-level problem -- the
 * library never writes the user's x, and the user's array is independent of
 * the plan (mutating it post-plan does not disturb execute). */
void Y(check_nplan_user_x_pristine)(void);

/* The test-only raceable permuting AWAKE_ZERO
 * solver (tests/nplan_perm.c) proves the plan-level AWAKE_ZERO mutate/
 * restore contract end-to-end -- direct-drive restore, the restore guard
 * predicate catching a broken restore, the guard passing inside a real
 * measured race with >= 2 candidates, and AWAKE_ZERO being an internal-only
 * state (a fresh winner is SLEEPY; precompute awakens it to AWAKE). */
void Y(check_nplan_awake_zero_restore)(void);
void Y(check_nplan_restore_guard_fires)(void);
void Y(check_nplan_in_race_guard_passes)(void);
void Y(check_nplan_awake_zero_internal)(void);

/* A unit axis (N_t==1) is elided at construction, so the surviving
 * axes go through the fast algorithm (NFFT_NO_DIRECT builds only if elision
 * happened). */
void Y(check_nplan_elides_unit_axes)(void);

/* All axes unit -> rank-0 exact base-case solver nfft_solver_const_0d
 * (forward broadcast, adjoint reduce). Ungated: plans under NFFT_NO_DIRECT. */
void Y(check_nplan_rank0_solver)(void);

/* Position-agnostic correctness -- four unit-axis shapes
 * (leading/middle/trailing/two-unit), forward AND adjoint compared against an
 * in-test direct NDFT oracle over the full d-tensor under the window-aware
 * bound, plus an NFFT_NO_DIRECT non-NULL assertion per shape (the fast path
 * only engages if the unit axis was elided). */
void Y(check_nplan_unit_axis_correct)(void);

/* New-array execute (execute_on/execute_adjoint_on) on a native-fast plan:
 * the swapped problem pointers must reach the DECONV and CONV children. */
void Y(check_nplan_newarray_native_fast)(void);

/* New-array execute (execute_on) on a unit-axis plan -- a
 * full-d-layout replacement f_hat must read correctly through the compressed
 * (unit-axis-elided) strides; compared to the direct NDFT oracle. */
void Y(check_nplan_unit_axis_execute_on)(void);

/* The guru must return NULL (FFTW's "decline to plan" contract)
 * on NULL/zero required args in RELEASE builds, where A() is a no-op. */
void Y(check_nplan_guru_rejects_null_args)(void);

/* The guru returns NULL on non-positive per-axis geometry (N[t]<=0 /
 * n[t]<=0) and on a bad M, m or window ordinal. */
void Y(check_nplan_guru_rejects_bad_geometry)(void);

/* Entry lines in an exported wisdom file are indented; the preamble and the
 * closing paren are not. Declared here because the wisdom-only case reuses
 * it. */
int Y(count_wisdom_entries)(void);

/* After planning, an export must describe the whole winning tree: the NFFT
 * solution and its DECONV and CONV children. */
void Y(check_nplan_blesses_whole_tree)(void);

/* Wisdom-only answers from the store or fails outright, at every patience
 * level and with children included; a miss disturbs nothing and never
 * poisons the next ordinary call. */
void Y(check_nplan_wisdom_only)(void);

#endif
