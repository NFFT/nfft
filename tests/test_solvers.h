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

#ifndef TEST_SOLVERS_H
#define TEST_SOLVERS_H

#include "infft.h"
#include "iplanner.h"

/* Two test-only NFFT solvers (tests/test_solvers.c), registered only through
 * the functions below, never by the library roster. */

/* The permuting solver.  A coreless plan that borrows problem_nfft.x: awaking
 * to AWAKE_ZERO reverses the M node-blocks of x in place, returning to SLEEPY
 * reverses them back.  apply/apply_adjoint are no-ops -- the plan exists to
 * exercise the wakefulness and restore machinery, not to compute anything.
 * Set Y(nfft_perm_break_restore) to skip the restore, so a test can prove the
 * md5 restore guard actually fires. */
extern int Y(nfft_perm_break_restore);
plan *Y(nfft_perm_test_mkplan)(const problem *p);
void Y(nfft_solver_perm_test_register)(planner *pl);

/* The slow solver.  Accepts every NFFT problem with pcost 1e18, so it always
 * joins the candidate set and must always be pruned before timing.
 * Y(nfft_slow_test_applies) counts how often it was actually run: a correct
 * planner leaves it at zero. */
extern INT Y(nfft_slow_test_applies);
void Y(nfft_solver_slow_test_register)(planner *pl);

#endif
