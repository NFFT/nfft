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

#ifndef TUNE_NG_TEST_H
#define TUNE_NG_TEST_H

#include "infft.h"

/* The m that Y(tune) returns for a met goal really delivers that accuracy:
 * a Kaiser-Bessel plan_ng is built at that m and its error, in the same
 * measure tune predicts, is compared against the goal. */
void Y(check_tune_meets_goal)(void);
/* A goal below the minimum of the error curve is reported as not met, with
 * the argmin m and the attainable error handed back. */
void Y(check_tune_unreachable)(void);
/* Shape of the answer over the geometry grid: m within 1..n/2-1, monotone
 * in the goal, and non-increasing as the oversampling grows. */
void Y(check_tune_geometries)(void);
/* Invalid arguments are refused with -1 and leave the outputs alone. */
void Y(check_tune_bad_args)(void);
/* The n that Y(tune_sigma) returns admits a solution: Y(tune) reports the
 * goal met there, and one grid step below it does not. */
void Y(check_tune_sigma_agrees)(void);
/* A goal beyond reach at sigma = 4 is reported as such, and a real transform
 * at the recommended n meets a reachable goal. */
void Y(check_tune_sigma_limits)(void);
/* The (n, m) pair Y(tune_plan) hands back is self-consistent: n is even and
 * 5-smooth, Y(tune) agrees the goal is met at that n, and a real Kaiser-Bessel
 * transform run at the pair lands inside the goal. */
void Y(check_tune_plan)(void);
/* Goals below the reachable floor are capped rather than refused, and the
 * capped answer is still honoured. */
void Y(check_tune_plan_capped)(void);

/* The pair is the cheapest of those that reach the goal. */
void Y(check_tune_plan_cost)(void);

/* The measured refinement only ever takes cut-offs off, and never the one the
 * goal needs. */
void Y(check_tune_refine)(void);

#endif
