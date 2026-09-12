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

/* Test-only introspection hooks for the opaque plan_ng handle.  Kept OUT of
 * iplanner.h so the internal planner header carries no test scaffolding; the
 * implementation lives in kernel/nfft/plan.c (the only place with the
 * plan_ng_s struct body).  Included by both that definer and the acceptance
 * tests (tests/nplan.c) so the prototype has a single home. */

#ifndef PLAN_NG_TEST_H
#define PLAN_NG_TEST_H

#include "infft.h" /* Y() name mangling */
#include "iplanner.h" /* opcnt */

struct Y(plan_ng_s); /* opaque handle */

/* The winning forward plan's wakefulness: SLEEPY right after construction,
 * AWAKE after precompute. */
int Y(plan_ng_test_awake_state)(const struct Y(plan_ng_s) *p);

/* Fills the winning forward plan's own counts and the three parts they are
 * composed from. Returns 0 when the winner is not the native fast solver, so a
 * caller can skip rather than assert against a plan that has no children. */
int Y(plan_ng_test_native_parts)(const struct Y(plan_ng_s) *p, opcnt *top,
                              opcnt *deconv, opcnt *conv, opcnt *fft);

#endif
