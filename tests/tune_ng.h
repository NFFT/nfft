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

/* Y(tune_plan_dyadic) returns one of the three rungs
 * n = 2^j * next_power_of_2(N), and a real transform at the pair meets the
 * goal. */
void Y(check_tune_dyadic_plan)(void);
/* The rung returned is the cheapest of those reaching the goal. Rung 1 is the
 * legacy grid, so this also pins that the answer is never rated dearer than
 * the legacy choice. */
void Y(check_tune_dyadic_cost)(void);
/* Goals below the reachable floor are capped rather than refused, illegal
 * rungs and bad arguments are refused, and a looser goal never costs more. */
void Y(check_tune_dyadic_capped)(void);

#endif
