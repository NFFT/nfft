/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "infft.h"
#include "iplanner.h"

struct Y(plan_ng_s);

void X(check_log2i)(void);
void X(check_next_power_of_2)(void);

/* Read a tests/refgen reference case file: d, N[d], M, x[d*M] node-major
 * (x[j*d+t]), f_hat[prod N] then f[M] as "re im" pairs. On success returns 1
 * and four arrays the caller releases with Y(free); on failure returns 0 with
 * every out pointer NULL. */
int Y(test_read_case)(const char *rel, int *d, INT **N, INT *NN, INT *M, R **x,
                      C **f_hat, C **f);

/* max|a - b| over max|b|, or max|a - b| when b is zero. */
R Y(test_rel_max_err)(const C *a, const C *b, INT len);

/* Relative accuracy bound for the fast NFFT pipeline, per window. */
R Y(test_err_bound)(int window, R m, R s);

/* Assert the printed plan tree of p contains needle. */
void Y(test_assert_plan_names)(struct Y(plan_ng_s) *p, const char *needle);

/* Wisdom round-trip through a string. The export is malloc'd; release it with
 * Y(free). */
char *Y(test_wisdom_export)(planner *pl);
int Y(test_wisdom_import)(planner *pl, const char *s);
