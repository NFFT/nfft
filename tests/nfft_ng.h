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

#ifndef NFFT_NG_TEST_H
#define NFFT_NG_TEST_H

#include "infft.h"

/* Planner-API (plan_ng) NFFT accuracy suite: the tests/nfft.c roster with the
 * tests/nfft.c bounds, plus the per-axis type-I/type-II, odd-N and unit-axis
 * geometries only this API can express.  See tests/nfft_ng.c. */

/* Reference-file reader, shared with tests/nplan_data.c.  Reads d, N[d], M,
 * x[M*d], f_hat[NN], f[M] from a tests/refgen case file named relative to
 * tests/, at the build precision (__FI__) -- the files carry 64 significant
 * digits, and reading them through double would cap every case at double
 * round-off.  Returns 0 and leaves every out pointer NULL on failure; the
 * caller frees N/x/f_hat/f. */
int Y(nfft_ng_read_case)(const char *rel, int *d, INT **N, INT *NN, INT *M,
                         R **x, C **f_hat, C **f);

/* Reference-file cases per dimension: direct NDFT, composed fast NFFT, and the
 * planner's own choice. */
void Y(check_nfft_ng_1d_file)(void);
void Y(check_nfft_ng_2d_file)(void);
void Y(check_nfft_ng_3d_file)(void);
void Y(check_nfft_ng_4d_file)(void);
void Y(check_nfft_ng_adjoint_1d_file)(void);
void Y(check_nfft_ng_adjoint_2d_file)(void);
void Y(check_nfft_ng_adjoint_3d_file)(void);
void Y(check_nfft_ng_adjoint_4d_file)(void);

/* The tests/nfft.c online grid: fast NFFT against the direct NDFT. */
void Y(check_nfft_ng_1d_online)(void);
void Y(check_nfft_ng_2d_online)(void);
void Y(check_nfft_ng_3d_online)(void);
void Y(check_nfft_ng_4d_online)(void);
void Y(check_nfft_ng_adjoint_1d_online)(void);
void Y(check_nfft_ng_adjoint_2d_online)(void);
void Y(check_nfft_ng_adjoint_3d_online)(void);
void Y(check_nfft_ng_adjoint_4d_online)(void);

/* Fast NFFT at online sizes over per-axis type-I, type-II, odd and mixed
 * geometries.  Forward and adjoint, oracle is the direct NDFT. */
void Y(check_nfft_ng_fast_variants_1d)(void);
void Y(check_nfft_ng_fast_variants_2d)(void);
void Y(check_nfft_ng_fast_variants_3d)(void);
void Y(check_nfft_ng_fast_variants_4d)(void);

/* The same, with unit axes elided, down to the all-unit (rank 0) case. */
void Y(check_nfft_ng_fast_unit_axes_1d)(void);
void Y(check_nfft_ng_fast_unit_axes_2d)(void);
void Y(check_nfft_ng_fast_unit_axes_3d)(void);
void Y(check_nfft_ng_fast_unit_axes_4d)(void);

#endif /* NFFT_NG_TEST_H */
