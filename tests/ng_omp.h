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

#ifndef NG_OMP_TEST_H
#define NG_OMP_TEST_H

#include "infft.h"

/* The add-on library's roster contract now that it carries a threaded solver. */
void Y(check_ng_omp_roster)(void);

/* The threaded DECONV solver agrees with the serial solver of the same rank. */
void Y(check_ng_omp_deconv_1d)(void);

/* Ranks 2, 3 and 4 against the serial solver of the same rank. */
void Y(check_ng_omp_deconv_2d)(void);
void Y(check_ng_omp_deconv_3d)(void);
void Y(check_ng_omp_deconv_4d)(void);

/* An odd axis and a type-II axis both move the slot split. */
void Y(check_ng_omp_deconv_variants)(void);

/* The results do not move with the team size. */
void Y(check_ng_omp_thread_count_invariance)(void);

/* A one-thread budget makes the threaded solver decline. */
void Y(check_ng_omp_declines_serial_budget)(void);

/* The threaded CONV forward reproduces the serial one to the last bit. */
void Y(check_ng_omp_conv_forward)(void);

/* The threaded CONV adjoint reproduces the serial one to the last bit. */
void Y(check_ng_omp_conv_adjoint)(void);

/* A grid too narrow to band (nb == 1) still matches the serial solver. */
void Y(check_ng_omp_conv_narrow_grid)(void);

/* Nodes clustered into a slice of the grid still match the serial solver. */
void Y(check_ng_omp_conv_clustered)(void);

/* A band count that does not divide the grid, wrapping at the boundary. */
void Y(check_ng_omp_conv_wrap)(void);

/* A measured CONV race holds one candidate's tables at a time. */
void Y(check_ng_omp_race_memory)(void);

/* With both rosters live and two threads, a child problem's two applicable
 * solvers are raced by measurement. */
void Y(check_ng_omp_child_race)(void);

/* End to end through the public API with the count raised, at every patience
 * level: both threaded children answer, and results match one thread. */
void Y(check_ng_omp_nfft_all_levels)(void);

/* NFFT_NO_FAST_NATIVE leaves nothing to plan below NFFT_PATIENT. */
void Y(check_ng_omp_no_fast_native_refuses)(void);

/* The sorted threaded CONV solver declines on presorted nodes. */
void Y(check_ng_omp_conv_sorted_declines)(void);

#endif /* NG_OMP_TEST_H */
