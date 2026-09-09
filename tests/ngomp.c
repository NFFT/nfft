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

#include <complex.h> /* before nfft3.h so fftw_complex is C-compatible */
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "ngomp.h"

/* The add-on library's contract while its roster is empty. init_threads reports
 * failure, so plan_with_nthreads leaves the count at 1 and no caller can reach
 * a state with more than one thread and nothing to run on them. When the first
 * threaded solver is added to kernel/threads/conf.c these expectations invert,
 * deliberately and visibly. */
void Y(check_ngomp_empty_roster)(void)
{
  CU_ASSERT_EQUAL(Y(nfft_threads_roster_size)(), 0);
  CU_ASSERT_EQUAL(NFFT(init_threads)(), 0);

  NFFT(plan_with_nthreads)(4);
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);

  NFFT(cleanup_threads)(); /* safe when init never succeeded */
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);
}

/* Planning still works at every patience level through the add-on library. */
void Y(check_ngomp_plans_serially)(void)
{
  const INT N = 64, n = 128, M = 100;
  const unsigned levels[3] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  int i;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  NFFT(plan_with_nthreads)(4); /* refused; the count stays 1 */
  for (i = 0; i < 3; i++) {
    Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                       NFFT(get_window_id)(), x,
                              (FC *)f_hat, (FC *)f, 0u, levels[i]);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p) {
      NFFT(precompute)(p);
      NFFT(execute)(p);
      NFFT(plan_ng_destroy)(p);
    }
  }

  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
