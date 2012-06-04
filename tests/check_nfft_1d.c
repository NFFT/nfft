/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: simple_test.c 3509 2010-05-25 19:00:59Z keiner $ */

#include <CUnit/CUnit.h>

#include "check_nfft.h"
#include "infft.h"

static const char* filenames[] =
{
  "data/nfft_1d_1_1.txt",
  /*"data/nfft_1d_1_10.txt",
  "data/nfft_1d_1_20.txt",
  "data/nfft_1d_1_50.txt",
  "data/nfft_1d_2_1.txt",
  "data/nfft_1d_2_10.txt",
  "data/nfft_1d_2_20.txt",
  "data/nfft_1d_2_50.txt",
  "data/nfft_1d_4_1.txt",
  "data/nfft_1d_4_10.txt",
  "data/nfft_1d_4_20.txt",
  "data/nfft_1d_4_50.txt",
  "data/nfft_1d_10_1.txt",
  "data/nfft_1d_10_10.txt",
  "data/nfft_1d_10_20.txt",
  "data/nfft_1d_10_50.txt",
  "data/nfft_1d_20_1.txt",
  "data/nfft_1d_20_10.txt",
  "data/nfft_1d_20_20.txt",
  "data/nfft_1d_20_50.txt",
  "data/nfft_1d_50_1.txt",
  "data/nfft_1d_50_10.txt",
  "data/nfft_1d_50_20.txt",
  "data/nfft_1d_50_50.txt",*/
};

static const testcase_delegate_file_t nfft_1d_1_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_1.txt"};
static const testcase_delegate_file_t nfft_1d_1_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_10.txt"};
static const testcase_delegate_file_t nfft_1d_1_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_20.txt"};
static const testcase_delegate_file_t nfft_1d_1_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_50.txt"};
static const testcase_delegate_file_t nfft_1d_2_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_1.txt"};
static const testcase_delegate_file_t nfft_1d_2_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_10.txt"};
static const testcase_delegate_file_t nfft_1d_2_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_20.txt"};
static const testcase_delegate_file_t nfft_1d_2_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_50.txt"};
static const testcase_delegate_file_t nfft_1d_4_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_1.txt"};
static const testcase_delegate_file_t nfft_1d_4_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_10.txt"};
static const testcase_delegate_file_t nfft_1d_4_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_20.txt"};
static const testcase_delegate_file_t nfft_1d_4_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_50.txt"};
static const testcase_delegate_file_t nfft_1d_10_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_1.txt"};
static const testcase_delegate_file_t nfft_1d_10_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_10.txt"};
static const testcase_delegate_file_t nfft_1d_10_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_20.txt"};
static const testcase_delegate_file_t nfft_1d_10_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_50.txt"};
static const testcase_delegate_file_t nfft_1d_20_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_1.txt"};
static const testcase_delegate_file_t nfft_1d_20_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_10.txt"};
static const testcase_delegate_file_t nfft_1d_20_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_20.txt"};
static const testcase_delegate_file_t nfft_1d_20_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_50.txt"};
static const testcase_delegate_file_t nfft_1d_50_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_1.txt"};
static const testcase_delegate_file_t nfft_1d_50_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_10.txt"};
static const testcase_delegate_file_t nfft_1d_50_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_20.txt"};
static const testcase_delegate_file_t nfft_1d_50_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_50.txt"};
static const testcase_delegate_file_t *testcases[] = {
  &nfft_1d_1_1, &nfft_1d_1_10, &nfft_1d_1_20, &nfft_1d_1_50,
  &nfft_1d_2_1, &nfft_1d_2_10, &nfft_1d_2_20, &nfft_1d_2_50,
  &nfft_1d_4_1, &nfft_1d_4_10, &nfft_1d_4_20, &nfft_1d_4_50,
  &nfft_1d_10_1, &nfft_1d_10_10, &nfft_1d_10_20, &nfft_1d_10_50,
  &nfft_1d_20_1, &nfft_1d_20_10, &nfft_1d_20_20, &nfft_1d_20_50,
  &nfft_1d_50_1, &nfft_1d_50_10, &nfft_1d_50_20, &nfft_1d_50_50,
};

/* Initializers. */
static void init_1d_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  X(init_1d)(p, N[0], M);
}

static void init_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  X(init)(p, d, N, M);
}

static void init_advanced_pre_psi_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  int *n = malloc(d*sizeof(int));
  int i;
  for (i = 0; i < d; i++)
    n[i] = 2*X(next_power_of_2)(N[i]);
  X(init_guru)(p, d, N, M, n, ego->m, ego->nfft_flags, ego->fftw_flags);
  free(n);
}

/*static void init_advanced_pre_full_psi_(X(plan) *p, const int d, const int *N, const int M)
{
  X(init_guru)(p, d, N, M, PRE_PHI_HUT | PRE_FULL_PSI | MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT | FFT_OUT_OF_PLACE, 0U);
}*/

#define DEFAULT_NFFT_FLAGS MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT | FFT_OUT_OF_PLACE
#define DEFAULT_FFTW_FLAGS FFTW_ESTIMATE | FFTW_DESTROY_INPUT

static init_delegate_t init_1d = {"init_1d", init_1d_};
static init_delegate_t init = {"init", init_};
static init_delegate_t init_advanced_pre_psi = {"init_guru (PRE PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
static init_delegate_t init_advanced_pre_full_psi = {"init_guru (PRE FULL PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_FULL_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
static init_delegate_t init_advanced_pre_lin_psi = {"init_guru (PRE LIN PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_LIN_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
#if defined(GAUSSIAN)
static init_delegate_t init_advanced_pre_fg_psi = {"init_guru (PRE FG PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | FG_PSI | PRE_FG_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
#endif
/* Transformations. */
static double trafo_direct_cost_factor = 1.0E-6;

static double trafo_direct_cost(X(plan) *p)
{
  if (trafo_direct_cost_factor == 0.0)
  {
    int M, d, Nd, x = 0;
    for (d = 1; d <= 4; d++)
    {
      for (Nd = 4; Nd < 128; Nd *= 2)
      {
        for (M = 4; M <= 128; M *= 2)
        {
          X(plan) p;
          int *N = malloc(d*sizeof(int)), i;
          for (i = 0; i < d; i++)
          {
            N[i] = Nd;
          }
          X(init)(&p, d, N, M);
          for (i = 0; i < M; i++)
            p.x[i] = K(0.0);
          if(p.nfft_flags & PRE_ONE_PSI)
            X(precompute_one_psi)(&p);
          for (i = 0; i < d*Nd; i++)
          {
            p.f_hat[i] = K(0.0) + K(0.0) * I;
          }
          {
            double r;
            ticks t0, t1;
            t0 = getticks();
            X(trafo_direct)(&p);
            t1 = getticks();
            r = X(elapsed_seconds)(t1, t0)/M;
            for (i = 0; i < d; i++)
              r = r / Nd;
            trafo_direct_cost_factor += r;
            printf("%E\n", r);
            x += 1;
          }
          X(finalize)(&p);
          free(N);
        }
      }
    }
    trafo_direct_cost_factor = trafo_direct_cost_factor/((double)x);
    printf("--> %E\n", trafo_direct_cost_factor);
  }

  {
    int c = p->M_total, i;

    for (i = 0; i < p->d; i++)
      c *= p->N[i];

    return trafo_direct_cost_factor * c;
  }
}

static trafo_delegate_t trafo_direct = {"trafo_direct", X(trafo_direct), 0, X(trafo_direct_cost), X(err_trafo_direct)};
static trafo_delegate_t trafo = {"trafo", X(trafo), X(check), 0, X(err_trafo)};

/* Delegates. */
static const init_delegate_t* initializers[] = {
    &init_1d,
    &init,
    &init_advanced_pre_psi,
    &init_advanced_pre_full_psi,
    &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
    &init_advanced_pre_fg_psi,
#endif
};
static const trafo_delegate_t* trafos[] = {&trafo_direct, &trafo};

static void check_nfft_1d_file()
{
  X(check_many)(SIZE(testcases), SIZE(initializers), SIZE(trafos), testcases,
    initializers, trafos);
}

int main(void)
{
  CU_pSuite s;
  CU_initialize_registry();
  CU_set_output_filename("nfft");
  s = CU_add_suite("nfft_1d", 0, 0);
  CU_add_test(s, "nfft_1d_file", check_nfft_1d_file);
  CU_automated_run_tests();
  /*CU_basic_run_tests();*/
  CU_cleanup_registry();
  return 0;
}
