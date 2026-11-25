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

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "nnfft.h"
#include "accuracy_log.h"
#include "util.h"

#define SEED 1234567890L

#if defined(KAISER_BESSEL)

typedef struct testcase_s
{
  int d;
  int N;   /* bandwidth */
  int N1;  /* oversampled bandwidth, sigma = N1 / N */
  int NV;  /* number of frequency nodes v */
  int M;   /* number of spatial nodes x */
} testcase_t;

/* NNFFT wraps an inner NFFT on the aN1 grid in a window convolution in v, so
 * the bound adds two NFFT bounds, outer with sigma = N1 / N and inner with the
 * inner plan's parameters, and the direct transform's round-off. */
static R err_trafo(nnfft_plan *p)
{
  INT N[p->d], N1[p->d];
  int t;

  for (t = 0; t < p->d; t++)
  {
    N[t] = p->N[t];
    N1[t] = p->N1[t];
  }

  return Y(kb_err_trafo)(p->d, p->m, p->sigma, N, N1, 0)
      + Y(kb_err_trafo)(p->direct_plan->d, p->direct_plan->m,
            p->direct_plan->sigma, p->direct_plan->N, p->direct_plan->n, 0)
      + K(56.0) * Y(float_property)(NFFT_EPSILON);
}

/* A unit coefficient at a node just below v = 1/2, whose run reaches the far
 * end of the inner grid, or a unit sample at one node for the adjoint. With
 * l1 norm 1 the error is the sup norm the bound describes. */
static int check_single(const testcase_t *tc, const int adjoint)
{
  nnfft_plan p;
  int j, t, ok, N[tc->d], N1[tc->d];
  C *ref;
  R err = K(0.0), bound;
  const int len = adjoint ? tc->NV : tc->M;

  for (t = 0; t < tc->d; t++)
  {
    N[t] = tc->N;
    N1[t] = tc->N1;
  }

  nnfft_init_guru(&p, tc->d, tc->NV, tc->M, N, N1, WINDOW_HELP_ESTIMATE_m,
      PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_V | MALLOC_F | MALLOC_F_HAT);

  printf("%-31s d = %-1d, N = %-5d, N1 = %-5d, NV = %-5d, M = %-5d, m = %2d",
      adjoint ? "nnfft_adjoint_online_edge" : "nnfft_online_edge",
      tc->d, tc->N, tc->N1, tc->NV, tc->M, p.m);

  Y(srand48)(SEED);
  for (j = 0; j < tc->M * tc->d; j++)
    p.x[j] = Y(drand48)() - K(0.5);
  for (j = 0; j < tc->NV * tc->d; j++)
    p.v[j] = Y(drand48)() - K(0.5);
  for (t = 0; t < tc->d; t++)
    p.v[t] = K(0.5) - K(1.0) / K(1048576.0);

  nnfft_precompute_one_psi(&p);

  ref = Y(malloc)((size_t)len * sizeof(C));

  if (adjoint)
  {
    for (j = 0; j < tc->M; j++)
      p.f[j] = IF(j == 0, K(1.0), K(0.0));
    nnfft_adjoint_direct(&p);
    for (j = 0; j < len; j++)
      ref[j] = p.f_hat[j];
    nnfft_adjoint(&p);
    for (j = 0; j < len; j++)
      err = FMAX(err, CABS(p.f_hat[j] - ref[j]));
  }
  else
  {
    for (j = 0; j < tc->NV; j++)
      p.f_hat[j] = IF(j == 0, K(1.0), K(0.0));
    nnfft_trafo_direct(&p);
    for (j = 0; j < len; j++)
      ref[j] = p.f[j];
    nnfft_trafo(&p);
    for (j = 0; j < len; j++)
      err = FMAX(err, CABS(p.f[j] - ref[j]));
  }

  bound = err_trafo(&p);
  ok = IF(err < bound, 1, 0);
  printf(" -> %-4s " __FE__ " (" __FE__ ")\n", IF(ok == 0, "FAIL", "OK"), err, bound);
  accuracy_log_append("nnfft", "online_edge", tc->d, N, tc->M, "init_guru (PRE PSI)",
      adjoint ? "adjoint" : "trafo", (long double)err, (long double)bound, ok);

  Y(free)(ref);
  nnfft_finalize(&p);
  return ok;
}

static const testcase_t testcases_1d_online_edge[] =
{
  {1, 32, 64, 50, 50},
};

static void check_many(const testcase_t *tcs, const size_t n, const int adjoint)
{
  size_t i;
  int ok = 1, r;

  for (i = 0; i < n; i++)
  {
    r = check_single(&tcs[i], adjoint);
    ok = MIN(ok, r);
  }
  CU_ASSERT(ok);
}

void nnfft_check_1d_online_edge(void)
{
  printf("check_1d_online_edge:\n");
  check_many(testcases_1d_online_edge, SIZE(testcases_1d_online_edge), 0);
}

void nnfft_check_adjoint_1d_online_edge(void)
{
  printf("check_adjoint_1d_online_edge:\n");
  check_many(testcases_1d_online_edge, SIZE(testcases_1d_online_edge), 1);
}

#endif
