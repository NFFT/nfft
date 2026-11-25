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

typedef struct testcase_s
{
  int d;
  int N;   /* bandwidth */
  int N1;  /* oversampled bandwidth, sigma = N1 / N */
  int NV;  /* number of frequency nodes v */
  int M;   /* number of spatial nodes x */
  R v0;    /* first frequency node, on every axis */
} testcase_t;

/* Per-axis truncation term and round-off amplification phi_hat(0) /
 * phi_hat(N/2) on grid n of one NFFT with oversampling sigma. Truncation
 * terms are those of the NFFT tests, with the actual sigma. */
#if defined(KAISER_BESSEL)
static R trunc_w(const int m, const R sigma)
{
  const R z = K2PI * ((R)m + K(1.0)) * SQRT(K(1.0) - K(1.0) / sigma);
  return K(3.5) * EXP(-z) / Y(bessel_i0_exp_scaled)(z);
}

static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R b = KPI * (K(2.0) - K(1.0) / sigma), r = (R)m + K(1.0);
  return Y(kb_phi_hut)(b, K(1.0), r, n, K(0.0)) / Y(kb_phi_hut)(b, K(1.0), r, n, N / K(2.0));
}
#elif defined(GAUSSIAN)
static R trunc_w(const int m, const R sigma)
{
  return K(0.41) * EXP(-(R)m * KPI * (K(1.0) - K(1.0) / (K(2.0) * sigma - K(1.0))));
}

static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R b = K(2.0) * sigma / (K(2.0) * sigma - K(1.0)) * ((R)m + K(1.0)) / KPI;
  const R x = KPI * N / (K(2.0) * n);
  return EXP(b * x * x);
}
#elif defined(B_SPLINE)
static R trunc_w(const int m, const R sigma)
{
  return K(12000.0) * POW(K(2.0) * sigma - K(1.0), -K(2.0) * (R)m);
}

static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R x = KPI * N / (K(2.0) * n);
  UNUSED(sigma);
  return POW(x / SIN(x), K(2.0) * (R)m + K(2.0));
}
#elif defined(SINC_POWER)
static R trunc_w(const int m, const R sigma)
{
  return K(4e-3) / ((R)m - K(1.0)) * (K(2.0) * POW(sigma, -K(2.0) * (R)m)
      + POW(sigma / (K(2.0) * sigma - K(1.0)), K(2.0) * (R)m));
}

/* phi_hat is the cardinal B-spline of order 2m, centred, scaled by w. */
static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R w = (sigma - SINCPOW_A) / ((R)m * sigma);
  return Y(bsplines)(2 * m, (R)m) / Y(bsplines)(2 * m, (R)m + N / (K(2.0) * w * n));
}
#endif

/* Round-off of the direct transform: it forms the phase v x N in double
 * without reduction, so the error grows with N. */
static R err_direct(const nnfft_plan *p)
{
  R sum_N = K(0.0);
  int t;

  for (t = 0; t < p->d; t++)
    sum_N += (R)p->N[t];

  return (K(56.0) + KPI / K(2.0) * sum_N) * Y(float_property)(NFFT_EPSILON);
}

/* NNFFT = outer window convolution onto the aN1 grid, inner NFFT on that
 * grid, deconvolution at x. Outer truncation plus inner truncation amplified
 * by the outer deconvolution (Kircheis, Potts, Tasche 2023, Thm 3.5, with
 * the sharper factor of Elbel, Steidl 1998), plus the round-off of both
 * deconvolutions and of the direct transform. */
static R err_trafo(const nnfft_plan *p)
{
  const nfft_plan *in = p->direct_plan;
  R t1 = K(0.0), t2 = K(0.0), a1 = K(1.0), a2 = K(1.0);
  int t;

  for (t = 0; t < p->d; t++)
  {
    const R s1 = (R)p->N1[t] / (R)p->N[t], s2 = (R)in->n[t] / (R)p->aN1[t];
    t1 += trunc_w(p->m, s1);
    t2 += trunc_w(p->m, s2);
    a1 *= amp_w(p->m, s1, (R)p->N[t], (R)p->N1[t]);
    a2 *= amp_w(p->m, s2, (R)p->aN1[t], (R)in->n[t]);
  }

  return t1 + a1 * t2 + K(3.5) * (a1 + a2) * Y(float_property)(NFFT_EPSILON)
      + err_direct(p);
}

/* The Gaussian edge error stays below the bound from m = 9 on, with or
 * without the full run on the grid. B-spline and sinc-power windows have
 * compact support and show no edge error at any m. */
#if defined(GAUSSIAN)
#define EDGE_M 4
#else
#define EDGE_M WINDOW_HELP_ESTIMATE_m
#endif

/* A unit coefficient at the first frequency node v0, at a band edge where its
 * run reaches the end of the inner grid, or a unit sample at one node for the
 * adjoint. With l1 norm 1 the error is the sup norm the bound describes. */
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

  nnfft_init_guru(&p, tc->d, tc->NV, tc->M, N, N1, EDGE_M,
      PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_V | MALLOC_F | MALLOC_F_HAT);

  printf("%-31s d = %-1d, N = %-5d, N1 = %-5d, NV = %-5d, M = %-5d, v0 = % .7f, m = %2d",
      adjoint ? "nnfft_adjoint_online_edge" : "nnfft_online_edge",
      tc->d, tc->N, tc->N1, tc->NV, tc->M, (double)tc->v0, p.m);

  Y(srand48)(SEED);
  for (j = 0; j < tc->M * tc->d; j++)
    p.x[j] = Y(drand48)() - K(0.5);
  for (j = 0; j < tc->NV * tc->d; j++)
    p.v[j] = Y(drand48)() - K(0.5);
  for (t = 0; t < tc->d; t++)
    p.v[t] = tc->v0;

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
  {1, 32, 64, 50, 50, K(0.5) - K(1.0) / K(1048576.0)},
  {1, 32, 64, 50, 50, K(-0.5)},
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
