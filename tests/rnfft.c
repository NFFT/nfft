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

/* The fast RNFFT against the direct one for every precomputation scheme, the
 * direct RNFFT against the complex NFFT on Hermitian spectra, the adjoint
 * identity of both pairs, and the c2r solver. */

#include "config.h"

#include <stdio.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"
#include "util.h"
#include "rnfft.h"
#include "accuracy_log.h"

#define SEED 1234567890L

#define DIRECT_ERR (K(56.0) * Y(float_property)(NFFT_EPSILON))

#define BASE_FLAGS (MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT)
#define FFTW_FLAGS (FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

typedef struct
{
  const char *name;
  int guru; /* 0: X(init); 1: X(init_guru), default m; 2: half m */
  unsigned flags;
} init_t;

static const init_t inits[] = {
  {"init", 0, 0},
  {"init_guru (PRE PSI)", 1, PRE_PHI_HUT | PRE_PSI},
  {"init_guru (PRE PSI, out of place)", 1,
    PRE_PHI_HUT | PRE_PSI | FFT_OUT_OF_PLACE},
  {"init_guru (PRE PSI, sorted)", 1, PRE_PHI_HUT | PRE_PSI | NFFT_SORT_NODES},
  {"init_guru (PRE FULL PSI)", 1, PRE_PHI_HUT | PRE_FULL_PSI},
  {"init_guru (PRE FULL PSI, blockwise)", 1,
    PRE_FULL_PSI | NFFT_OMP_BLOCKWISE_ADJOINT},
  /* The table for the default m takes GBs at d = 3. */
  {"init_guru (PRE LIN PSI, m=half)", 2, PRE_PHI_HUT | PRE_LIN_PSI},
  {"init_guru (no PRE)", 1, 0},
  {"init_guru (no PRE, blockwise)", 1, NFFT_OMP_BLOCKWISE_ADJOINT},
#if defined(GAUSSIAN)
  {"init_guru (PRE FG PSI)", 1, PRE_PHI_HUT | PRE_FG_PSI},
  {"init_guru (FG PSI, sorted)", 1, FG_PSI | NFFT_SORT_NODES},
#endif
};

/* Error bound of the fast transforms with sup norm over the l1 norm of the
 * weighted input. */
static R err_trafo(const X(plan) *p)
{
#if defined(KAISER_BESSEL)
  return Y(kb_err_trafo)(p->d, p->m, p->sigma, p->N, p->n, 0);
#else
  const R m = (R)p->m, eps = Y(float_property)(NFFT_EPSILON);
  R s = p->sigma[0], a, b, err;
  INT t;

  for (t = 1; t < p->d; t++)
    s = FMIN(s, p->sigma[t]);
#if defined(GAUSSIAN)
  a = K(0.41);
  b = MANT_DIG == 24 ? K(2000.0) : K(50.0);
  err = EXP(-m * KPI * (K(1.0) - K(1.0) / (K(2.0) * K(2.0) - K(1.0))));
#elif defined(B_SPLINE)
  a = K(1.0);
  b = K(2000.0);
  err = K(12000.0) * POW(K(1.0) / (K(2.0) * s - K(1.0)), K(2.0) * m);
#elif defined(SINC_POWER)
  a = K(1.0);
  b = K(2000.0);
  err = (K(1.0) / (m - K(1.0))) * ((K(2.0) / (POW(s, K(2.0) * m)))
    + POW(s / (K(2.0) * s - K(1.0)), K(2.0) * m));
#else
#error Unsupported window function.
#endif
  return FMAX(FMAX(a * err, b * eps), DIRECT_ERR);
#endif
}

/* Weight c_k of coefficient i: 1 at k_{d-1} in {0, N_{d-1}/2}, 2 otherwise. */
static R weight(const X(plan) *p, const INT i)
{
  const INT Nh = p->N[p->d - 1] / 2, l = i % (Nh + 1);

  return l == 0 || l == Nh ? K(1.0) : K(2.0);
}

static void init_plan(X(plan) *p, const init_t *ini, const INT d,
  const int *N, const int M)
{
  int n[4];
  INT t;

  if (!ini->guru)
  {
    X(init)(p, (int)d, (int *)N, M);
    return;
  }
  for (t = 0; t < d; t++)
    n[t] = 2 * (int)Y(next_power_of_2)(N[t]);
  X(init_guru)(p, (int)d, (int *)N, M, n,
    ini->guru == 2 ? WINDOW_HELP_ESTIMATE_m / 2 : WINDOW_HELP_ESTIMATE_m,
    BASE_FLAGS | ini->flags, FFTW_FLAGS);
}

static void print_case(const char *what, const INT d, const int *N,
  const int M, const char *name)
{
  INT t;

  printf("%-14s d = %d, N = [", what, (int)d);
  for (t = 0; t < d; t++)
    printf("%s%d", t ? ", " : "", N[t]);
  printf("], M = %d, %-36s", M, name);
}

static void fill_coefficients(X(plan) *p)
{
  INT i;

  for (i = 0; i < p->N_total; i++)
    p->f_hat[i] = (Y(drand48)() - K(0.5)) + (Y(drand48)() - K(0.5)) * II;
}

static void fill_samples(X(plan) *p)
{
  INT j;

  for (j = 0; j < p->M_total; j++)
    p->f[j] = Y(drand48)() - K(0.5);
}

/* Fast against direct. Trafo error is sup |f - f_direct| over
 * sum c_k |f_hat_k|, adjoint error sup |f_hat - f_hat_direct| / c_k over
 * sum |f_j|; both are the errors of the complex NFFT of the Hermitian
 * extension. */
static int check_single(const INT d, const int *N, const int M,
  const init_t *ini, const int adjoint)
{
  X(plan) p;
  const char *reason;
  R err = K(0.0), norm = K(0.0), bound;
  INT i, t;
  int ok;

  print_case(adjoint ? "rnfft_adjoint" : "rnfft", d, N, M, ini->name);
  init_plan(&p, ini, d, N, M);
  for (i = 0; i < d * M; i++)
    p.x[i] = Y(drand48)() - K(0.5);

  reason = X(check)(&p);
  for (t = 0; !reason && t < d; t++)
    if (p.m >= p.N[t])
      reason = "Cut-off m not below N";
  if (reason)
  {
    printf(" -> OK   (%s)\n", reason);
    X(finalize)(&p);
    return 1;
  }

  if (p.flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&p);

  if (!adjoint)
  {
    R *ref = Y(malloc)((size_t)M * sizeof(R));
    C *in = Y(malloc)((size_t)p.N_total * sizeof(C));

    fill_coefficients(&p);
    for (i = 0; i < p.N_total; i++)
      in[i] = p.f_hat[i];
    X(trafo_direct)(&p);
    for (i = 0; i < M; i++)
      ref[i] = p.f[i];
    X(trafo)(&p);
    for (i = 0; i < M; i++)
      err = FMAX(err, FABS(p.f[i] - ref[i]));
    for (i = 0; i < p.N_total; i++)
      norm += weight(&p, i) * CABS(in[i]);
    Y(free)(in);
    Y(free)(ref);
  }
  else
  {
    C *ref = Y(malloc)((size_t)p.N_total * sizeof(C));

    fill_samples(&p);
    X(adjoint_direct)(&p);
    for (i = 0; i < p.N_total; i++)
      ref[i] = p.f_hat[i];
    X(adjoint)(&p);
    for (i = 0; i < p.N_total; i++)
      err = FMAX(err, CABS(p.f_hat[i] - ref[i]) / weight(&p, i));
    for (i = 0; i < M; i++)
      norm += FABS(p.f[i]);
    Y(free)(ref);
  }

  err /= norm;
  bound = err_trafo(&p) + DIRECT_ERR;
  ok = err < bound;
  printf(", m = %2d -> %-4s " __FE__ " (" __FE__ ")\n", (int)p.m,
    ok ? "OK" : "FAIL", err, bound);
  {
    int Nl[4];

    for (t = 0; t < d; t++)
      Nl[t] = N[t];
    accuracy_log_append("rnfft", "online", (int)d, Nl, M, ini->name,
      adjoint ? "adjoint" : "trafo", (long double)err, (long double)bound, ok);
  }

  X(finalize)(&p);
  return ok;
}

static void check_sizes(const INT d, const int *N, const int *M, const int nc,
  const int adjoint)
{
  int ok = 1, c;
  size_t i;

  Y(srand48)(SEED);
  for (c = 0; c < nc; c++)
    for (i = 0; i < sizeof(inits) / sizeof(inits[0]); i++)
      ok &= check_single(d, N + c * d, M[c], &inits[i], adjoint);
  CU_ASSERT(ok);
}

/* d = 1, 2, 3 with sizes that hit both even and odd row counts, small and
 * large N relative to m, and in-place padding. */
static const int N1[] = {16, 50, 256}, M1[] = {50, 80, 300};
static const int N2[] = {16, 16, 12, 30, 32, 24}, M2[] = {100, 200, 150};
static const int N3[] = {12, 12, 12, 16, 10, 14}, M3[] = {100, 150};
static const int N4[] = {16, 16, 16, 16, 14, 18, 16, 20}, M4[] = {60, 60};

static void check_online_all(const int adjoint)
{
  check_sizes(1, N1, M1, 3, adjoint);
  check_sizes(2, N2, M2, 3, adjoint);
  check_sizes(3, N3, M3, 2, adjoint);
}

void X(check_online)(void)
{
  check_online_all(0);
}

void X(check_adjoint_online)(void)
{
  check_online_all(1);
}

void X(check_online_4d)(void)
{
  check_sizes(4, N4, M4, 2, 0);
}

void X(check_adjoint_online_4d)(void)
{
  check_sizes(4, N4, M4, 2, 1);
}

/* Index in the full spectrum of the complex NFFT of half spectrum index i,
 * negated in every dimension if neg; -1 if outside the full spectrum. */
static INT full_index(const X(plan) *p, INT i, const int neg)
{
  const INT d = p->d, H = p->N[d - 1] / 2 + 1;
  INT k[4], t, r = 0;

  k[d - 1] = i % H;
  i /= H;
  for (t = d - 2; t >= 0; t--)
  {
    k[t] = i % p->N[t] - p->N[t] / 2;
    i /= p->N[t];
  }
  for (t = 0; t < d; t++)
  {
    const INT kt = neg ? -k[t] : k[t];

    if (kt < -p->N[t] / 2 || kt >= p->N[t] / 2)
      return -1;
    r = r * p->N[t] + kt + p->N[t] / 2;
  }
  return r;
}

static void hermitian_plans(NFFT(plan) *c, X(plan) *p, const INT d,
  const int *N, const int M)
{
  INT j;

  NFFT(init)(c, (int)d, (int *)N, M);
  X(init)(p, (int)d, (int *)N, M);
  for (j = 0; j < d * M; j++)
    c->x[j] = p->x[j] = Y(drand48)() - K(0.5);
}

/* A random Hermitian spectrum f_hat_{-k} = conj(f_hat_k) of the complex NFFT,
 * zero where -k leaves the spectrum, is given to the RNFFT by its half. Both
 * direct trafos agree. */
static void check_hermitian_single(const INT d, const int *N, const int M,
  int *ok)
{
  NFFT(plan) c;
  X(plan) p;
  R err = K(0.0), im = K(0.0), norm = K(0.0);
  INT i, j;

  print_case("rnfft_herm", d, N, M, "init");
  hermitian_plans(&c, &p, d, N, M);

  for (i = 0; i < c.N_total; i++)
    c.f_hat[i] = K(0.0);
  for (i = 0; i < p.N_total; i++)
  {
    const INT a = full_index(&p, i, 0), b = full_index(&p, i, 1);
    C v = (Y(drand48)() - K(0.5)) + (Y(drand48)() - K(0.5)) * II;

    if (a < 0 || b < 0)
    {
      p.f_hat[i] = K(0.0);
      continue;
    }
    if (a == b)
      v = CREAL(v);
    c.f_hat[a] = v;
    c.f_hat[b] = CONJ(v);
  }
  for (i = 0; i < p.N_total; i++)
  {
    const INT a = full_index(&p, i, 0), b = full_index(&p, i, 1);

    p.f_hat[i] = (a < 0 || b < 0) ? K(0.0) : c.f_hat[a];
  }

  NFFT(trafo_direct)(&c);
  X(trafo_direct)(&p);
  for (j = 0; j < M; j++)
  {
    err = FMAX(err, FABS(CREAL(c.f[j]) - p.f[j]));
    im = FMAX(im, FABS(CIMAG(c.f[j])));
  }
  for (i = 0; i < c.N_total; i++)
    norm += CABS(c.f_hat[i]);
  err = FMAX(err, im) / norm;

  *ok &= err < DIRECT_ERR;
  printf(" -> %-4s " __FE__ " (" __FE__ ")\n", err < DIRECT_ERR ? "OK" : "FAIL",
    err, DIRECT_ERR);
  NFFT(finalize)(&c);
  X(finalize)(&p);
}

/* For real samples the RNFFT adjoint is c_k times the complex adjoint on the
 * half spectrum. */
static void check_adjoint_hermitian_single(const INT d, const int *N,
  const int M, int *ok)
{
  NFFT(plan) c;
  X(plan) p;
  R err = K(0.0), norm = K(0.0);
  INT i, j;

  print_case("rnfft_adj_herm", d, N, M, "init");
  hermitian_plans(&c, &p, d, N, M);

  for (j = 0; j < M; j++)
  {
    p.f[j] = Y(drand48)() - K(0.5);
    c.f[j] = p.f[j];
    norm += FABS(p.f[j]);
  }
  NFFT(adjoint_direct)(&c);
  X(adjoint_direct)(&p);
  for (i = 0; i < p.N_total; i++)
  {
    const INT a = full_index(&p, i, 0);

    if (a >= 0)
      err = FMAX(err, CABS(weight(&p, i) * c.f_hat[a] - p.f_hat[i])
        / weight(&p, i));
  }
  err /= norm;

  *ok &= err < DIRECT_ERR;
  printf(" -> %-4s " __FE__ " (" __FE__ ")\n", err < DIRECT_ERR ? "OK" : "FAIL",
    err, DIRECT_ERR);
  NFFT(finalize)(&c);
  X(finalize)(&p);
}

static const int NH[] = {64, 16, 24, 12, 8, 16}, MH[] = {100, 120, 80};

void X(check_direct_hermitian)(void)
{
  int ok = 1;

  Y(srand48)(SEED);
  check_hermitian_single(1, NH, MH[0], &ok);
  check_hermitian_single(2, NH + 1, MH[1], &ok);
  check_hermitian_single(3, NH + 3, MH[2], &ok);
  CU_ASSERT(ok);
}

void X(check_adjoint_direct_hermitian)(void)
{
  int ok = 1;

  Y(srand48)(SEED);
  check_adjoint_hermitian_single(1, NH, MH[0], &ok);
  check_adjoint_hermitian_single(2, NH + 1, MH[1], &ok);
  check_adjoint_hermitian_single(3, NH + 3, MH[2], &ok);
  CU_ASSERT(ok);
}

/* sum_j (A f_hat)_j y_j = Re sum_k f_hat_k conj((A^T y)_k), direct and fast,
 * relative to |A f_hat| |y|. */
static int check_inner_product_single(const INT d, const int *N, const int M,
  const init_t *ini)
{
  X(plan) p;
  C *fh;
  R *y, *Af, lhs[2], rhs[2], scale[2], bound = K(1000.0)
    * Y(float_property)(NFFT_EPSILON);
  INT i;
  int fast, ok = 1;

  print_case("rnfft_ip", d, N, M, ini->name);
  init_plan(&p, ini, d, N, M);
  for (i = 0; i < d * M; i++)
    p.x[i] = Y(drand48)() - K(0.5);
  if (p.flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&p);

  fh = Y(malloc)((size_t)p.N_total * sizeof(C));
  y = Y(malloc)((size_t)M * sizeof(R));
  Af = Y(malloc)((size_t)M * sizeof(R));
  fill_coefficients(&p);
  for (i = 0; i < p.N_total; i++)
    fh[i] = p.f_hat[i];
  for (i = 0; i < M; i++)
    y[i] = Y(drand48)() - K(0.5);

  for (fast = 0; fast < 2; fast++)
  {
    R na = K(0.0), ny = K(0.0);

    for (i = 0; i < p.N_total; i++)
      p.f_hat[i] = fh[i];
    if (fast)
      X(trafo)(&p);
    else
      X(trafo_direct)(&p);
    for (i = 0; i < M; i++)
      Af[i] = p.f[i];
    for (i = 0; i < M; i++)
      p.f[i] = y[i];
    if (fast)
      X(adjoint)(&p);
    else
      X(adjoint_direct)(&p);

    lhs[fast] = rhs[fast] = K(0.0);
    for (i = 0; i < M; i++)
    {
      lhs[fast] += Af[i] * y[i];
      na += Af[i] * Af[i];
      ny += y[i] * y[i];
    }
    for (i = 0; i < p.N_total; i++)
      rhs[fast] += CREAL(fh[i] * CONJ(p.f_hat[i]));
    scale[fast] = SQRT(na * ny);
    ok &= FABS(lhs[fast] - rhs[fast]) < bound * scale[fast];
  }

  printf(" -> %-4s " __FE__ ", " __FE__ " (" __FE__ ")\n", ok ? "OK" : "FAIL",
    FABS(lhs[0] - rhs[0]) / scale[0], FABS(lhs[1] - rhs[1]) / scale[1], bound);
  Y(free)(Af);
  Y(free)(y);
  Y(free)(fh);
  X(finalize)(&p);
  return ok;
}

void X(check_inner_product)(void)
{
  int ok = 1;
  size_t i;

  Y(srand48)(SEED);
  for (i = 0; i < sizeof(inits) / sizeof(inits[0]); i++)
  {
    ok &= check_inner_product_single(1, N1 + 1, M1[1], &inits[i]);
    ok &= check_inner_product_single(2, N2 + 2, M2[1], &inits[i]);
    ok &= check_inner_product_single(3, N3, M3[0], &inits[i]);
  }
  CU_ASSERT(ok);
}

/* Minimum-norm solution for data in the range of the trafo: in the plane
 * k_{d-1} = 0 the trafo sees only f_hat_k + conj(f_hat_{-k}), so the
 * projection keeps half of it at k and its conjugate at -k. */
static void project(const X(plan) *p, C *fh)
{
  const INT H = p->N[p->d - 1] / 2 + 1;
  INT i;

  for (i = 0; i < p->N_total; i += H)
  {
    const INT a = full_index(p, i, 0), b = full_index(p, i, 1);
    INT ib;

    if (b < 0)
      continue;
    if (a == b)
    {
      fh[i] = CREAL(fh[i]);
      continue;
    }
    ib = i + (b - a) / p->N[p->d - 1] * H;
    if (ib > i)
    {
      const C s = (fh[i] + CONJ(fh[ib])) * K(0.5);

      fh[i] = s;
      fh[ib] = CONJ(s);
    }
  }
}

/* From exact samples, CGNR on an overdetermined system recovers the
 * minimum-norm coefficients, and CGNE on an underdetermined one fits the
 * samples. */
static int check_solver_single(const INT d, const int *N, const int M,
  const unsigned flags, const int iter)
{
  X(plan) p;
  SOLVER(plan_c2r) ip;
  C *fh;
  R err = K(0.0), norm = K(0.0), bound = K(1.0E4)
    * Y(float_property)(NFFT_EPSILON);
  INT i;
  int l, ok;

  print_case(flags & CGNR ? "solver_c2r_cgnr" : "solver_c2r_cgne", d, N, M,
    "init");
  X(init)(&p, (int)d, (int *)N, M);
  for (i = 0; i < d * M; i++)
    p.x[i] = Y(drand48)() - K(0.5);
  X(precompute_one_psi)(&p);

  fh = Y(malloc)((size_t)p.N_total * sizeof(C));
  fill_coefficients(&p);
  project(&p, p.f_hat);
  for (i = 0; i < p.N_total; i++)
    fh[i] = p.f_hat[i];
  X(trafo)(&p);

  SOLVER(init_advanced_c2r)(&ip, (NFFT(mv_plan_c2r) *)&p, flags);
  for (i = 0; i < M; i++)
    ip.y[i] = p.f[i];
  for (i = 0; i < p.N_total; i++)
    ip.f_hat_iter[i] = K(0.0);
  SOLVER(before_loop_c2r)(&ip);
  for (l = 0; l < iter; l++)
    SOLVER(loop_one_step_c2r)(&ip);

  if (flags & CGNR)
  {
    for (i = 0; i < p.N_total; i++)
    {
      err = FMAX(err, CABS(ip.f_hat_iter[i] - fh[i]));
      norm = FMAX(norm, CABS(fh[i]));
    }
    err /= norm;
  }
  else
  {
    for (i = 0; i < M; i++)
      norm += ip.y[i] * ip.y[i];
    err = SQRT(ip.dot_r_iter / norm);
  }
  ok = err < bound;
  printf(" -> %-4s " __FE__ " (" __FE__ ")\n", ok ? "OK" : "FAIL", err, bound);

  SOLVER(finalize_c2r)(&ip);
  Y(free)(fh);
  X(finalize)(&p);
  return ok;
}

void X(check_solver_c2r)(void)
{
  static const int Ns1[] = {32}, Ns2[] = {12, 16};
  int ok = 1;

  Y(srand48)(SEED);
  ok &= check_solver_single(1, Ns1, 200, CGNR, 60);
  ok &= check_solver_single(2, Ns2, 800, CGNR, 100);
  ok &= check_solver_single(1, Ns1, 12, CGNE, 30);
  ok &= check_solver_single(2, Ns2, 60, CGNE, 60);
  CU_ASSERT(ok);
}
