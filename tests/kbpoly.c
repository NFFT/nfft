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

#include <stdio.h>
#include <string.h>
#include <CUnit/CUnit.h>
#include "nfft3.h"
#include "infft.h"
#include "kbpoly.h"

/* The window macros in infft.h fall through to Kaiser-Bessel. */
#if !defined(DIRAC_DELTA) && !defined(GAUSSIAN) && !defined(B_SPLINE) \
    && !defined(SINC_POWER)
#define WINDOW_IS_KAISER_BESSEL 1
#endif

#if defined(WINDOW_IS_KAISER_BESSEL)

/* Errors are absolute against the window peak, which is what the psi table
 * owes the transform. Relative error is meaningless in a tail 1E-19 of the
 * peak. */
static R kb_peak(const R b, const R lg_tail, const R peak_inv, const R mr)
{
  return Y(kb_phi)(b, lg_tail, peak_inv, mr, K(0.0));
}

/* Expected fit error at the degree Y(kb_poly_degree) picks: the tolerance, or
 * 7.5 eps (rounding plateau up to 6 eps, 1 eps truncation, margin). */
/* 5% of the transform error c / I0(a R), c = 1.10, a = 2 pi sqrt(1 - 1/sigma),
 * or 2 eps where that is larger: the table rule of gen/kbpolygen. */
static R pick_target(const R b, const R reach)
{
  const R sigma = KPI / (K2PI - b);
  const R ar = K2PI * SQRT(K(1.0) - K(1.0) / sigma) * reach;
  const R alias = K(1.10) * EXP(-ar) / Y(bessel_i0_exp_scaled)(ar);
  const R tol = K(0.05) * alias;
  const R lowest = K(2.0) * Y(float_property)(NFFT_EPSILON);

  return IF(tol > lowest, tol, lowest);
}

#define KB_FN_NXI 40
#define KB_FN_NT 32

/* The transform's error functional of a table, on a grid of frequencies xi
 * in [0, 1/(2 sigma)] and offsets t in [0, 1): J_inf = max |D / C - 1| with D
 * the sum of the 2m+2 taps times exp(-2 pi i xi (position)), C = PHI_HUT; and
 * the polynomial's own part max |D_poly - D_phi| / C, D_phi the same sum of
 * the closed-form window. Tap l of a run anchored at the grid point below the
 * node sits at nx0 - l, nx0 = m + t. */
static void functional(const R *coef, const INT m, const INT deg, const R b,
    const R lt, const R pki, const R i0e_pki, R *j_inf, R *poly_part)
{
  const INT w = 2 * m + 2;
  const R reach = (R)m + K(1.0), sigma = KPI / (K2PI - b);
  R *got = (R*) Y(malloc)((size_t)(2 * w) * sizeof(R)), *ref = got + w;
  INT i, q, l;

  *j_inf = *poly_part = K(0.0);

  for (i = 0; i < KB_FN_NXI; i++)
  {
    const R xi = (R)i / (R)(KB_FN_NXI - 1) / (K(2.0) * sigma);
    const R C = Y(kb_phi_hut)(b, i0e_pki, reach, K(1.0), xi);

    for (q = 0; q < KB_FN_NT; q++)
    {
      const R t = ((R)q + K(0.5)) / (R)KB_FN_NT, nx0 = (R)m + t;
      R dre = K(0.0), dim = K(0.0), ere = K(0.0), eim = K(0.0), e;

      Y(kb_poly_run)(got, coef, m, deg, nx0);
      Y(kb_phi_run)(ref, b, lt, pki, reach, m, nx0);

      for (l = 0; l < w; l++)
      {
        const R ph = -K2PI * xi * (nx0 - (R)l), c = COS(ph), sn = SIN(ph);

        dre += got[l] * c;
        dim += got[l] * sn;
        ere += (got[l] - ref[l]) * c;
        eim += (got[l] - ref[l]) * sn;
      }

      e = SQRT((dre / C - K(1.0)) * (dre / C - K(1.0)) + dim * dim / (C * C));
      *j_inf = IF(e > *j_inf, e, *j_inf);
      e = SQRT(ere * ere + eim * eim) / C;
      *poly_part = IF(e > *poly_part, e, *poly_part);
    }
  }

  Y(free)(got);
}

#define KB_N 120
static const INT kb_n[] = {132, 150, 180, 210, 240, 300, 360, 480};
#define KB_NSIGMA (sizeof(kb_n) / sizeof(kb_n[0]))

/* Up to where the transform reaches the precision floor at sigma 1.1. */
#if MANT_DIG >= 113
  #define KB_M_MAX 40
#elif MANT_DIG >= 64
  #define KB_M_MAX 24
#elif MANT_DIG >= 53
  #define KB_M_MAX 20
#else
  #define KB_M_MAX 8
#endif

static INT kb_m_next(const INT m)
{
  return IF(m < 12, m + 1, m + 4);
}

void X(check_kaiser_bessel_poly)(void)
{
  unsigned int s;
  INT m;

  printf("KAISER-BESSEL POLYNOMIAL\n------------------------\n");

  for (s = 0; s < KB_NSIGMA; s++)
    for (m = 2; m <= KB_M_MAX; m = kb_m_next(m))
    {
      const R sigma = (R)kb_n[s] / (R)KB_N;
      const R b = KPI * (K(2.0) - K(1.0) / sigma);
      const R reach = (R)m + K(1.0);
      const R lt = Y(bessel_i0_logtail)(reach * b);
      const R pki = EXP(-reach * b - lt);
      const R i0e_pki = K(1.0) / Y(bessel_i0_exp_scaled)(reach * b);
      const R peak = kb_peak(b, lt, pki, reach);
      const INT deg = Y(kb_poly_degree)(b, reach, m);
      const INT cols = KB_POLY_COLS(m);
      const R tol = pick_target(b, reach);
      /* The functional is formed in R from taps that each carry a few ulp,
       * divided by C, which is smallest at the band edge. */
      const R cmin = Y(kb_phi_hut)(b, i0e_pki, reach, K(2.0) * sigma, K(1.0));
      const R slack = K(4.0) * (R)(2 * m + 2) * Y(float_property)(NFFT_EPSILON)
          / cmin;
      const R bound = K(1.1) * tol + slack;
      R *coef = (R*) Y(malloc)((size_t)(KB_POLY_ROWS(deg) * cols) * sizeof(R));
      R j_inf, part, worst = K(0.0);
      INT i;
      int ok;

      Y(kb_poly_fit)(coef, b, lt, pki, reach, m, deg);
      functional(coef, m, deg, b, lt, pki, i0e_pki, &j_inf, &part);

      /* The single-point evaluator, out to |nx| = m + 2 where PRE_LIN_PSI
       * reads it, against the closed form. Its error is the interpolation
       * error itself, up to two orders above its share of the transform
       * error, which cancels across the taps; PRE_LIN_PSI's own
       * interpolation error is 1e-4 and above, so this is a sanity bound,
       * with 1e-12 for the formats whose tolerance is below any use of a
       * single point. */
      for (i = 0; i <= 512; i++)
      {
        const R nx = (R)(m + 2) * ((R)i / K(256.0) - K(1.0));
        const R err = ABS(Y(kb_poly_phi)(coef, m, deg, nx)
            - Y(kb_phi)(b, lt, pki, reach, nx)) / peak;

        worst = IF(err > worst, err, worst);
      }

      ok = IF(part < bound && worst < K(100.0) * bound + K(1E-12)
          && deg <= KB_POLY_DEG_MAX, 1, 0);
      printf("poly[sigma=" __FE__ ", m=%2td, deg=%2td] part = " __FE__
          " phi = " __FE__ " %-2s " __FE__ ", J_inf = " __FE__ " -> %-4s\n",
          sigma, (ptrdiff_t)m, (ptrdiff_t)deg, part, worst,
          IF(ok == 0, ">=", "<"), bound, j_inf, IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);

      Y(free)(coef);
    }

  printf("\n");
}

/* The least-squares correction lowers the functional's J_inf, or leaves the
 * table alone where its guard says the error is at the rounding floor. */
void X(check_kaiser_bessel_poly_opt)(void)
{
  unsigned int s;
  INT m;

  printf("KAISER-BESSEL POLYNOMIAL OPT\n----------------------------\n");

  for (s = 0; s < KB_NSIGMA; s++)
    for (m = 2; m <= KB_M_MAX; m = kb_m_next(m))
    {
      const R sigma = (R)kb_n[s] / (R)KB_N;
      const R b = KPI * (K(2.0) - K(1.0) / sigma);
      const R reach = (R)m + K(1.0);
      const R lt = Y(bessel_i0_logtail)(reach * b);
      const R pki = EXP(-reach * b - lt);
      const R i0e_pki = K(1.0) / Y(bessel_i0_exp_scaled)(reach * b);
      const INT deg = Y(kb_poly_degree)(b, reach, m);
      const INT n = KB_POLY_ROWS(deg) * KB_POLY_COLS(m);
      const int admitted = Y(kb_poly_correct_ok)(b, m);
      const R cmin = Y(kb_phi_hut)(b, i0e_pki, reach, K(2.0) * sigma, K(1.0));
      const R slack = K(4.0) * (R)(2 * m + 2) * Y(float_property)(NFFT_EPSILON)
          / cmin;
      R *coef = (R*) Y(malloc)((size_t)(2 * n) * sizeof(R)), *plain = coef + n;
      R j0, j1, p0, p1;
      INT i;
      int ok;

      Y(kb_poly_fit)(coef, b, lt, pki, reach, m, deg);
      memcpy(plain, coef, (size_t)n * sizeof(R));
      functional(plain, m, deg, b, lt, pki, i0e_pki, &j0, &p0);

      if (admitted)
        Y(kb_poly_correct)(coef, b, i0e_pki, m, deg);

      functional(coef, m, deg, b, lt, pki, i0e_pki, &j1, &p1);

      if (admitted)
        ok = IF(j1 <= K(1.02) * j0 + slack, 1, 0);
      else
      {
        ok = 1;

        for (i = 0; i < n; i++)
          ok = IF(coef[i] == plain[i], ok, 0);
      }

      printf("opt[sigma=" __FE__ ", m=%2td, deg=%2td] %s J_inf " __FE__
          " -> " __FE__ " (%.2fx) -> %-4s\n", sigma, (ptrdiff_t)m,
          (ptrdiff_t)deg, admitted ? "corrected" : "guarded  ", j0, j1,
          (double)(j0 / j1), IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);

      Y(free)(coef);
    }

  printf("\n");
}

/* Worst |f_j - exp(-2 pi i k x_j)| over the nodes and the band-edge
 * frequencies, f_hat = e_k. The nodes are multiples of 2^-16, so k x_j is
 * exact and the reference phase reduces exactly. */
static R edge_err(NFFT(plan) *p)
{
  const INT N = p->N[0], ks[] = {-N / 2, -N / 2 + 1, -N / 2 + 2, N / 2 - 1};
  C *f_hat = (C*)p->f_hat, *f = (C*)p->f;
  R worst = K(0.0);
  unsigned int i;
  INT j;

  for (i = 0; i < sizeof(ks) / sizeof(ks[0]); i++)
  {
    for (j = 0; j < N; j++)
      f_hat[j] = K(0.0);

    f_hat[ks[i] + N / 2] = K(1.0);
    NFFT(trafo)(p);

    for (j = 0; j < p->M_total; j++)
    {
      const R ph = (R)ks[i] * p->x[j];
      const R e = CABS(f[j] - CEXP(-K2PI * I * (ph - ROUND(ph))));

      worst = IF(e > worst, e, worst);
    }
  }

  return worst;
}

/* X(check) rejects sigma off [1.1, 4], where the degree table ends. */
void X(check_kaiser_bessel_poly_sigma)(void)
{
  const unsigned flags = PRE_PHI_HUT | MALLOC_X | MALLOC_F | MALLOC_F_HAT
      | FFTW_INIT | FFT_OUT_OF_PLACE;
  const int ns[] = {124, 132, 480, 540};
  unsigned int i;

  printf("KAISER-BESSEL POLYNOMIAL SIGMA\n-----------------------------\n");

  for (i = 0; i < sizeof(ns) / sizeof(ns[0]); i++)
  {
    int N[1] = {KB_N}, n[1] = {ns[i]};
    const int in = ns[i] >= 132 && ns[i] <= 480;
    NFFT(plan) p;
    const char *msg;
    int ok;

    NFFT(init_guru)(&p, 1, N, 8, n, 4, flags,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    p.x[0] = p.x[1] = p.x[2] = p.x[3] = K(0.0);
    p.x[4] = p.x[5] = p.x[6] = p.x[7] = K(0.25);
    msg = NFFT(check)(&p);
    ok = IF(in, msg == 0, msg != 0);
    printf("check[sigma=%5.3f] %s -> %s\n", (double)ns[i] / KB_N,
        msg ? msg : "accepted", ok ? "OK" : "FAIL");
    CU_ASSERT(ok);
    NFFT(finalize)(&p);
  }

  printf("\n");
}

/* Overwrites a 1-d PRE_PSI table with the closed-form window, the reference
 * the polynomial is held to. n is the grid the module anchors on, u the
 * anchor of each run as in the module's uo(): floor(n x) - m. */
static void analytic_psi(R *psi, const R *x, const R *b, const INT M,
    const INT n, const INT m)
{
  INT j;

  for (j = 0; j < M; j++)
  {
    const R xn = x[j] * (R)n;
    const INT u = (INT)FLOOR(xn) - m;

    Y(kb_phi_run)(psi + j * (2 * m + 2), b[0], b[1], b[3], (R)m + K(1.0), m,
        xn - (R)u);
  }
}

#if defined(HAVE_NFCT) || defined(HAVE_NFST)
/* The same for the cosine and sine series, f_hat = e_k at the three highest
 * frequencies, nodes in [0, 1/2). Index l of f_hat is frequency l + OFF. */
#define TRIG_EDGE_ERR(MOD, TRIG, OFF) \
static R MOD##_edge_err(MOD(plan) *p) \
{ \
  const INT N = p->N[0]; \
  R worst = K(0.0); \
  INT k, j; \
\
  for (k = N - 3; k < N; k++) \
  { \
    for (j = 0; j < N - (OFF); j++) \
      p->f_hat[j] = K(0.0); \
\
    p->f_hat[k - (OFF)] = K(1.0); \
    MOD(trafo)(p); \
\
    for (j = 0; j < p->M_total; j++) \
    { \
      const R ph = (R)k * p->x[j]; \
      const R e = ABS(p->f[j] - TRIG(K2PI * (ph - ROUND(ph)))); \
\
      worst = IF(e > worst, e, worst); \
    } \
  } \
\
  return worst; \
}

/* n = sigma N + NADJ, OFF the lowest frequency, NANCH the anchor grid in n. */
#define TRIG_POLY_TRAFO(MOD, NAME, OFF, NADJ, NANCH) \
static void MOD##_poly_trafo(void) \
{ \
  const INT M = 2000; \
  const unsigned flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F \
      | MALLOC_F_HAT | FFTW_INIT | FFT_OUT_OF_PLACE; \
  unsigned int s; \
  INT m, j; \
  int opt; \
\
  for (s = 0; s < KB_NSIGMA; s++) \
    for (m = 2; m <= KB_M_MAX; m = kb_m_next(m)) \
    for (opt = 0; opt < 2; opt++) \
    { \
      int N[1] = {KB_N}, n[1] = {(int)kb_n[s] + (NADJ)}; \
      MOD(plan) pp, pa; \
      R ep, ea, uk, bound; \
      int ok; \
\
      MOD(init_guru)(&pp, 1, N, (int)M, n, (int)m, \
          flags | (opt ? OPT_WINDOW : 0), FFTW_ESTIMATE | FFTW_DESTROY_INPUT); \
      MOD(init_guru)(&pa, 1, N, (int)M, n, (int)m, flags, \
          FFTW_ESTIMATE | FFTW_DESTROY_INPUT); \
\
      for (j = 0; j < M; j++) \
        pp.x[j] = pa.x[j] = (R)(INT)(Y(drand48)() * K(32768.0)) \
            / K(65536.0); \
\
      MOD(precompute_one_psi)(&pp); \
      analytic_psi(pa.psi, pa.x, pa.b, M, NANCH(n[0]), m); \
\
      ep = MOD##_edge_err(&pp); \
      ea = MOD##_edge_err(&pa); \
      uk = Y(float_property)(NFFT_EPSILON) \
          * pa.c_phi_inv[0][KB_N - 1 - (OFF)] / pa.c_phi_inv[0][0]; \
      bound = K(1.05) * ea + K(2.0) * uk; \
      ok = IF(ep <= bound, 1, 0); \
\
      printf(NAME "[sigma=%5.2f, m=%2td, deg=%2td%s] poly %10.3e analytic " \
          "%10.3e ratio %6.3f, in u*kappa %7.2f %7.2f -> %-4s\n", \
          (double)((R)kb_n[s] / (R)KB_N), (ptrdiff_t)m, \
          (ptrdiff_t)(INT)pp.spline_coeffs[0], opt ? ", opt" : "", \
          (double)ep, (double)ea, \
          (double)(ep / ea), (double)(ep / uk), (double)(ea / uk), \
          IF(ok == 0, "FAIL", "OK")); \
      CU_ASSERT(ok == 1); \
\
      MOD(finalize)(&pa); \
      MOD(finalize)(&pp); \
    } \
}
#endif

#ifdef HAVE_NFCT
TRIG_EDGE_ERR(NFCT, COS, 0)
#define NFCT_ANCH(n) (2 * ((n) - 1))
TRIG_POLY_TRAFO(NFCT, "nfct", 0, 1, NFCT_ANCH)
#endif

#ifdef HAVE_NFST
TRIG_EDGE_ERR(NFST, SIN, 1)
#define NFST_ANCH(n) (2 * ((n) + 1))
TRIG_POLY_TRAFO(NFST, "nfst", 1, -1, NFST_ANCH)
#endif

void X(check_kaiser_bessel_poly_trafo)(void)
{
  const INT M = 2000;
  const unsigned flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F
      | MALLOC_F_HAT | FFTW_INIT | FFT_OUT_OF_PLACE;
  unsigned int s;
  INT m, j;
  int opt;

  printf("KAISER-BESSEL POLYNOMIAL TRAFO\n------------------------------\n");
  Y(srand48)(1234567890L);

  for (s = 0; s < KB_NSIGMA; s++)
    for (m = 2; m <= KB_M_MAX; m = kb_m_next(m))
    for (opt = 0; opt < 2; opt++)
    {
      int N[1] = {KB_N}, n[1] = {(int)kb_n[s]};
      NFFT(plan) pp, pa;
      R ep, ea;

      NFFT(init_guru)(&pp, 1, N, (int)M, n, (int)m,
          flags | (opt ? OPT_WINDOW : 0), FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      NFFT(init_guru)(&pa, 1, N, (int)M, n, (int)m, flags,
          FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

      for (j = 0; j < M; j++)
        pp.x[j] = pa.x[j] = (R)(INT)(Y(drand48)() * K(65536.0)) / K(65536.0)
            - K(0.5);

      NFFT(precompute_one_psi)(&pp);
      analytic_psi(pa.psi, pa.x, pa.b, M, n[0], m);

      ep = edge_err(&pp);
      ea = edge_err(&pa);

      {
        /* The polynomial may add 5% where truncation dominates, and two u
         * kappa of float Horner round-off, kappa the spread of c_phi_inv. The
         * round-off varies with FMA contraction, so with the build CPU. */
        const R uk = Y(float_property)(NFFT_EPSILON) * pa.c_phi_inv[0][0]
            / pa.c_phi_inv[0][KB_N / 2];
        const R bound = K(1.05) * ea + K(2.0) * uk;
        const int ok = IF(ep <= bound, 1, 0);

        printf("trafo[sigma=%5.2f, m=%2td, deg=%2td%s] poly %10.3e analytic "
            "%10.3e ratio %6.3f, in u*kappa %7.2f %7.2f -> %-4s\n",
            (double)((R)kb_n[s] / (R)KB_N), (ptrdiff_t)m,
            (ptrdiff_t)(INT)pp.spline_coeffs[0], opt ? ", opt" : "",
            (double)ep, (double)ea, (double)(ep / ea), (double)(ep / uk),
            (double)(ea / uk), IF(ok == 0, "FAIL", "OK"));
        CU_ASSERT(ok == 1);
      }

      NFFT(finalize)(&pa);
      NFFT(finalize)(&pp);
    }

#ifdef HAVE_NFCT
  NFCT_poly_trafo();
#endif
#ifdef HAVE_NFST
  NFST_poly_trafo();
#endif

  printf("\n");
}

#else

void X(check_kaiser_bessel_poly)(void)
{
  printf("KAISER-BESSEL POLYNOMIAL\n------------------------\nskipped: "
      "window is not Kaiser-Bessel\n\n");
  CU_ASSERT(1);
}

void X(check_kaiser_bessel_poly_trafo)(void)
{
  printf("KAISER-BESSEL POLYNOMIAL TRAFO\n------------------------------\n"
      "skipped: window is not Kaiser-Bessel\n\n");
  CU_ASSERT(1);
}

void X(check_kaiser_bessel_poly_opt)(void)
{
  printf("KAISER-BESSEL POLYNOMIAL OPT\n----------------------------\n");
  printf("skipped: not the Kaiser-Bessel window\n\n");
  CU_ASSERT(1);
}

void X(check_kaiser_bessel_poly_sigma)(void)
{
  printf("KAISER-BESSEL POLYNOMIAL SIGMA\n-----------------------------\n"
      "skipped: window is not Kaiser-Bessel\n\n");
  CU_ASSERT(1);
}

#endif
