/*
 * Copyright (c) 2026 The NFFT Authors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
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
static R pick_target(const R b, const R reach)
{
  const R tol = K(0.05) * KPI * EXP(-b * reach);
  const R lowest = K(7.5) * Y(float_property)(NFFT_EPSILON);

  return IF(tol > lowest, tol, lowest);
}

/* sigma = n / N with N = 120 and n even. */
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
      const R peak = kb_peak(b, lt, pki, reach);
      const INT deg = Y(kb_poly_degree)(b, reach, m), w = 2 * m + 2;
      const INT cols = KB_POLY_COLS(m);
      const R bound = K(2.0) * pick_target(b, reach);
      R *coef = (R*) Y(malloc)((size_t)(KB_POLY_ROWS(deg) * cols) * sizeof(R));
      R *got = (R*) Y(malloc)((size_t)w * sizeof(R));
      R worst = K(0.0);
      INT i, l;
      int ok;

      Y(kb_poly_fit)(coef, b, lt, pki, reach, m, deg);

      /* Offsets off the interpolation nodes, where the residual is zero by
       * construction. Runs start in [m, m + 1], anchored at the grid point
       * below the node; points reach |nx| = m + 2, the PRE_LIN_PSI table. */
      for (i = 0; i <= 512; i++)
      {
        const R nx0 = (R)m + (R)i / K(512.0);
        const R nx = (R)(m + 2) * ((R)i / K(256.0) - K(1.0));
        R err;

        Y(kb_poly_run)(got, coef, m, deg, nx0);

        for (l = 0; l < w; l++)
        {
          const R ref = Y(kb_phi)(b, lt, pki, reach, nx0 - (R)l);

          err = ABS(got[l] - ref) / peak;

          if (err > worst)
            worst = err;
        }

        err = ABS(Y(kb_poly_phi)(coef, m, deg, nx)
            - Y(kb_phi)(b, lt, pki, reach, nx)) / peak;

        if (err > worst)
          worst = err;
      }

      ok = IF(worst < bound && deg <= KB_POLY_DEG_MAX, 1, 0);
      printf("poly[sigma=" __FE__ ", m=%2td, deg=%2td] err_peak = " __FE__
          " %-2s " __FE__ " -> %-4s\n", sigma, (ptrdiff_t)m,
          (ptrdiff_t)deg, worst, IF(ok == 0, ">=", "<"), bound,
          IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);

      Y(free)(got);
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
\
  for (s = 0; s < KB_NSIGMA; s++) \
    for (m = 2; m <= KB_M_MAX; m = kb_m_next(m)) \
    { \
      int N[1] = {KB_N}, n[1] = {(int)kb_n[s] + (NADJ)}; \
      MOD(plan) pp, pa; \
      R ep, ea, uk, bound; \
      int ok; \
\
      MOD(init_guru)(&pp, 1, N, (int)M, n, (int)m, flags, \
          FFTW_ESTIMATE | FFTW_DESTROY_INPUT); \
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
      printf(NAME "[sigma=%5.2f, m=%2td, deg=%2td] poly %10.3e analytic " \
          "%10.3e ratio %6.3f, in u*kappa %7.2f %7.2f -> %-4s\n", \
          (double)((R)kb_n[s] / (R)KB_N), (ptrdiff_t)m, \
          (ptrdiff_t)(INT)pp.spline_coeffs[0], (double)ep, (double)ea, \
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

  printf("KAISER-BESSEL POLYNOMIAL TRAFO\n------------------------------\n");
  Y(srand48)(1234567890L);

  for (s = 0; s < KB_NSIGMA; s++)
    for (m = 2; m <= KB_M_MAX; m = kb_m_next(m))
    {
      int N[1] = {KB_N}, n[1] = {(int)kb_n[s]};
      NFFT(plan) pp, pa;
      R ep, ea;

      NFFT(init_guru)(&pp, 1, N, (int)M, n, (int)m, flags,
          FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
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

        printf("trafo[sigma=%5.2f, m=%2td, deg=%2td] poly %10.3e analytic "
            "%10.3e ratio %6.3f, in u*kappa %7.2f %7.2f -> %-4s\n",
            (double)((R)kb_n[s] / (R)KB_N), (ptrdiff_t)m,
            (ptrdiff_t)(INT)pp.spline_coeffs[0], (double)ep, (double)ea,
            (double)(ep / ea), (double)(ep / uk), (double)(ea / uk),
            IF(ok == 0, "FAIL", "OK"));
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

void X(check_kaiser_bessel_poly_sigma)(void)
{
  printf("KAISER-BESSEL POLYNOMIAL SIGMA\n-----------------------------\n"
      "skipped: window is not Kaiser-Bessel\n\n");
  CU_ASSERT(1);
}

#endif
