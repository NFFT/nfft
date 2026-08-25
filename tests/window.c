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

#include "nfft3.h"
#include "infft.h"
#include "window.h"

/* The window macros in infft.h fall through to Kaiser-Bessel. */
#if !defined(DIRAC_DELTA) && !defined(GAUSSIAN) && !defined(B_SPLINE) \
    && !defined(SINC_POWER)
#define WINDOW_IS_KAISER_BESSEL 1
#endif

#define ERR(x,y) IF(ABS(y) == K(0.0), ABS((x) - (y)), ABS((x) - (y)) / ABS(y))

/* The test build carries -ffast-math, which rules out isnan and isinf, so bound
 * the magnitude instead. Every peak-normalized window value is at most one. */
#define FINITE(v) (ABS(v) < K(1.0E30))

static R bound(void)
{
  return K(64.0) * Y(float_property)(NFFT_EPSILON);
}

/* Arguments straddling NFFT_I0_ASYMP_SPLIT in every precision, up to a peak
 * far beyond what I0 itself can represent. */
static const R i0_args[] =
{
  K(0.0), K(0.5), K(1.0), K(3.0), K(7.0), K(14.5), K(15.0), K(15.5),
  K(24.0), K(25.5), K(40.0), K(80.0), K(200.0)
};

/* Largest argument for which I0 is representable in single precision. */
#define I0_REPRESENTABLE K(80.0)

void X(check_bessel_i0_scaling)(void)
{
  const R b = bound();
  unsigned int j;

  printf("BESSEL I0 SCALING\n-----------------\n");

  for (j = 0; j < sizeof(i0_args) / sizeof(i0_args[0]); j++)
  {
    const R x = i0_args[j];
    const R lg = Y(bessel_i0_log)(x);
    const R lt = Y(bessel_i0_logtail)(x);
    const R sc = Y(bessel_i0_scaled)(x, lg);
    /* x - lg is formed inside bessel_i0_scaled, so its rounding grows with x */
    const R tol = b * (K(1.0) + x);
    R err;
    int ok;

    /* I0(x) exp(-log I0(x)) is one whatever x, and never forms I0(x) */
    err = ERR(sc, K(1.0));
    ok = IF(err < tol && FINITE(sc), 1, 0);
    printf("i0_scaled[" __FE__ "] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", x, sc, err, IF(ok == 0, ">=", "<"), tol,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);

    /* the tail agrees with the cancelling difference it replaces */
    err = ABS(lt - (lg - x));
    ok = IF(err < tol && FINITE(lt), 1, 0);
    printf("i0_logtail[" __FE__ "] = " __FE__ " err_abs = " __FE__ " %-2s " __FE__
      " -> %-4s\n", x, lt, err, IF(ok == 0, ">=", "<"), tol,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);

    if (x <= I0_REPRESENTABLE)
    {
      const R i0 = Y(bessel_i0)(x);

      err = ERR(Y(bessel_i0_scaled)(x, K(0.0)), i0);
      ok = IF(err < b, 1, 0);
      CU_ASSERT(ok == 1);

      err = ERR(EXP(lg), i0);
      ok = IF(err < tol, 1, 0);
      printf("i0_log[" __FE__ "] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
        " -> %-4s\n", x, lg, err, IF(ok == 0, ">=", "<"), tol,
        IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);
    }
  }
}

/* Oversampling factors and cut-offs. The large cut-offs put I0(m b) far above
 * the overflow threshold of every precision. */
static const R kb_sigma[] = {K(1.25), K(1.5), K(2.0)};
static const int kb_m[] = {2, 4, 6, 8, 12, 16, 20, 24};

void X(check_kaiser_bessel_peak)(void)
{
  const R b = bound();
  const INT N = 64;
  unsigned int s, i;

  printf("KAISER-BESSEL PEAK\n------------------\n");

  for (s = 0; s < sizeof(kb_sigma) / sizeof(kb_sigma[0]); s++)
  {
    for (i = 0; i < sizeof(kb_m) / sizeof(kb_m[0]); i++)
    {
      const int m = kb_m[i];
      const INT n = 2 * (INT)(kb_sigma[s] * (R)N / K(2.0));
      const R sh = KPI * (K(2.0) - (R)N / (R)n);
      const R lg = Y(bessel_i0_log)((R)m * sh);
      const R lt = Y(bessel_i0_logtail)((R)m * sh);
      const R peak = Y(kb_phi_hut)(sh, lg, lt, (R)m, (R)n, K(0.0));
      const R err = ERR(peak, K(1.0));
      const int ok = IF(err < b, 1, 0);
      R prev = peak;
      INT k;
      int l;

      printf("phi_hut[sigma=" __FE__ ", m=%2d, k=0] = " __FE__ " err_rel = "
        __FE__ " %-2s " __FE__ " -> %-4s\n", (R)n / (R)N, m, peak, err,
        IF(ok == 0, ">=", "<"), b, IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);

      /* phi_hut decays from the peak across the band and stays finite */
      for (k = 1; k <= N / 2; k++)
      {
        const R v = Y(kb_phi_hut)(sh, lg, lt, (R)m, (R)n, (R)k);
        const int okk = IF(FINITE(v) && v > K(0.0)
            && v <= prev * (K(1.0) + b), 1, 0);
        if (okk == 0)
          printf("phi_hut[sigma=" __FE__ ", m=%2d, k=" __D__ "] = " __FE__
            " -> FAIL\n", (R)n / (R)N, m, k, v);
        CU_ASSERT(okk == 1);
        prev = v;
      }

      /* phi peaks at the node and decays to zero across its support; at the
       * far edge the peak-normalized value is below the smallest number the
       * precision can hold, so it may be zero */
      prev = Y(kb_phi)(sh, lg, lt, (R)m, (R)n, K(0.0));
      CU_ASSERT(FINITE(prev) && prev > K(0.0));
      for (l = 1; l <= m; l++)
      {
        const R v = Y(kb_phi)(sh, lg, lt, (R)m, (R)n, (R)l / (R)n);
        const int okl = IF(FINITE(v) && v >= K(0.0)
            && v <= prev * (K(1.0) + b), 1, 0);
        if (okl == 0)
          printf("phi[sigma=" __FE__ ", m=%2d, l=%d] = " __FE__ " -> FAIL\n",
            (R)n / (R)N, m, l, v);
        CU_ASSERT(okl == 1);
        prev = v;
      }
    }
  }
}

void X(check_kaiser_bessel_reference)(void)
{
  /* m b stays well inside the range of I0 and sinh in single precision, so the
   * textbook forms below are usable as a reference. */
  const int m = 6;
  const INT N = 64, n = 128;
  const R sh = KPI * (K(2.0) - (R)N / (R)n);
  const R lg = Y(bessel_i0_log)((R)m * sh);
  const R lt = Y(bessel_i0_logtail)((R)m * sh);
  const R peak = Y(bessel_i0)((R)m * sh);
  const R b = K(4.0) * bound();
  INT k;
  int l;

  printf("KAISER-BESSEL REFERENCE\n-----------------------\n");

  for (k = 0; k <= N / 2; k++)
  {
    const R t = K(2.0) * KPI * (R)k / (R)n;
    const R a = (R)m * SQRT(sh * sh - t * t);
    const R ref = Y(bessel_i0)(a) / peak;
    const R v = Y(kb_phi_hut)(sh, lg, lt, (R)m, (R)n, (R)k);
    const R err = ERR(v, ref);
    const int ok = IF(err < b, 1, 0);
    printf("phi_hut[k=" __D__ "] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", k, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }

  for (l = 0; l <= m; l++)
  {
    const R x = (R)l / (R)n;
    const R s = (R)m * (R)m - ((R)n * x) * ((R)n * x);
    const R ref = IF(s > K(0.0),
        SINH(sh * SQRT(ABS(s))) / (KPI * SQRT(ABS(s))) / peak,
        (sh / KPI) / peak);
    const R v = Y(kb_phi)(sh, lg, lt, (R)m, (R)n, x);
    const R err = ERR(v, ref);
    const int ok = IF(err < b, 1, 0);
    printf("phi[l=%d] = " __FE__ " err_rel = " __FE__ " %-2s " __FE__
      " -> %-4s\n", l, v, err, IF(ok == 0, ">=", "<"), b,
      IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);
  }
}

void X(check_kaiser_bessel_nfft)(void)
{
#if defined(WINDOW_IS_KAISER_BESSEL)
  /* Cut-offs for which the unnormalized I0(m b) raised to the power d leaves
   * the range of single precision. */
  static const int cut_off[] = {6, 8, 10, 12, 14};
  const int d = 3, M = 24;
  int N[3] = {16, 16, 16}, n[3] = {32, 32, 32};
  unsigned int i;

  printf("KAISER-BESSEL NFFT\n------------------\n");

  for (i = 0; i < sizeof(cut_off) / sizeof(cut_off[0]); i++)
  {
    const R sigma = (R)n[0] / (R)N[0];
    /* roundoff plus the truncation error of the window itself */
    const R b = K(1000.0) * Y(float_property)(NFFT_EPSILON)
        + K(4.0) * EXP(K(-2.0) * KPI * (R)cut_off[i]
            * SQRT(K(1.0) - K(1.0) / sigma));
    Y(plan) p;
    C *ref;
    R err;
    int j, ok;

    Y(init_guru)(&p, d, N, M, n, cut_off[i],
        PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    Y(vrand_shifted_unit_double)(p.x, p.d * p.M_total);
    Y(vrand_unit_complex)((C*)p.f_hat, p.N_total);
    Y(precompute_one_psi)(&p);

    ref = (C*)Y(malloc)((size_t)p.M_total * sizeof(C));
    Y(trafo_direct)(&p);
    for (j = 0; j < p.M_total; j++)
      ref[j] = ((C*)p.f)[j];

    Y(trafo)(&p);
    err = Y(error_l_infty_1_complex)(ref, (C*)p.f, p.M_total, (C*)p.f_hat,
        p.N_total);
    ok = IF(FINITE(err) && err < b, 1, 0);
    printf("nfft_3d[m=%2d] err = " __FE__ " %-2s " __FE__ " -> %-4s\n",
      cut_off[i], err, IF(ok == 0, ">=", "<"), b, IF(ok == 0, "FAIL", "OK"));
    CU_ASSERT(ok == 1);

    Y(free)(ref);
    Y(finalize)(&p);
  }
#else
  printf("KAISER-BESSEL NFFT\n------------------\nskipped: %s window\n",
    Y(get_window_name)());
#endif
}
