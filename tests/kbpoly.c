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

/* The polynomial's contract is absolute against the window peak, which is what
 * the psi table owes the transform. Relative error is not the yardstick out in
 * a tail that is 1E-19 of the peak, and holding it there is what sank the
 * earlier attempt at this. */
static R kb_peak(const R b, const R lg_tail, const R peak_inv, const R mr)
{
  return Y(kb_phi)(b, lg_tail, peak_inv, mr, K(0.0));
}

/* Approximation error at degree m + 6, measured in double. It is a property of
 * the interpolant, not of the precision, so it bounds every precision from
 * above; where the precision runs out first the eps floor takes over. */
static R approx_bound(const INT m)
{
  static const R tab[] =
  {
    K(0.0), K(0.0),
    K(7.1E-8),  /* m = 2 */
    K(4.2E-10), /* m = 3 */
    K(3.3E-11), /* m = 4 */
    K(2.1E-13), /* m = 5 */
    K(1.7E-14), /* m = 6 */
    K(6.0E-16), /* m = 7 */
    K(4.5E-16)  /* m = 8 */
  };

  return K(4.0) * tab[m];
}

static R eps_floor(const INT deg)
{
  /* The fit sum runs over deg + 1 terms and Horner adds a couple of units in
   * the last place on top. */
  return K(4.0) * ((R)deg + K(3.0)) * Y(float_property)(NFFT_EPSILON);
}

static const R kb_sigma[] = {K(1.25), K(2.0)};

void X(check_kaiser_bessel_poly)(void)
{
  unsigned int s;
  INT m;

  printf("KAISER-BESSEL POLYNOMIAL\n------------------------\n");

  for (s = 0; s < sizeof(kb_sigma) / sizeof(kb_sigma[0]); s++)
    for (m = 2; m <= 8; m++)
    {
      const INT N = 64, n = 2 * (INT)(kb_sigma[s] * (R)N / K(2.0));
      const R b = KPI * (K(2.0) - (R)N / (R)n);
      const R reach = (R)m + K(1.0);
      const R lt = Y(bessel_i0_logtail)(reach * b);
      const R pki = EXP(-reach * b - lt);
      const R peak = kb_peak(b, lt, pki, reach);
      const INT deg = Y(kb_poly_degree)(m), w = 2 * m + 2;
      const INT cols = KB_POLY_COLS(m);
      const R bound = IF(approx_bound(m) > eps_floor(deg), approx_bound(m),
          eps_floor(deg));
      R *coef = (R*) Y(malloc)((size_t)((deg + 1) * cols) * sizeof(R));
      R *got = (R*) Y(malloc)((size_t)w * sizeof(R));
      R worst = K(0.0);
      INT i, l;
      int ok;

      Y(kb_poly_fit)(coef, b, lt, pki, reach, m, deg);

      /* 513 offsets, off the interpolation nodes where the residual is zero
       * by construction. Runs start in [m - 1, m + 1], the reach of both run
       * anchors; points reach |nx| = m + 2, the PRE_LIN_PSI table. */
      for (i = 0; i <= 512; i++)
      {
        const R nx0 = (R)(m - 1) + (R)i / K(256.0);
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

      ok = IF(worst < bound, 1, 0);
      printf("poly[sigma=" __FE__ ", m=%2td, deg=%2td] err_peak = " __FE__
          " %-2s " __FE__ " -> %-4s\n", kb_sigma[s], (ptrdiff_t)m,
          (ptrdiff_t)deg, worst, IF(ok == 0, ">=", "<"), bound,
          IF(ok == 0, "FAIL", "OK"));
      CU_ASSERT(ok == 1);

      Y(free)(got);
      Y(free)(coef);
    }

  printf("\n");
}

#else

void X(check_kaiser_bessel_poly)(void)
{
  printf("KAISER-BESSEL POLYNOMIAL\n------------------------\nskipped: "
      "window is not Kaiser-Bessel\n\n");
  CU_ASSERT(1);
}

#endif
