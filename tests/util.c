/*
 * Copyright (c) 2016-2017, 2025-2026 The NFFT Authors
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

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"

static INT _log2i(const INT m)
{
  INT l = 0;
  INT mm = m;

  if (m <= 0)
    return -1;

  while (mm > (INT)(0))
  {
    mm = (mm >> 1);
    l++;
  }

  return (l-1);
}

void X(check_log2i)(void)
{
    INT i;
    INT j;

    {
      INT r = Y(log2i)(-1);
      int ok = r == -1;
      printf("log2i("__D__") = "__D__" -> %s\n", (INT)(-1), r, ok ? "OK" : "FAIL");
      CU_ASSERT(ok)
    }

    {
      INT r = Y(log2i)(0);
      int ok = r == -1;
      printf("log2i("__D__") = "__D__" -> %s\n", (INT)(0), r, ok ? "OK" : "FAIL");
      CU_ASSERT(ok)
    }

    {
      INT r = Y(log2i)(1);
      int ok = r == 0;
      printf("log2i("__D__") = "__D__" -> %s\n", (INT)(1), r, ok ? "OK" : "FAIL");
      CU_ASSERT(ok)
    }

    for (i = 0, j = 1; i < 8 * SIZEOF_PTRDIFF_T - 2; i++)
    {
      j <<= 1;
      {
        INT r = Y(log2i)(j);
        INT r2 = _log2i(j);
        int ok = r == r2;
        printf("log2i("__D__") = "__D__" -> %s\n", j, r, ok ? "OK" : "FAIL");
        CU_ASSERT(ok)
      }
      {
        INT r = Y(log2i)(j - 1);
        INT r2 = _log2i(j - 1);
        int ok = r == r2;
        printf("log2i("__D__") = "__D__" -> %s\n", j - 1, r, ok ? "OK" : "FAIL");
        CU_ASSERT(ok)
      }
    }
}

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
static INT _next_power_of_2(const INT N)
{
  INT n,i,logn;
  INT N_is_not_power_of_2=0;

  if (N == 0)
    return 1;
  else if (N == 1)
    return 2;
  else
  {
    n = N;
    logn = 0;
    while (n != 1)
    {
      if (n%2 == 1)
        N_is_not_power_of_2=1;
      n = n/2;
      logn++;
    }

    if (!N_is_not_power_of_2)
      logn--;

    for (i = 0; i <= logn; i++)
      n = n*2;

    return n;
  }
}

void X(check_next_power_of_2)(void)
{
    INT i;
    INT j;

    {
      INT r = Y(next_power_of_2)(-1);
      int ok = r == -1;
      printf("next_power_of_2("__D__") = "__D__" -> %s\n", (INT)(-1), r, ok ? "OK" : "FAIL");
      CU_ASSERT(ok)
    }

    {
      INT r = Y(next_power_of_2)(0);
      int ok = r == 1;
      printf("next_power_of_2("__D__") = "__D__" -> %s\n", (INT)(0), r, ok ? "OK" : "FAIL");
      CU_ASSERT(ok)
    }

    {
      INT r = Y(next_power_of_2)(1);
      int ok = r == 2;
      printf("next_power_of_2("__D__") = "__D__" -> %s\n", (INT)(1), r, ok ? "OK" : "FAIL");
      CU_ASSERT(ok)
    }

    for (i = 0, j = 1; i < 8 * SIZEOF_PTRDIFF_T - 2; i++)
    {
      j <<= 1;
      {
        INT r = Y(next_power_of_2)(j);
        INT r2 = _next_power_of_2(j);
        int ok = r == r2;
        printf("next_power_of_2("__D__") = "__D__" -> %s\n", j, r, ok ? "OK" : "FAIL");
        CU_ASSERT(ok)
      }
      {
        INT r = Y(next_power_of_2)(j - 1);
        INT r2 = _next_power_of_2(j - 1);
        int ok = r == r2;
        printf("next_power_of_2("__D__") = "__D__" -> %s\n", j - 1, r, ok ? "OK" : "FAIL");
        CU_ASSERT(ok)
      }
    }
}


/* Error bound for NFFT/NFCT/NFST using the Kaiser-bessel window.
 *
 * The error is window truncation plus round-off
 *
 * nfct and nfst run on the even or odd extension of f_hat with grid length
 * 2 (n + ext) and modes up to N - 1, against n and N / 2 for NFFT. */
R Y(kb_err_trafo)(INT d, INT m, const R *sigma, const INT *N, const INT *n,
    int ext)
{
  const R r = (R)m + K(1.0); // Stencil reach.
  const R eps = Y(float_property)(NFFT_EPSILON);
  R trunc = K(0.0);
  R spread = K(1.0); // Round-off amplification due to deconvolution.
  INT i;

  // All dimensions contribute.
  for (i = 0; i < d; i++)
  {
    const R b = KPI * (K(2.0) - K(1.0) / sigma[i]);
    const R z = K2PI * r * SQRT(K(1.0) - K(1.0) / sigma[i]); // Exponent of aliasing term.
    const R nfull = (R)(ext ? 2 * (n[i] + ext) : n[i]);
    const R kmax = (R)(ext ? N[i] - 1 : N[i] / 2); // Frequency index of band edge.

    // Truncation error.
    trunc += EXP(-z) / Y(bessel_i0_exp_scaled)(z);

    // Round-off amplification.
    spread *= Y(kb_phi_hut)(b, K(1.0), r, nfull, K(0.0))
        / Y(kb_phi_hut)(b, K(1.0), r, nfull, kmax);
  }

  // Total error term with some slack.
  return K(3.5) * (trunc + spread * eps);
}
