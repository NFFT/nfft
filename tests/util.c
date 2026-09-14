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
#include <complex.h> /* before nfft3.h so fftw_complex is C-compatible */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <CUnit/CUnit.h>

#include "config.h" /* ABS_SRCDIR */
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "util.h"

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

int Y(test_read_case)(const char *rel, int *d, INT **N, INT *NN, INT *M, R **x,
                      C **f_hat, C **f)
{
  char path[4096];
  FILE *fp;
  int t;
  long v;
  INT j, nn = 1;

  *N = NULL;
  *x = NULL;
  *f_hat = NULL;
  *f = NULL;

  snprintf(path, sizeof path, "%s/tests/%s", ABS_SRCDIR, rel);
  fp = fopen(path, "r");
  if (!fp)
    return 0;

  if (fscanf(fp, "%d", d) != 1 || *d < 1)
    goto fail;
  *N = (INT *)Y(malloc)((size_t)*d * sizeof(INT));
  for (t = 0; t < *d; t++)
  {
    if (fscanf(fp, "%ld", &v) != 1)
      goto fail;
    (*N)[t] = (INT)v;
    nn *= (INT)v;
  }
  if (fscanf(fp, "%ld", &v) != 1)
    goto fail;
  *M = (INT)v;
  *NN = nn;

  *x = (R *)Y(malloc)((size_t)(*d * *M) * sizeof(R));
  for (j = 0; j < *d * *M; j++)
  {
    if (fscanf(fp, __FI__, &(*x)[j]) != 1)
      goto fail;
  }
  *f_hat = (C *)Y(malloc)((size_t)nn * sizeof(C));
  for (j = 0; j < nn; j++)
  {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2)
      goto fail;
    (*f_hat)[j] = re + II * im;
  }
  *f = (C *)Y(malloc)((size_t)*M * sizeof(C));
  for (j = 0; j < *M; j++)
  {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2)
      goto fail;
    (*f)[j] = re + II * im;
  }
  fclose(fp);
  return 1;

fail:
  fclose(fp);
  Y(free)(*N);
  Y(free)(*x);
  Y(free)(*f_hat);
  Y(free)(*f);
  *N = NULL;
  *x = NULL;
  *f_hat = NULL;
  *f = NULL;
  return 0;
}

R Y(test_rel_max_err)(const C *a, const C *b, INT len)
{
  R num = (R)0, den = (R)0;
  INT j;
  for (j = 0; j < len; j++)
  {
    R e = CABS(a[j] - b[j]);
    if (e > num)
      num = e;
    if (CABS(b[j]) > den)
      den = CABS(b[j]);
  }
  return den > (R)0 ? num / den : num;
}

/* The a/b calibration follows err_trafo in tests/nfft.c; m and sigma are
 * runtime, a and b are precision-fixed. The 56*eps floor covers the fast
 * pipeline's DECONV+CONV round-off. */
R Y(test_err_bound)(int window, R m, R s)
{
  R eps = Y(float_property)(NFFT_EPSILON), a, b, err;
  switch (window)
  {
  case NFFT_WINDOW_GAUSSIAN:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(0.41);
    b = K(50.0);
#else
    a = K(0.95);
    b = K(50.0);
#endif
    err = EXP(-m * KPI * (K(1.0) - K(1.0) / (K(2.0) * K(2.0) - K(1.0))));
    break;
  case NFFT_WINDOW_B_SPLINE:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(1.0);
    b = K(2000.0);
#else
    a = K(0.3);
    b = K(50.0);
#endif
    err = K(3000.0) * K(4.0) * POW(K(1.0) / (K(2.0) * s - K(1.0)), K(2.0) * m);
    break;
  case NFFT_WINDOW_SINC_POWER:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(1.0);
    b = K(2000.0);
#else
    a = K(0.3);
    b = K(50.0);
#endif
    err = (K(1.0) / (m - K(1.0))) * ((K(2.0) / POW(s, K(2.0) * m)) + POW(s / (K(2.0) * s - K(1.0)), K(2.0) * m));
    break;
  case NFFT_WINDOW_KAISER_BESSEL:
  default:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    /* The eps floor must stay above the error check_nfast_float_accuracy
     * measures on the 3D reference case. */
    a = K(0.3);
    b = K(3000.0);
#else
    a = K(1.5);
    b = K(50.0);
#endif
    err = KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0) / K(2.0))) * EXP(-K2PI * m * SQRT(K(1.0) - K(1.0) / K(2.0)));
    break;
  }
  return FMAX(FMAX(a * err, b * eps), K(56.0) * eps);
}

void Y(test_assert_plan_names)(struct Y(plan_ng_s) *p, const char *needle)
{
  FILE *tmp = tmpfile();
  long len;
  char *buf;

  CU_ASSERT_PTR_NOT_NULL_FATAL(tmp);
  Y(fprint_plan)((Y(plan_ng) *)p, tmp);
  fseek(tmp, 0, SEEK_END);
  len = ftell(tmp);
  rewind(tmp);
  buf = (char *)Y(malloc)((size_t)len + 1);
  if (len > 0)
    CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tmp), (size_t)len);
  buf[len] = '\0';
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, needle));
  Y(free)(buf);
  fclose(tmp);
}

char *Y(test_wisdom_export)(planner *pl)
{
  size_t cnt;
  char *s;
  printer *p = Y(printer_create_cnt)(&cnt);
  Y(planner_export)(pl, p);
  Y(printer_destroy)(p);
  s = (char *)Y(malloc)(cnt + 1);
  p = Y(printer_create_str)(s);
  Y(planner_export)(pl, p);
  Y(printer_destroy)(p);
  return s;
}

int Y(test_wisdom_import)(planner *pl, const char *s)
{
  scanner *sc = Y(scanner_create_str)(s);
  int ret = Y(planner_import)(pl, sc);
  Y(scanner_destroy)(sc);
  return ret;
}

