/* Acceptance of the fractional minimax correction (Y(kb_poly_correct)).
 *   kbfrac rows               the reference rows, from scratch
 *   kbfrac cell sigma m reps  n = sigma 4096 even, library degree: from
 *                             scratch against the cell start of the lookup
 * hut: the transform's functional of tests/kbpoly.c on 401 x 64, C from the
 * hut table. exact: D and C in long double from the stored coefficients, C
 * the window's own spectrum by Gauss-Legendre. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "nfft3.h"
#include "infft.h"

typedef long double L;

static double now(void)
{
  struct timespec t;

  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

#define KB_FN_NXI 401
#define KB_FN_NT 64

/* tests/kbpoly.c, finer grid. */
static void functional(const R *coef, const R *hut, const INT m,
    const INT deg, const R b, const R lt, const R pki, const R i0e_pki,
    R *j_inf, R *j_2, R *poly_part)
{
  const INT w = 2 * m + 2;
  R ss = K(0.0);
  const R reach = (R)m + K(1.0), sigma = KPI / (K2PI - b);
  R *got = (R*) Y(malloc)((size_t)(2 * w) * sizeof(R)), *ref = got + w;
  INT i, q, l;

  *j_inf = *poly_part = K(0.0);

  for (i = 0; i < KB_FN_NXI; i++)
  {
    const R xi = (R)i / (R)(KB_FN_NXI - 1) / (K(2.0) * sigma);
    const R C = Y(kb_poly_phi_hut)(hut, b, i0e_pki, reach, K(1.0), xi);

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
      ss += e * e;
      e = SQRT(ere * ere + eim * eim) / C;
      *poly_part = IF(e > *poly_part, e, *poly_part);
    }
  }

  *j_2 = SQRT(ss / (R)(KB_FN_NXI * KB_FN_NT));
  Y(free)(got);
}

#define NG 40
static L gx[NG], gw[NG];

static void gauleg(void)
{
  const L pi = 3.14159265358979323846264338327950288L;
  int i, j, k;

  for (i = 0; i < NG; i++)
  {
    L x = cosl(pi * (i + 0.75L) / (NG + 0.5L)), p0, p1, dp;

    for (k = 0; k < 100; k++)
    {
      L dx;

      p0 = 1.0L;
      p1 = x;
      for (j = 2; j <= NG; j++)
      {
        const L p2 = ((2 * j - 1) * x * p1 - (j - 1) * p0) / j;

        p0 = p1;
        p1 = p2;
      }
      dp = NG * (x * p1 - p0) / (x * x - 1.0L);
      dx = p1 / dp;
      x -= dx;
      if (fabsl(dx) < 1e-33L)
        break;
    }
    gx[i] = x / 2;
    gw[i] = 1.0L / ((1.0L - x * x) * dp * dp);
  }
}

/* phi(|x|) in cell column c, s from the cell centre. */
static L cell_poly(const R *coef, const INT m, const INT deg, const INT c,
    const L s)
{
  const INT cs = KB_POLY_COLS(m);
  L v = 0.0L;
  INT k;

  for (k = KB_POLY_ROWS(deg) - 1; k >= 0; k--)
    v = v * s + (L)coef[k * cs + c];

  return v;
}

static L phi_l(const R *coef, const INT m, const INT deg, const L x)
{
  const L a = fabsl(x);
  const INT f = (INT)a;

  if (f > m)
    return 0.0L;

  return cell_poly(coef, m, deg, m + 1 - f, a - (L)f - 0.5L);
}

static R exact_max(const R *coef, const INT m, const INT deg, const R b)
{
  const L pi2 = 6.28318530717958647692528676655900577L;
  const L sigma = (L)KPI / ((L)K2PI - (L)b);
  L worst = 0.0L;
  INT i, q, l, c, g;

  for (i = 0; i < KB_FN_NXI; i++)
  {
    const L xi = (L)i / (L)(KB_FN_NXI - 1) / (2.0L * sigma);
    L C = 0.0L;

    for (c = 1; c <= m + 1; c++)
      for (g = 0; g < NG; g++)
        C += 2.0L * gw[g] * cell_poly(coef, m, deg, c, gx[g])
            * cosl(pi2 * xi * ((L)(m + 1 - c) + 0.5L + gx[g]));

    for (q = 0; q < KB_FN_NT; q++)
    {
      const L nx0 = (L)m + ((L)q + 0.5L) / (L)KB_FN_NT;
      L dre = 0.0L, dim = 0.0L, e;

      for (l = 0; l < 2 * m + 2; l++)
      {
        const L v = phi_l(coef, m, deg, nx0 - (L)l), ph = -pi2 * xi
            * (nx0 - (L)l);

        dre += v * cosl(ph);
        dim += v * sinl(ph);
      }

      e = sqrtl((dre / C - 1.0L) * (dre / C - 1.0L) + dim * dim / (C * C));
      worst = e > worst ? e : worst;
    }
  }

  return (R)worst;
}

typedef struct
{
  R b, lt, pki, i0e;
  INT m, deg, n;
  R *plain, *coef, *hut;
} win;

static void win_init(win *wn, const R sigma, const INT m, const INT deg)
{
  const R reach = (R)m + K(1.0);

  wn->b = KPI * (K(2.0) - K(1.0) / sigma);
  wn->lt = Y(bessel_i0_logtail)(reach * wn->b);
  wn->pki = EXP(-reach * wn->b - wn->lt);
  wn->i0e = K(1.0) / Y(bessel_i0_exp_scaled)(reach * wn->b);
  wn->m = m;
  wn->deg = deg;
  wn->n = KB_POLY_ROWS(deg) * KB_POLY_COLS(m);
  wn->plain = (R*) Y(malloc)((size_t)(2 * wn->n + KB_POLY_HUT_LEN)
      * sizeof(R));
  wn->coef = wn->plain + wn->n;
  wn->hut = wn->coef + wn->n;
  Y(kb_poly_fit)(wn->plain, wn->b, wn->lt, wn->pki, reach, m, deg);
}

/* One correction from the plain fit; returns the wall time in seconds. */
static double win_correct(win *wn, const kb_poly_opt_cell *cell,
    const kb_poly_opt_cell *above, R *ret)
{
  double t0;

  memcpy(wn->coef, wn->plain, (size_t)wn->n * sizeof(R));
  wn->hut[0] = K(0.0);
  t0 = now();
  *ret = Y(kb_poly_correct)(wn->coef, wn->hut, wn->b, wn->m, wn->deg, 1,
      cell, above);

  return now() - t0;
}

static R win_hut(const win *wn)
{
  R j_inf, j_2, part;

  functional(wn->coef, wn->hut, wn->m, wn->deg, wn->b, wn->lt, wn->pki,
      wn->i0e, &j_inf, &j_2, &part);

  return j_inf;
}

static int cmp(const void *a, const void *b)
{
  const double x = *(const double*)a, y = *(const double*)b;

  return (x > y) - (x < y);
}

/* Median, min, max in ms of reps corrections. */
static void timed(win *wn, const kb_poly_opt_cell *cell,
    const kb_poly_opt_cell *above, const int reps, double *st, R *ret)
{
  double *t = (double*) malloc((size_t)reps * sizeof(double));
  int r;

  for (r = 0; r < reps; r++)
    t[r] = 1e3 * win_correct(wn, cell, above, ret);

  qsort(t, (size_t)reps, sizeof(double), cmp);
  st[0] = reps % 2 ? t[reps / 2] : 0.5 * (t[reps / 2 - 1] + t[reps / 2]);
  st[1] = t[0];
  st[2] = t[reps - 1];
  free(t);
}

/* sigma, paper m, deg, e of the reference. */
static const double rows[][4] = {
  {1.125, 2, 8, 1.987e-02}, {1.125, 3, 9, 1.260e-03},
  {1.125, 4, 10, 7.857e-05}, {1.125, 5, 11, 4.885e-06},
  {1.125, 6, 12, 3.004e-07}, {1.125, 7, 13, 1.966e-08},
  {1.25, 2, 8, 7.567e-03}, {1.25, 3, 9, 2.737e-04},
  {1.25, 4, 10, 1.058e-05}, {1.25, 5, 11, 3.453e-07},
  {1.25, 6, 12, 1.309e-08}, {1.25, 7, 13, 4.869e-10},
  {1.5, 2, 8, 2.135e-03}, {1.5, 3, 9, 4.016e-05},
  {1.5, 4, 10, 6.496e-07}, {1.5, 5, 11, 1.183e-08},
  {1.5, 6, 12, 2.174e-10}, {1.5, 7, 13, 3.955e-12},
  {1.5, 6, 8, 1.593e-08},
  {2.0, 2, 8, 4.055e-04}, {2.0, 3, 9, 3.590e-06},
  {2.0, 4, 10, 2.183e-08}, {2.0, 5, 11, 1.613e-10},
  {2.0, 6, 12, 1.260e-12}};

static void run_rows(void)
{
  unsigned i;

  printf("%6s %2s %3s %10s %10s %6s %10s %6s %10s %9s\n", "sigma", "pm",
      "deg", "ref e", "hut", "hut/e", "exact", "ex/e", "return", "ms");

  for (i = 0; i < sizeof(rows) / sizeof(rows[0]); i++)
  {
    win wn;
    R ret, hut, ex;
    double t;

    win_init(&wn, (R)rows[i][0], (INT)rows[i][1] - 1, (INT)rows[i][2]);
    t = win_correct(&wn, NULL, NULL, &ret);
    hut = win_hut(&wn);
    ex = exact_max(wn.coef, wn.m, wn.deg, wn.b);
    printf("%6.3f %2d %3d %10.3e %10.3e %6.3f %10.3e %6.3f %10.3e %9.1f\n",
        rows[i][0], (int)rows[i][1], (int)rows[i][2], rows[i][3],
        (double)hut, (double)hut / rows[i][3], (double)ex,
        (double)ex / rows[i][3], (double)ret, 1e3 * t);
    fflush(stdout);
    Y(free)(wn.plain);
  }
}

static void run_cell(const double sig, const INT m, const int reps)
{
  const INT N = 4096;
  INT n = (INT)(sig * N);
  R sigma, b, ret;
  INT deg;
  const kb_poly_opt_cell *exact, *below, *above;
  double st[3];
  win wn;

  n += n % 2;
  sigma = (R)n / (R)N;
  b = KPI * (K(2.0) - K(1.0) / sigma);
  deg = Y(kb_poly_degree)(b, (R)m + K(1.0), m);
  win_init(&wn, sigma, m, deg);
  Y(kb_poly_table)(b, m, deg, &exact, &below, &above);

  timed(&wn, NULL, NULL, reps, st, &ret);
  printf("sigma=%.5f m=%d deg=%d ok=%d scratch: median %.2f min %.2f max %.2f"
      " ms, return %.3e hut %.3e exact %.3e", (double)sigma, (int)m, (int)deg,
      Y(kb_poly_correct_ok)(b, m), st[0], st[1], st[2], (double)ret,
      (double)win_hut(&wn), (double)exact_max(wn.coef, m, deg, b));

  if (exact != NULL)
    printf(" | table e %.3e scratch/e %.3f\n", (double)exact->e,
        (double)win_hut(&wn) / (double)exact->e);
  else if (below != NULL)
  {
    const double sc = st[0];

    timed(&wn, below, above, reps, st, &ret);
    printf(" | start below %.3f above %.3f: median %.2f min %.2f max %.2f ms"
        " (%.2fx), return %.3e hut %.3e\n", (double)below->sigma,
        above ? (double)above->sigma : 0.0, st[0], st[1], st[2], st[0] / sc,
        (double)ret, (double)win_hut(&wn));
  }
  else
    printf(" | no cell\n");

  Y(free)(wn.plain);
}

int main(int argc, char **argv)
{
  gauleg();

  if (argc == 2 && strcmp(argv[1], "rows") == 0)
    run_rows();
  else if (argc == 5 && strcmp(argv[1], "cell") == 0)
    run_cell(atof(argv[2]), (INT)atoi(argv[3]), atoi(argv[4]));
  else
  {
    fprintf(stderr, "usage: kbfrac rows | kbfrac cell sigma m reps\n");
    return 1;
  }

  return 0;
}
