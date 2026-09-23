/* nfft_sigma_m_1d against the legacy grid, 1-D, double, random input.
 *
 * Per "N M goal" line: (n_h, m_h) from nfft_sigma_m_1d at sigma in
 * [1.25, 2.5], then two arms on the same random nodes and coefficients, both
 * at m_h, flags as nfft_init_1d:
 *   helper   n = n_h
 *   legacy   n = 2 * next_power_of_2(N), the grid nfft_init_1d picks
 * Error: max_j |f_j - s_j| / ||f_hat||_1 over NCHECK nodes. Nodes are
 * multiples of 2^-30, so (k - N/2) q reduces exactly in 64-bit integers and
 * the double reference needs no long double. Time: best of 5 means of trafo
 * and of adjoint, each mean over at least 50 ms.
 *
 * usage: hostcheck < cases
 * rows:  N,M,goal,arm,n,m,sigma,err,t_trafo,t_adjoint */
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "nfft3.h"

#define NCHECK 100
#define QBITS 30
#define TWO_PI 6.283185307179586476925286766559

static double now(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

static double timeit(void (*f)(nfft_plan *), nfft_plan *p)
{
  double best = 1e30;
  int reps = 1, r, i;

  for (;;)
  {
    double t0 = now();
    for (i = 0; i < reps; i++)
      f(p);
    if (now() - t0 > 0.05)
      break;
    reps *= 2;
  }
  for (r = 0; r < 5; r++)
  {
    double t0 = now(), t;
    for (i = 0; i < reps; i++)
      f(p);
    t = (now() - t0) / reps;
    if (t < best)
      best = t;
  }
  return best;
}

static uint64_t rs = 12345u;
static uint64_t next(void)
{
  rs = rs * 6364136223846793005u + 1442695040888963407u;
  return rs >> 11;
}
static double urand(void)
{
  return (double)next() / 9007199254740992.0;
}

/* Neumaier summation, the running sum and its correction. */
static void add(double *s, double *c, double v)
{
  const double t = *s + v;
  *c += (fabs(*s) >= fabs(v)) ? (*s - t) + v : (v - t) + *s;
  *s = t;
}

static void arm(const char *name, int N, int M, double goal, int n, int m,
  const double *x, const double complex *fh, const int *chk,
  const double complex *ref, double l1)
{
  nfft_plan p;
  double err = 0.0;
  int j;

  nfft_init_guru(&p, 1, &N, M, &n, m,
      PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT
        | FFT_OUT_OF_PLACE,
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  for (j = 0; j < M; j++)
    p.x[j] = x[j];
  nfft_precompute_one_psi(&p);
  for (j = 0; j < N; j++)
    p.f_hat[j] = fh[j];
  nfft_trafo(&p);
  for (j = 0; j < NCHECK; j++)
  {
    const double e = cabs(ref[j] - p.f[chk[j]]);
    if (e > err)
      err = e;
  }
  printf("%d,%d,%.0e,%s,%d,%d,%.4f,%.3e,%.4e,%.4e\n", N, M, goal, name, n, m,
      (double)n / N, err / l1, timeit(nfft_trafo, &p),
      timeit(nfft_adjoint, &p));
  fflush(stdout);
  nfft_finalize(&p);
}

int main(void)
{
  const int64_t Q = (int64_t)1 << QBITS;
  int N, M;
  double goal;

  printf("# fftw %s\n", fftw_version);
  printf("N,M,goal,arm,n,m,sigma,err,t_trafo,t_adjoint\n");
  while (scanf("%d %d %lg", &N, &M, &goal) == 3)
  {
    double *x = malloc(sizeof(double) * (size_t)M);
    int64_t *q = malloc(sizeof(int64_t) * (size_t)M);
    double complex *fh = malloc(sizeof(double complex) * (size_t)N);
    double complex ref[NCHECK];
    int chk[NCHECK], n_h, m_h, n_l = 1, j, k;
    double l1 = 0.0;

    for (j = 0; j < M; j++)
    {
      q[j] = (int64_t)(next() % (uint64_t)Q) - Q / 2;
      x[j] = (double)q[j] / (double)Q;
    }
    for (k = 0; k < N; k++)
    {
      fh[k] = urand() + urand() * I;
      l1 += cabs(fh[k]);
    }
    for (j = 0; j < NCHECK; j++)
      chk[j] = (int)(next() % (uint64_t)M);

    /* s_j = sum_k fh_k exp(-2 pi i (k - N/2) x_j), phase index exact */
    for (j = 0; j < NCHECK; j++)
    {
      double sr = 0, cr = 0, si = 0, ci = 0;
      for (k = 0; k < N; k++)
      {
        const int64_t r = (((int64_t)(k - N / 2) * q[chk[j]]) % Q + Q) % Q;
        const double a = -TWO_PI * ((double)r / (double)Q);
        const double c = cos(a), s = sin(a);
        add(&sr, &cr, creal(fh[k]) * c - cimag(fh[k]) * s);
        add(&si, &ci, creal(fh[k]) * s + cimag(fh[k]) * c);
      }
      ref[j] = (sr + cr) + (si + ci) * I;
    }

    while (n_l < N)
      n_l *= 2;
    n_l *= 2;
    nfft_sigma_m_1d(N, M, goal, 1.25, 2.5, &n_h, &m_h);
    arm("helper", N, M, goal, n_h, m_h, x, fh, chk, ref, l1);
    arm("legacy", N, M, goal, n_l, m_h, x, fh, chk, ref, l1);

    free(x);
    free(q);
    free(fh);
  }
  return 0;
}
