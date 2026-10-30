/* Measure the NFFT error with the stock polynomial Kaiser-Bessel window and
 * with replacement coefficient tables written by kbopt5.py / kbopt7.py.
 *
 *   kbgain2 N M sigma m seed [coef.txt ...]
 *
 * The reference is a long double direct sum with exact argument reduction
 * (k x is exact in long double), so it is good to ~1e-18 relative and the
 * NFFT's own rounding floor is visible. For each table: random f_hat
 * (max |f - f_direct| / sum |f_hat| and the RMS estimate of J_2) and single
 * frequencies k (max over nodes of |f_j - e^{-2 pi i k x_j}| = max_t |E(k/n, t)|).
 *
 * A coefficient file may end with a block "C <N>" followed by N values of
 * psi_hat(k/n), k = -N/2 .. N/2-1; they replace the plan's 1/PHI_HUT table.
 * Compile with -DFLOATPREC against libnfft3f for the float build. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "nfft3.h"

#ifdef FLOATPREC
  typedef float R;
  #define PLAN nfftf_plan
  #define INIT nfftf_init_guru
  #define PRE nfftf_precompute_one_psi
  #define TRAFO nfftf_trafo
  #define FINAL nfftf_finalize
#else
  typedef double R;
  #define PLAN nfft_plan
  #define INIT nfft_init_guru
  #define PRE nfft_precompute_one_psi
  #define TRAFO nfft_trafo
  #define FINAL nfft_finalize
#endif

static const long double TWOPI = 6.283185307179586476925286766559L;

static long double complex ref_exp(int k, R x)
{
  long double kx = (long double)k * (long double)x;
  long double f = kx - floorl(kx);
  return cosl(TWOPI * f) - I * sinl(TWOPI * f);
}

static void measure(PLAN *p, const char *label, int N)
{
  int j, k, i;
  long double sum_hat = 0, sum_hat2 = 0, err_inf = 0, err2 = 0;
  long double complex *fd = (long double complex*) malloc((size_t)p->M_total * sizeof(long double complex));
  const int ks[] = {0, N / 32, N / 16, N / 8, N / 4, 3 * N / 8, N / 2 - 1};
  const int nks = (int)(sizeof ks / sizeof ks[0]);

  PRE(p);
  srand(12345);
  for (k = 0; k < N; k++)
  {
    p->f_hat[k] = (R)((double)rand() / RAND_MAX - 0.5) + I * (R)((double)rand() / RAND_MAX - 0.5);
    sum_hat += cabsl(p->f_hat[k]);
    sum_hat2 += cabsl(p->f_hat[k]) * cabsl(p->f_hat[k]);
  }
  for (j = 0; j < p->M_total; j++)
  {
    long double complex s = 0;
    for (k = 0; k < N; k++)
      s += (long double complex)p->f_hat[k] * ref_exp(k - N / 2, p->x[j]);
    fd[j] = s;
  }
  TRAFO(p);
  for (j = 0; j < p->M_total; j++)
  {
    long double e = cabsl((long double complex)p->f[j] - fd[j]);
    if (e > err_inf) err_inf = e;
    err2 += e * e;
  }
  printf("%-8s random: max|f-fd|/sum|fhat| = %.3Le   rms est. of J_2 = %.3Le\n",
      label, err_inf / sum_hat, sqrtl(err2 / (p->M_total * sum_hat2)));
  printf("%-8s single k:", label);
  for (i = 0; i < nks; i++)
  {
    long double mx = 0;
    for (k = 0; k < N; k++)
      p->f_hat[k] = 0;
    p->f_hat[ks[i] + N / 2] = 1;
    TRAFO(p);
    for (j = 0; j < p->M_total; j++)
    {
      long double e = cabsl((long double complex)p->f[j] - ref_exp(ks[i], p->x[j]));
      if (e > mx) mx = e;
    }
    printf("  k=%5d: %.3Le", ks[i], mx);
  }
  printf("\n");
  free(fd);
}

int main(int argc, char **argv)
{
  PLAN p;
  int N = atoi(argv[1]), M = atoi(argv[2]), m = atoi(argv[4]), a, j;
  double sigma = atof(argv[3]);
  unsigned seed = (unsigned)atoi(argv[5]);
  int n = (int)lrint(sigma * N);
  int deg, rows, cols;

  if (n % 2) n++;
  INIT(&p, 1, &N, M, &n, m,
      PRE_PHI_HUT | PRE_PSI | (getenv("OPT") ? OPT_WINDOW : 0) | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  srand(seed);
  for (j = 0; j < p.M_total; j++)
    p.x[j] = (R)((double)rand() / RAND_MAX - 0.5);
  deg = (int)p.spline_coeffs[0];
  rows = 2 * ((deg + 1) / 2) + 1;
  cols = m + 2;
  printf("N=%d M=%d n=%d sigma=%g m=%d  table degree %d (rows %d, cols %d)  %s\n",
      N, M, n, sigma, m, deg, rows, cols, sizeof(R) == 4 ? "float" : "double");
  measure(&p, "stock", N);

  for (a = 6; a < argc; a++)
  {
    FILE *f = fopen(argv[a], "r");
    int fdeg, fm, k, c, fN;
    char tag[8];
    if (!f || fscanf(f, "%d %d", &fdeg, &fm) != 2 || fm != m)
    {
      fprintf(stderr, "bad file %s\n", argv[a]);
      return 1;
    }
    if (2 * ((fdeg + 1) / 2) + 1 != rows)
    {
      rows = 2 * ((fdeg + 1) / 2) + 1;
      free(p.spline_coeffs);
      p.spline_coeffs = (R*) calloc((size_t)(1 + rows * cols), sizeof(R));
    }
    p.spline_coeffs[0] = (R)fdeg;
    for (k = 0; k < rows; k++)
      for (c = 1; c <= m + 1; c++)
      {
        double v;
        if (fscanf(f, "%lf", &v) != 1)
        {
          fprintf(stderr, "short file %s\n", argv[a]);
          return 1;
        }
        p.spline_coeffs[1 + k * cols + c] = (R)v;
      }
    if (fscanf(f, "%7s %d", tag, &fN) == 2 && !strcmp(tag, "C") && fN == N)
    {
      for (k = 0; k < N; k++)
      {
        double v;
        if (fscanf(f, "%lf", &v) != 1)
        {
          fprintf(stderr, "short C block %s\n", argv[a]);
          return 1;
        }
        p.c_phi_inv[0][k] = (R)(1.0 / v);
      }
      printf("--- %s (degree %d, own deconvolution)\n", argv[a], fdeg);
    }
    else
      printf("--- %s (degree %d)\n", argv[a], fdeg);
    fclose(f);
    {
      const char *base = strrchr(argv[a], '/');
      char label[64];
      base = base ? base + 1 : argv[a];
      snprintf(label, sizeof label, "%.*s", (int)strcspn(base + 5, "_"), base + 5);
      measure(&p, label, N);
    }
  }
  FINAL(&p);
  return 0;
}
