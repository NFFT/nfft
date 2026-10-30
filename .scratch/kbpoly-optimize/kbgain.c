/* Measure the NFFT error with the stock polynomial Kaiser-Bessel window and
 * with a replacement coefficient table written by kbopt5.py.
 *
 *   kbgain N M sigma m seed [coef.txt ...]
 *
 * For each table: random f_hat (max |f - f_direct| / sum |f_hat|, and the
 * RMS estimate of J_2) and single frequencies k (max over nodes of
 * |f_j - e^{2 pi i k x_j}|, which is max_t |E(k/n, t)|). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "nfft3.h"

static void nodes(nfft_plan *p, unsigned seed)
{
  int j;
  srand(seed);
  for (j = 0; j < p->M_total; j++)
    p->x[j] = (double)rand() / RAND_MAX - 0.5;
}

static void measure(nfft_plan *p, const char *label, int N)
{
  int j, k, i;
  double sum_hat = 0, sum_hat2 = 0, err_inf = 0, err2 = 0;
  double complex *fd = (double complex*) malloc((size_t)p->M_total * sizeof(double complex));
  const int ks[] = {0, N / 8, N / 4, N / 2 - 1};

  nfft_precompute_one_psi(p);

  /* random f_hat */
  srand(12345);
  for (k = 0; k < N; k++)
  {
    p->f_hat[k] = ((double)rand() / RAND_MAX - 0.5) + I * ((double)rand() / RAND_MAX - 0.5);
    sum_hat += cabs(p->f_hat[k]);
    sum_hat2 += cabs(p->f_hat[k]) * cabs(p->f_hat[k]);
  }
  nfft_trafo_direct(p);
  memcpy(fd, p->f, (size_t)p->M_total * sizeof(double complex));
  nfft_trafo(p);
  for (j = 0; j < p->M_total; j++)
  {
    double e = cabs(p->f[j] - fd[j]);
    if (e > err_inf) err_inf = e;
    err2 += e * e;
  }
  printf("%-8s random: max|f-fd|/sum|fhat| = %.3e   rms est. of J_2 = %.3e\n",
      label, err_inf / sum_hat, sqrt(err2 / (p->M_total * sum_hat2)));

  /* single frequencies: the error at node j is exactly |E(k/n, t_j)| */
  printf("%-8s single k:", label);
  for (i = 0; i < 4; i++)
  {
    double mx = 0;
    for (k = 0; k < N; k++)
      p->f_hat[k] = 0;
    p->f_hat[ks[i] + N / 2] = 1;           /* f_hat index k + N/2 is frequency k */
    nfft_trafo(p);
    for (j = 0; j < p->M_total; j++)
    {
      double e = cabs(p->f[j] - cexp(-2 * M_PI * I * ks[i] * p->x[j]));
      if (e > mx) mx = e;
    }
    printf("  k=%5d: %.3e", ks[i], mx);
  }
  printf("\n");
  free(fd);
}

int main(int argc, char **argv)
{
  nfft_plan p;
  int N = atoi(argv[1]), M = atoi(argv[2]), m = atoi(argv[4]), a;
  double sigma = atof(argv[3]);
  unsigned seed = (unsigned)atoi(argv[5]);
  int n = (int)lrint(sigma * N);
  int deg, rows, cols;

  if (n % 2) n++;
  nfft_init_guru(&p, 1, &N, M, &n, m,
      PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  nodes(&p, seed);
  deg = (int)p.spline_coeffs[0];
  rows = 2 * ((deg + 1) / 2) + 1;
  cols = m + 2;
  printf("N=%d M=%d n=%d sigma=%g m=%d  table degree %d (rows %d, cols %d)\n",
      N, M, n, sigma, m, deg, rows, cols);
  measure(&p, "stock", N);

  for (a = 6; a < argc; a++)
  {
    FILE *f = fopen(argv[a], "r");
    int fdeg, fm, k, c;
    if (!f || fscanf(f, "%d %d", &fdeg, &fm) != 2 || fm != m)
    {
      fprintf(stderr, "bad file %s\n", argv[a]);
      return 1;
    }
    if (2 * ((fdeg + 1) / 2) + 1 != rows)
    {
      /* a different degree: rebuild the table with the file's degree */
      rows = 2 * ((fdeg + 1) / 2) + 1;
      free(p.spline_coeffs);
      p.spline_coeffs = (double*) calloc((size_t)(1 + rows * cols), sizeof(double));
      /* column 0 (cell [m+1, m+2)) is outside the run: leave it zero */
    }
    p.spline_coeffs[0] = (double)fdeg;
    for (k = 0; k < rows; k++)
      for (c = 1; c <= m + 1; c++)
        if (fscanf(f, "%lf", &p.spline_coeffs[1 + k * cols + c]) != 1)
        {
          fprintf(stderr, "short file %s\n", argv[a]);
          return 1;
        }
    fclose(f);
    {
      const char *base = strrchr(argv[a], '/');
      char label[64];
      base = base ? base + 1 : argv[a];
      snprintf(label, sizeof label, "%.*s", (int)strcspn(base + 5, "_"), base + 5);
      printf("--- %s (degree %d)\n", argv[a], fdeg);
      measure(&p, label, N);
    }
  }
  nfft_finalize(&p);
  return 0;
}
