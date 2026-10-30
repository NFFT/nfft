/* Time nfft_trafo / nfft_adjoint and nfft_precompute_one_psi with the stock
 * polynomial window table and with replacement tables of lower degree.
 *
 *   kbtime d N M sigma m mode reps [coef.txt ...]     mode: otf | pre
 *
 * otf: window evaluated on the fly in the transform (no PRE_PSI).
 * pre: PRE_PSI; the window is evaluated in precompute_one_psi only.
 * Prints the best of reps wall-clock times. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "nfft3.h"

static double now(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void timeit(nfft_plan *p, const char *label, int reps)
{
  int r, j;
  double tp = 1e9, tt = 1e9, ta = 1e9;
  for (r = 0; r < reps; r++)
  {
    double t0 = now();
    nfft_precompute_one_psi(p);
    double t1 = now();
    nfft_trafo(p);
    double t2 = now();
    nfft_adjoint(p);
    double t3 = now();
    if (t1 - t0 < tp) tp = t1 - t0;
    if (t2 - t1 < tt) tt = t2 - t1;
    if (t3 - t2 < ta) ta = t3 - t2;
  }
  for (j = 0; j < 3; j++)
    ;
  printf("%-8s degree %2d: precompute %8.2f ms   trafo %8.2f ms   adjoint %8.2f ms\n",
      label, (int)p->spline_coeffs[0], 1e3 * tp, 1e3 * tt, 1e3 * ta);
}

int main(int argc, char **argv)
{
  nfft_plan p;
  int d = atoi(argv[1]), N1 = atoi(argv[2]), M = atoi(argv[3]), m = atoi(argv[5]);
  double sigma = atof(argv[4]);
  int reps = atoi(argv[7]), a, t, j, k;
  int N[3], n[3];
  unsigned flags = MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE | PRE_PHI_HUT;
  int rows, cols;

  if (!strcmp(argv[6], "pre"))
    flags |= PRE_PSI;
  for (t = 0; t < d; t++)
  {
    N[t] = N1;
    n[t] = (int)lrint(sigma * N1);
    if (n[t] % 2) n[t]++;
  }
  nfft_init_guru(&p, d, N, M, n, m, flags, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  srand(1);
  for (j = 0; j < M * d; j++)
    p.x[j] = (double)rand() / RAND_MAX - 0.5;
  for (k = 0; k < p.N_total; k++)
    p.f_hat[k] = ((double)rand() / RAND_MAX - 0.5) + I * ((double)rand() / RAND_MAX - 0.5);
  rows = 2 * (((int)p.spline_coeffs[0] + 1) / 2) + 1;
  cols = m + 2;
  printf("d=%d N=%d M=%d n=%d sigma=%g m=%d mode=%s\n", d, N1, M, n[0], sigma, m, argv[6]);
  timeit(&p, "stock", reps);

  for (a = 8; a < argc; a++)
  {
    FILE *f = fopen(argv[a], "r");
    int fdeg, fm, c, frows;
    double *one;
    if (!f || fscanf(f, "%d %d", &fdeg, &fm) != 2 || fm != m)
    {
      fprintf(stderr, "bad file %s\n", argv[a]);
      return 1;
    }
    frows = 2 * ((fdeg + 1) / 2) + 1;
    one = (double*) calloc((size_t)(frows * cols), sizeof(double));
    for (k = 0; k < frows; k++)
      for (c = 1; c <= m + 1; c++)
        if (fscanf(f, "%lf", &one[k * cols + c]) != 1)
        {
          fprintf(stderr, "short file %s\n", argv[a]);
          return 1;
        }
    fclose(f);
    free(p.spline_coeffs);
    p.spline_coeffs = (double*) calloc((size_t)(1 + d * frows * cols), sizeof(double));
    p.spline_coeffs[0] = (double)fdeg;
    for (t = 0; t < d; t++)
      memcpy(p.spline_coeffs + 1 + t * frows * cols, one, (size_t)(frows * cols) * sizeof(double));
    free(one);
    rows = frows;
    {
      const char *base = strrchr(argv[a], '/');
      char label[64];
      base = base ? base + 1 : argv[a];
      snprintf(label, sizeof label, "%.*s", (int)strcspn(base + 5, "_"), base + 5);
      timeit(&p, label, reps);
    }
  }
  (void)rows;
  nfft_finalize(&p);
  return 0;
}
