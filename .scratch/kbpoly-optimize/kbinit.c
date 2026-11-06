/* Plan-time cost of OPT_WINDOW: nfft_init_guru with and without the flag,
 * median and minimum of reps (default 7). */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "nfft3.h"
static double now(void){struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec;}
static int cmp(const void *a, const void *b){double x = *(const double*)a, y = *(const double*)b; return (x > y) - (x < y);}
int main(int argc, char **argv)
{
  int d = atoi(argv[1]), N1 = atoi(argv[2]), m = atoi(argv[4]), reps = argc > 5 ? atoi(argv[5]) : 7, t, r, o;
  double sigma = atof(argv[3]), med[2], *ts = malloc(reps * sizeof(double));
  int N[3], n[3];
  for (t = 0; t < d; t++) { N[t] = N1; n[t] = (int)(sigma * N1); n[t] += n[t] % 2; }
  for (o = 0; o < 2; o++)
  {
    for (r = 0; r < reps; r++)
    {
      nfft_plan p; double t0 = now();
      nfft_init_guru(&p, d, N, 16, n, m, PRE_PHI_HUT | MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT | (o ? OPT_WINDOW : 0), FFTW_ESTIMATE);
      ts[r] = 1e3 * (now() - t0); nfft_finalize(&p);
    }
    qsort(ts, reps, sizeof(double), cmp);
    med[o] = reps % 2 ? ts[reps / 2] : 0.5 * (ts[reps / 2 - 1] + ts[reps / 2]);
    printf("d=%d N=%d sigma=%g m=%d %s: init median %.2f min %.2f ms\n", d, N1, sigma, m, o ? "OPT_WINDOW" : "plain     ", med[o], ts[0]);
  }
  printf("d=%d N=%d sigma=%g m=%d OPT_WINDOW - plain: %.2f ms\n", d, N1, sigma, m, med[1] - med[0]);
  free(ts);
  return 0;
}
