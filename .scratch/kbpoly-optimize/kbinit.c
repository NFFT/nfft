/* Plan-time cost of OPT_WINDOW: nfft_init_guru with and without the flag. */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "nfft3.h"
static double now(void){struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec;}
int main(int argc, char **argv)
{
  int d = atoi(argv[1]), N1 = atoi(argv[2]), m = atoi(argv[4]), t, r, o;
  double sigma = atof(argv[3]);
  int N[3], n[3];
  for (t = 0; t < d; t++) { N[t] = N1; n[t] = (int)(sigma * N1); n[t] += n[t] % 2; }
  for (o = 0; o < 2; o++)
  {
    double best = 1e9;
    for (r = 0; r < 5; r++)
    {
      nfft_plan p; double t0 = now();
      nfft_init_guru(&p, d, N, 16, n, m, PRE_PHI_HUT | MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT | (o ? OPT_WINDOW : 0), FFTW_ESTIMATE);
      double t1 = now(); nfft_finalize(&p);
      if (t1 - t0 < best) best = t1 - t0;
    }
    printf("d=%d N=%d sigma=%g m=%d %s: init %.2f ms (degree %d)\n", d, N1, sigma, m, o ? "OPT_WINDOW" : "plain     ", 1e3 * best, 0);
  }
  return 0;
}
