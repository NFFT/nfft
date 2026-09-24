#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "nfft3.h"
static double run(int n, int N, int m, unsigned extra)
{
  const int M = 3000;
  nfftl_plan p;
  nfftl_init_guru(&p, 1, (int[]){N}, M, (int[]){n}, m,
      PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT
      | FFT_OUT_OF_PLACE | extra, FFTW_ESTIMATE);
  srand(1);
  for (int j = 0; j < M; j++) p.x[j] = rand() / (RAND_MAX + 1.0) - 0.5;
  for (int k = 0; k < N; k++)
    p.f_hat[k] = (rand() / (RAND_MAX + 1.0) - 0.5) + I * (rand() / (RAND_MAX + 1.0) - 0.5);
  nfftl_precompute_one_psi(&p);
  nfftl_trafo_direct(&p);
  long double complex *d = malloc(M * sizeof *d);
  for (int j = 0; j < M; j++) d[j] = p.f[j];
  nfftl_trafo(&p);
  long double e = 0, r = 0;
  for (int j = 0; j < M; j++) { e += powl(cabsl(p.f[j] - d[j]), 2); r += powl(cabsl(d[j]), 2); }
  free(d); nfftl_finalize(&p);
  return (double)sqrtl(e / r);
}
int main(void)
{
  const int N = 200; double sig[] = {2.0, 4.0};
  for (int s = 0; s < 2; s++)
    for (int m = 6; m <= 14; m++)
    {
      int n = 2 * (int)lround(sig[s] * N / 2);
      double a = run(n, N, m, ANALYTIC_WINDOW), q = run(n, N, m, 0);
      printf("sigma=%.2f m=%d analytic=%.2e poly=%.2e ratio=%.2f\n", sig[s], m, a, q, q / a);
    }
  return 0;
}
