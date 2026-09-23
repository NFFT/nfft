/* Time one 1-D NFFT trafo and adjoint, and the FFT inside them, per (N, M, n, m).
 *
 * Flags as nfft_init_1d: PRE_PHI_HUT | PRE_PSI, FFTW_ESTIMATE. Precomputation
 * is outside the timed region: the cost is that of a transform on nodes that
 * are already set up. The FFT is timed through the plan's own FFTW plan, so
 * the rest of the trafo -- deconvolution plus convolution -- is the difference.
 *
 * usage: costprobe < configs      one "N M n m" per line
 * rows:  N,M,n,m,trafo,adjoint,fft   seconds, minimum over repeats
 */
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <nfft3.h>

static double now(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

/* Minimum over 5 rounds of the mean over enough calls to fill 20 ms. */
static double timeit(void (*f)(void *), void *arg)
{
  double best = 1e30;
  int reps = 1, r, i;

  for (;;)
  {
    double t0 = now();
    for (i = 0; i < reps; i++)
      f(arg);
    if (now() - t0 > 0.02)
      break;
    reps *= 2;
  }
  for (r = 0; r < 5; r++)
  {
    double t0 = now(), t;
    for (i = 0; i < reps; i++)
      f(arg);
    t = (now() - t0) / reps;
    if (t < best)
      best = t;
  }
  return best;
}

static void run_trafo(void *p) { nfft_trafo((nfft_plan *)p); }
static void run_adjoint(void *p) { nfft_adjoint((nfft_plan *)p); }
static void run_fft(void *p) { fftw_execute(((nfft_plan *)p)->my_fftw_plan1); }

int main(void)
{
  int N, M, n, m, j;

  srand(1);
  printf("N,M,n,m,trafo,adjoint,fft\n");
  while (scanf("%d %d %d %d", &N, &M, &n, &m) == 4)
  {
    nfft_plan p;
    double tt, ta, tf;

    nfft_init_guru(&p, 1, &N, M, &n, m,
                   PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
                       FFTW_INIT | FFT_OUT_OF_PLACE,
                   FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    for (j = 0; j < M; j++)
      p.x[j] = (double)rand() / RAND_MAX - 0.5;
    for (j = 0; j < N; j++)
      p.f_hat[j] = (double)rand() / RAND_MAX + I * ((double)rand() / RAND_MAX);
    nfft_precompute_one_psi(&p);

    tt = timeit(run_trafo, &p);
    ta = timeit(run_adjoint, &p);
    tf = timeit(run_fft, &p);
    printf("%d,%d,%d,%d,%.4e,%.4e,%.4e\n", N, M, n, m, tt, ta, tf);
    fflush(stdout);
    nfft_finalize(&p);
  }
  return 0;
}
