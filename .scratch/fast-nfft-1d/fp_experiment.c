/* Standalone floating-point error experiment for the fast 1D NFFT.
 *
 * Goal: measure how the *rounding* error of the fast forward (trafo_1d) and
 * adjoint (adjoint_1d) transforms scales with the bandwidth N, and contrast it
 * with the direct (NDFT) transforms.
 *
 * A long-double reference is computed independently of the library so that the
 * measured error is dominated by the double-precision transform under test, not
 * by the reference.  The window approximation error is held essentially
 * constant by fixing m and the oversampling factor sigma=2, so the *trend* with
 * N reflects rounding behaviour.
 *
 * Build (from repo root, after `./configure --enable-all && make`):
 *   gcc -O2 -I include scratch_fp_experiment.c \
 *       kernel/.libs/libnfft3.a -lfftw3 -lm -o /tmp/fp_exp
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "nfft3.h"

static double frand(void) { return (double)rand() / (double)RAND_MAX - 0.5; }

/* reference forward NDFT in long double: f[j] = sum_k fhat[k] exp(-2pi i k x[j]) */
static void ref_trafo(int N, int M, const double *x, const double complex *fhat,
                      long double complex *f)
{
  for (int j = 0; j < M; j++)
  {
    long double complex v = 0.0L;
    long double xj = (long double)x[j];
    for (int kk = 0; kk < N; kk++)
    {
      long double k = (long double)(kk - N / 2);
      long double ph = -2.0L * 3.141592653589793238462643383279502884L * k * xj;
      /* reduce phase to [-pi,pi] to keep the reference itself accurate */
      long double tp = 2.0L * 3.141592653589793238462643383279502884L;
      ph = ph - tp * floorl(ph / tp + 0.5L);
      v += (long double complex)fhat[kk] * (cosl(ph) + I * sinl(ph));
    }
    f[j] = v;
  }
}

/* reference adjoint NDFT: fhat[k] = sum_j f[j] exp(+2pi i k x[j]) */
static void ref_adjoint(int N, int M, const double *x, const double complex *f,
                        long double complex *fhat)
{
  for (int kk = 0; kk < N; kk++)
  {
    long double complex v = 0.0L;
    long double k = (long double)(kk - N / 2);
    for (int j = 0; j < M; j++)
    {
      long double xj = (long double)x[j];
      long double ph = 2.0L * 3.141592653589793238462643383279502884L * k * xj;
      long double tp = 2.0L * 3.141592653589793238462643383279502884L;
      ph = ph - tp * floorl(ph / tp + 0.5L);
      v += (long double complex)f[j] * (cosl(ph) + I * sinl(ph));
    }
    fhat[kk] = v;
  }
}

static double rel_err_f(int M, const double complex *approx,
                        const long double complex *ref)
{
  long double num = 0.0L, den = 0.0L;
  for (int j = 0; j < M; j++)
  {
    long double complex d = (long double complex)approx[j] - ref[j];
    num += creall(d) * creall(d) + cimagl(d) * cimagl(d);
    den += creall(ref[j]) * creall(ref[j]) + cimagl(ref[j]) * cimagl(ref[j]);
  }
  return (double)(sqrtl(num) / sqrtl(den));
}

int main(void)
{
  const int m = 12;           /* high cutoff: approx error << rounding floor  */
  printf("# m=%d sigma=2\n", m);
  printf("# %8s  %14s  %14s  %14s  %14s\n", "N", "fast_fwd", "direct_fwd",
         "fast_adj", "direct_adj");

  for (int e = 4; e <= 13; e++)
  {
    int N = 1 << e;
    int M = N;                /* M = N nodes */
    nfft_plan p;
    nfft_init_1d(&p, N, M);   /* default flags: PRE_PHI_HUT|PRE_PSI, m from window */
    /* override m via re-init using guru so we control the approximation error */
    nfft_finalize(&p);
    int Narr[1] = { N };
    int narr[1] = { 2 * (int)nfft_next_power_of_2(N) };
    unsigned flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
                     FFTW_INIT | FFT_OUT_OF_PLACE;
    nfft_init_guru(&p, 1, Narr, M, narr, m, flags,
                   FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    srand(1234 + e);
    for (int j = 0; j < M; j++)
      p.x[j] = frand();
    if (p.flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);

    /* ---- forward ---- */
    double complex *fhat = malloc(sizeof(double complex) * N);
    for (int k = 0; k < N; k++)
    {
      p.f_hat[k] = frand() + I * frand();
      fhat[k] = p.f_hat[k];
    }
    long double complex *fref = malloc(sizeof(long double complex) * M);
    ref_trafo(N, M, p.x, fhat, fref);

    nfft_trafo_1d(&p);
    double err_fast_fwd = rel_err_f(M, p.f, fref);
    nfft_trafo_direct(&p);
    double err_direct_fwd = rel_err_f(M, p.f, fref);

    /* ---- adjoint ---- */
    double complex *fin = malloc(sizeof(double complex) * M);
    for (int j = 0; j < M; j++)
    {
      p.f[j] = frand() + I * frand();
      fin[j] = p.f[j];
    }
    long double complex *fhatref = malloc(sizeof(long double complex) * N);
    ref_adjoint(N, M, p.x, fin, fhatref);

    nfft_adjoint_1d(&p);
    double err_fast_adj = rel_err_f(N, p.f_hat, fhatref);
    nfft_adjoint_direct(&p);
    double err_direct_adj = rel_err_f(N, p.f_hat, fhatref);

    printf("  %8d  %14.4e  %14.4e  %14.4e  %14.4e\n", N, err_fast_fwd,
           err_direct_fwd, err_fast_adj, err_direct_adj);

    free(fhat); free(fref); free(fin); free(fhatref);
    nfft_finalize(&p);
  }
  return 0;
}
