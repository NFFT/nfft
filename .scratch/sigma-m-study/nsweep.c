/* Fine-grained n sweep: is it worth pushing sigma past 2 to reach a size the
 * FFT likes?
 *
 * For each bandwidth N, every integer n from 2N to 4N is measured for both
 * accuracy and transform time, at a fixed cut-off m so only n varies. The
 * long-double NDFT reference depends on (N, M, x, f_hat) but not on n, so it
 * is computed once per N and reused across the whole row.
 *
 * Writes one CSV row per (N, n). Double precision only.
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "nfft3.h"

typedef double R;
typedef double _Complex C;
typedef long double LR;
typedef long double _Complex LC;

#define NF(name) nfft_##name
#define LPI 3.141592653589793238462643383279502884L

static unsigned long long rng_state;
static void rng_seed(unsigned long long s)
{
  rng_state = s ? s : 88172645463325252ULL;
}
static double rng_uniform(void)
{
  rng_state ^= rng_state >> 12;
  rng_state ^= rng_state << 25;
  rng_state ^= rng_state >> 27;
  return (double)((rng_state * 2685821657736338717ULL) >> 11) /
         9007199254740992.0;
}

static LC twiddle(LR k, LR x, int sign)
{
  LR t = fmodl(k * x, 1.0L);
  LR phi = 2.0L * LPI * t;
  return cosl(phi) + (LR)sign * I * sinl(phi);
}

static double now_seconds(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Largest prime factor of v, and whether v is a power of two. */
static int largest_prime_factor(int v)
{
  int f = 2, best = 1;
  while ((long)f * f <= (long)v)
  {
    while (v % f == 0)
    {
      best = f;
      v /= f;
    }
    f++;
  }
  if (v > 1)
    best = v;
  return best;
}

int main(int argc, char **argv)
{
  static const int Ns[] = {512, 509, 500, 486, 517, 1000};
  const int n_N = (int)(sizeof(Ns) / sizeof(Ns[0]));
  const int m = (argc > 1) ? atoi(argv[1]) : 6;
  const int reps = (argc > 2) ? atoi(argv[2]) : 20;
  int iN, nn, j, k, rep;

  printf("N,n,sigma,m,lpf,is_pow2,err_fwd,err_adj,t_fwd,t_adj,t_pre\n");

  for (iN = 0; iN < n_N; iN++)
  {
    const int N = Ns[iN];
    const int M = 2 * N;
    LR *xl;
    LC *fhat_l, *f_l, *ref_f, *ref_fhat;
    R *x;
    C *fhat_in, *f_in, *f_out, *fhat_out;
    LR nrm1_fhat = 0.0L, nrm1_f = 0.0L;

    xl = (LR *)malloc((size_t)M * sizeof(LR));
    fhat_l = (LC *)malloc((size_t)N * sizeof(LC));
    f_l = (LC *)malloc((size_t)M * sizeof(LC));
    ref_f = (LC *)malloc((size_t)M * sizeof(LC));
    ref_fhat = (LC *)malloc((size_t)N * sizeof(LC));
    x = (R *)NF(malloc)((size_t)M * sizeof(R));
    fhat_in = (C *)NF(malloc)((size_t)N * sizeof(C));
    f_in = (C *)NF(malloc)((size_t)M * sizeof(C));
    f_out = (C *)NF(malloc)((size_t)M * sizeof(C));
    fhat_out = (C *)NF(malloc)((size_t)N * sizeof(C));

    rng_seed(0xabcd0000ULL + (unsigned long long)N);
    for (j = 0; j < M; j++)
    {
      xl[j] = (LR)(rng_uniform() - 0.5);
      x[j] = (R)xl[j];
    }
    for (k = 0; k < N; k++)
    {
      R re = (R)(rng_uniform() - 0.5), im = (R)(rng_uniform() - 0.5);
      fhat_in[k] = re + I * im;
      fhat_l[k] = (LR)re + I * (LR)im;
      nrm1_fhat += cabsl(fhat_l[k]);
    }
    for (j = 0; j < M; j++)
    {
      R re = (R)(rng_uniform() - 0.5), im = (R)(rng_uniform() - 0.5);
      f_in[j] = re + I * im;
      f_l[j] = (LR)re + I * (LR)im;
      nrm1_f += cabsl(f_l[j]);
    }

    /* Reference, once per N: independent of n. */
    {
      const LR k0 = (LR)(-N / 2);
      for (j = 0; j < M; j++)
      {
        LC v = 0.0L;
        for (k = 0; k < N; k++)
          v += fhat_l[k] * twiddle(k0 + (LR)k, xl[j], -1);
        ref_f[j] = v;
      }
      for (k = 0; k < N; k++)
      {
        LC v = 0.0L;
        for (j = 0; j < M; j++)
          v += f_l[j] * twiddle(k0 + (LR)k, xl[j], +1);
        ref_fhat[k] = v;
      }
    }

    for (nn = 2 * N; nn <= 4 * N; nn++)
    {
      NFFT_INT Nl = (NFFT_INT)N, nl = (NFFT_INT)nn;
      NF(plan_ng) * p;
      LR e_fwd = 0.0L, e_adj = 0.0L;
      double t0, t_pre, t_fwd, t_adj;

      t0 = now_seconds();
      p = NF(plan_ng_guru)(1, &Nl, NULL, &nl, (NFFT_INT)M, m,
                           NFFT_WINDOW_KAISER_BESSEL, x, (void *)fhat_out,
                           (void *)f_out, FFTW_ESTIMATE,
                           NFFT_ESTIMATE | NFFT_NO_DIRECT);
      if (!p)
      {
        fprintf(stderr, "guru declined N=%d n=%d\n", N, nn);
        continue;
      }
      NF(precompute)(p);
      t_pre = now_seconds() - t0;

      /* Warm up, then time the mean over reps. */
      NF(execute_on)(p, (void *)fhat_in, (void *)f_out);
      t0 = now_seconds();
      for (rep = 0; rep < reps; rep++)
        NF(execute_on)(p, (void *)fhat_in, (void *)f_out);
      t_fwd = (now_seconds() - t0) / (double)reps;

      NF(execute_adjoint_on)(p, (void *)fhat_out, (void *)f_in);
      t0 = now_seconds();
      for (rep = 0; rep < reps; rep++)
        NF(execute_adjoint_on)(p, (void *)fhat_out, (void *)f_in);
      t_adj = (now_seconds() - t0) / (double)reps;

      NF(execute_on)(p, (void *)fhat_in, (void *)f_out);
      for (j = 0; j < M; j++)
      {
        LC a = (LR)creal(f_out[j]) + I * (LR)cimag(f_out[j]);
        LR e = cabsl(a - ref_f[j]);
        if (e > e_fwd)
          e_fwd = e;
      }
      NF(execute_adjoint_on)(p, (void *)fhat_out, (void *)f_in);
      for (k = 0; k < N; k++)
      {
        LC a = (LR)creal(fhat_out[k]) + I * (LR)cimag(fhat_out[k]);
        LR e = cabsl(a - ref_fhat[k]);
        if (e > e_adj)
          e_adj = e;
      }

      NF(plan_ng_destroy)(p);

      printf("%d,%d,%.10f,%d,%d,%d,%.10Le,%.10Le,%.9f,%.9f,%.9f\n", N, nn,
             (double)nn / (double)N, m, largest_prime_factor(nn),
             (nn & (nn - 1)) == 0 ? 1 : 0, e_fwd / nrm1_fhat, e_adj / nrm1_f,
             t_fwd, t_adj, t_pre);
    }
    fflush(stdout);

    free(xl);
    free(fhat_l);
    free(f_l);
    free(ref_f);
    free(ref_fhat);
    NF(free)(x);
    NF(free)(fhat_in);
    NF(free)(f_in);
    NF(free)(f_out);
    NF(free)(fhat_out);
  }

  return 0;
}
