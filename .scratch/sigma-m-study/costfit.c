/* Timing sweep for the cost weight in Y(tune_plan).
 *
 * The tuner ranks candidate pairs by
 *
 *     cost(n, m) = n*log2(n) + w * M * (2m+2)
 *
 * and only the ratio w matters. This driver measures one transform over a
 * grid of (N, n, m, M) that varies the two terms independently -- the
 * head-to-head data in compare-*.csv cannot separate them, since there n and
 * m move together.
 *
 * One CSV row per (N, n, m, M, direction). Built once per precision.
 * costfit.py does the fit.
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "nfft3.h"

#if SWEEP_PREC == 0
typedef float R;
typedef float _Complex C;
#define NF(name) nfftf_##name
#define PREC_NAME "float"
#elif SWEEP_PREC == 1
typedef double R;
typedef double _Complex C;
#define NF(name) nfft_##name
#define PREC_NAME "double"
#else
typedef long double R;
typedef long double _Complex C;
#define NF(name) nfftl_##name
#define PREC_NAME "long_double"
#endif

static unsigned long long rng_state = 88172645463325252ULL;
static double rng_uniform(void)
{
  rng_state ^= rng_state >> 12;
  rng_state ^= rng_state << 25;
  rng_state ^= rng_state >> 27;
  return (double)((rng_state * 2685821657736338717ULL) >> 11) /
         9007199254740992.0;
}

static double now_seconds(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Even 5-smooth sizes, so the FFT is the one the tuner would actually get. */
static int is_even_smooth5(NFFT_INT n)
{
  if (n < 2 || n % 2 != 0)
    return 0;
  while (n % 2 == 0)
    n /= 2;
  while (n % 3 == 0)
    n /= 3;
  while (n % 5 == 0)
    n /= 5;
  return n == 1;
}

int main(int argc, char **argv)
{
  static const NFFT_INT Ns[] = {243, 256, 500, 512};
  static const int ms[] = {1, 2, 3, 4, 6, 8, 12, 16};
  static const double m_factors[] = {0.125, 0.5, 2.0, 8.0, 32.0};
  const double budget = (argc > 1) ? atof(argv[1]) : 0.02;
  int iN, im, iM, dir;

  printf("prec,N,M,n,m,dir,reps,secs\n");

  for (iN = 0; iN < (int)(sizeof(Ns) / sizeof(Ns[0])); iN++)
    for (iM = 0; iM < (int)(sizeof(m_factors) / sizeof(m_factors[0])); iM++)
    {
      const NFFT_INT N = Ns[iN];
      NFFT_INT M = (NFFT_INT)(m_factors[iM] * (double)N);
      NFFT_INT n, j, k;
      R *x;
      C *in_fwd, *in_adj, *out, *sh, *sf;

      if (M < 1)
        M = 1;
      x = (R *)NF(malloc)((size_t)M * sizeof(R));
      in_fwd = (C *)NF(malloc)((size_t)N * sizeof(C));
      in_adj = (C *)NF(malloc)((size_t)M * sizeof(C));
      out = (C *)NF(malloc)((size_t)(N > M ? N : M) * sizeof(C));
      sh = (C *)NF(malloc)((size_t)N * sizeof(C));
      sf = (C *)NF(malloc)((size_t)M * sizeof(C));
      memset(sh, 0, (size_t)N * sizeof(C));
      memset(sf, 0, (size_t)M * sizeof(C));
      for (j = 0; j < M; j++)
      {
        x[j] = (R)(rng_uniform() - 0.5);
        in_adj[j] = (R)rng_uniform() + I * (R)rng_uniform();
      }
      for (k = 0; k < N; k++)
        in_fwd[k] = (R)rng_uniform() + I * (R)rng_uniform();

      /* Every even 5-smooth grid the tuner would consider, sigma in [5/4, 4]. */
      for (n = (5 * N + 3) / 4; n <= 4 * N; n++)
      {
        if (!is_even_smooth5(n))
          continue;
        for (im = 0; im < (int)(sizeof(ms) / sizeof(ms[0])); im++)
          for (dir = 0; dir <= 1; dir++)
          {
            const int m = ms[im];
            NFFT_INT Nl = N, nl = n;
            NF(plan_ng) * p;
            double t0, secs;
            int reps, r;

            if (m > (int)(n / 2 - 1))
              continue;
            p = NF(plan_ng_guru)(1, &Nl, NULL, &nl, M, m,
                                 NFFT_WINDOW_KAISER_BESSEL, x, (void *)sh,
                                 (void *)sf, FFTW_ESTIMATE,
                                 NFFT_ESTIMATE | NFFT_NO_DIRECT);
            if (!p)
              continue;
            NF(precompute)(p);

            /* Warm up, then time enough repetitions to fill the budget. */
            if (dir)
              NF(execute_adjoint_on)(p, (void *)out, (void *)in_adj);
            else
              NF(execute_on)(p, (void *)in_fwd, (void *)out);
            t0 = now_seconds();
            if (dir)
              NF(execute_adjoint_on)(p, (void *)out, (void *)in_adj);
            else
              NF(execute_on)(p, (void *)in_fwd, (void *)out);
            secs = now_seconds() - t0;
            reps = (int)(budget / (secs > 0.0 ? secs : 1e-9));
            if (reps < 3)
              reps = 3;
            if (reps > 20000)
              reps = 20000;

            t0 = now_seconds();
            for (r = 0; r < reps; r++)
            {
              if (dir)
                NF(execute_adjoint_on)(p, (void *)out, (void *)in_adj);
              else
                NF(execute_on)(p, (void *)in_fwd, (void *)out);
            }
            secs = (now_seconds() - t0) / (double)reps;
            NF(plan_ng_destroy)(p);

            printf("%s,%d,%d,%d,%d,%s,%d,%.12f\n", PREC_NAME, (int)N, (int)M,
                   (int)n, m, dir ? "adjoint" : "forward", reps, secs);
            fflush(stdout);
          }
      }

      NF(free)(x);
      NF(free)(in_fwd);
      NF(free)(in_adj);
      NF(free)(out);
      NF(free)(sh);
      NF(free)(sf);
    }

  return 0;
}
