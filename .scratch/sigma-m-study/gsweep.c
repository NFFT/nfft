/* Kaiser-Bessel accuracy sweep over the full geometry (sigma, N, M, m).
 *
 * sweep.c ties M to 2N, so N and M cannot be separated in the fit and the
 * model in kernel/nfft/tune.c carries neither. This one crosses N with M and
 * measures the same two errors as sweep.c and tests/nfft.c,
 *
 *   forward: max_j |f_j - s_j| / sum_k |f_hat_k|
 *   adjoint: max_k |fhat_k - s_k| / sum_j |f_j|
 *
 * against a long-double direct NDFT. The reference depends on (N, M, x, data)
 * but not on sigma or m, so it is computed once per (N, M, trial) and reused
 * across the whole sigma/m plane -- the reason this grid is affordable at all.
 *
 * One CSV row per (N, M, sigma, m, trial), same columns as sweep.c.
 * Build one binary per precision with -DSWEEP_PREC={0,1,2}.
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nfft3.h"

#if SWEEP_PREC == 0
typedef float R;
typedef float _Complex C;
#define NF(name) nfftf_##name
#define PREC_NAME "float"
#define CRE crealf
#define CIM cimagf
#elif SWEEP_PREC == 1
typedef double R;
typedef double _Complex C;
#define NF(name) nfft_##name
#define PREC_NAME "double"
#define CRE creal
#define CIM cimag
#elif SWEEP_PREC == 2
typedef long double R;
typedef long double _Complex C;
#define NF(name) nfftl_##name
#define PREC_NAME "long_double"
#define CRE creall
#define CIM cimagl
#else
#error set SWEEP_PREC
#endif

typedef long double LR;
typedef long double _Complex LC;

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

int main(int argc, char **argv)
{
  static const double sigmas[] = {1.25, 1.3,  1.35, 1.4, 1.5, 1.6, 1.75,
                                  2.0,  2.25, 2.5,  3.0, 3.5, 4.0};
  int Ns[16] = {32, 64, 128, 256, 512, 1024};
  /* M as a fraction of N, so the two are crossed rather than tied. */
  double m_factors[16] = {0.25, 1.0, 2.0, 8.0};
  const int n_sigma = (int)(sizeof(sigmas) / sizeof(sigmas[0]));
  int n_N = 6, n_M = 4;
  const int m_max_cap = (argc > 1) ? atoi(argv[1]) : 32;
  const int trials = (argc > 2) ? atoi(argv[2]) : 5;

  /* Optional comma-separated overrides, so the same binary can be pointed at
   * a validation grid outside the fitted box. */
  if (argc > 3)
  {
    char *t = strtok(argv[3], ",");
    for (n_N = 0; t && n_N < 16; t = strtok(0, ","))
      Ns[n_N++] = atoi(t);
  }
  if (argc > 4)
  {
    char *t = strtok(argv[4], ",");
    for (n_M = 0; t && n_M < 16; t = strtok(0, ","))
      m_factors[n_M++] = atof(t);
  }

  int is, iN, iM, m, t, j, k;

  printf("prec,sigma_req,sigma,N,n,M,m,trial,err_fwd,err_adj\n");

  for (iN = 0; iN < n_N; iN++)
    for (iM = 0; iM < n_M; iM++)
    {
      const int N = Ns[iN];
      int M = (int)(m_factors[iM] * (double)N + 0.5);
      LR *xl;
      LC *fhat_l, *f_l, *ref_f, *ref_fhat;
      R *x;
      C *fhat_in, *f_in, *f_out, *fhat_out;
      LR nrm1_fhat, nrm1_f;

      if (M < 1)
        M = 1;

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

      for (t = 0; t < trials; t++)
      {
        rng_seed(0x9eed0000ULL +
                 (unsigned long long)(iN * 10000 + iM * 100 + t));

        for (j = 0; j < M; j++)
        {
          xl[j] = (LR)(rng_uniform() - 0.5);
          x[j] = (R)xl[j];
        }
        nrm1_fhat = 0.0L;
        for (k = 0; k < N; k++)
        {
          R re = (R)(rng_uniform() - 0.5), im = (R)(rng_uniform() - 0.5);
          fhat_in[k] = re + I * im;
          fhat_l[k] = (LR)re + I * (LR)im;
          nrm1_fhat += cabsl(fhat_l[k]);
        }
        nrm1_f = 0.0L;
        for (j = 0; j < M; j++)
        {
          R re = (R)(rng_uniform() - 0.5), im = (R)(rng_uniform() - 0.5);
          f_in[j] = re + I * im;
          f_l[j] = (LR)re + I * (LR)im;
          nrm1_f += cabsl(f_l[j]);
        }

        /* Hoisted out of the sigma/m plane below. */
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

        for (is = 0; is < n_sigma; is++)
        {
          int n = (int)(sigmas[is] * N + 0.5);
          NFFT_INT Nl = (NFFT_INT)N, nl;
          int m_max;

          if (n % 2)
            n++;
          if (n <= N)
            n = N + 2;
          nl = (NFFT_INT)n;
          m_max = n / 2 - 1;
          if (m_max > m_max_cap)
            m_max = m_max_cap;

          for (m = 1; m <= m_max; m++)
          {
            NF(plan_ng) * p;
            LR worst_f = 0.0L, worst_fhat = 0.0L;

            p = NF(plan_ng_guru)(1, &Nl, NULL, &nl, (NFFT_INT)M, m,
                                 NFFT_WINDOW_KAISER_BESSEL, x,
                                 (void *)fhat_out, (void *)f_out, 0u,
                                 NFFT_ESTIMATE | NFFT_NO_DIRECT);
            if (!p)
              continue;
            NF(precompute)(p);

            NF(execute_on)(p, (void *)fhat_in, (void *)f_out);
            for (j = 0; j < M; j++)
            {
              LC a = (LR)CRE(f_out[j]) + I * (LR)CIM(f_out[j]);
              LR e = cabsl(a - ref_f[j]);
              if (e > worst_f)
                worst_f = e;
            }

            NF(execute_adjoint_on)(p, (void *)fhat_out, (void *)f_in);
            for (k = 0; k < N; k++)
            {
              LC a = (LR)CRE(fhat_out[k]) + I * (LR)CIM(fhat_out[k]);
              LR e = cabsl(a - ref_fhat[k]);
              if (e > worst_fhat)
                worst_fhat = e;
            }

            NF(plan_ng_destroy)(p);

            printf("%s,%.6f,%.10f,%d,%d,%d,%d,%d,%.10Le,%.10Le\n", PREC_NAME,
                   sigmas[is], (double)n / (double)N, N, n, M, m, t,
                   nrm1_fhat > 0.0L ? worst_f / nrm1_fhat : 0.0L,
                   nrm1_f > 0.0L ? worst_fhat / nrm1_f : 0.0L);
          }
        }
        fflush(stdout);
      }

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
