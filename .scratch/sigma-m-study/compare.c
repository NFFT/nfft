/* Tuned parameters against the legacy choice, head to head.
 *
 * DYADIC : (n, m) from Y(tune_plan_dyadic), which chooses among
 *          n = 2^j * next_power_of_2(N) for j in {0, 1, 2}.
 * LEGACY : n = 2 * next_power_of_2(N), the legacy default, with m found by
 *          searching upward for the smallest cut-off whose *measured* error
 *          meets the goal. The legacy API has no m estimator, so this gives
 *          it the best cut-off it could possibly have used -- an oracle.
 *          It is also rung 1 of the dyadic ladder.
 *
 * Both run through the same plan_ng path, so the only difference is the
 * parameter pair. Error is measured against a long-double direct NDFT, which
 * depends on (N, M, x) but not on n or m and is computed once per case.
 *
 * One CSV row per (case, goal, direction). Built once per precision.
 * argv[1] is the repetition count.
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
#define CRE crealf
#define CIM cimagf
#elif SWEEP_PREC == 1
typedef double R;
typedef double _Complex C;
#define NF(name) nfft_##name
#define PREC_NAME "double"
#define CRE creal
#define CIM cimag
#else
typedef long double R;
typedef long double _Complex C;
#define NF(name) nfftl_##name
#define PREC_NAME "long_double"
#define CRE creall
#define CIM cimagl
#endif

typedef long double LR;
typedef long double _Complex LC;
#define LPI 3.141592653589793238462643383279502884L
#define MAX_M 40

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

static NFFT_INT next_pow2(NFFT_INT x)
{
  NFFT_INT p = 1;
  while (p < x)
    p *= 2;
  return p;
}

/* Shared per-case state. */
static R *x_;
static C *in_fwd_, *in_adj_, *out_, *scratch_hat_, *scratch_f_;
static LC *ref_fwd_, *ref_adj_;
static LR nrm_fhat_, nrm_f_;
static NFFT_INT N_, M_;

/* Run one (n, m) and return the error; if secs is non-NULL, also time it. */
static LR run(NFFT_INT n, int m, int adjoint, double *secs, int reps)
{
  NFFT_INT Nl = N_, nl = n;
  NF(plan_ng) * p;
  const NFFT_INT out_len = adjoint ? N_ : M_;
  C *in = adjoint ? in_adj_ : in_fwd_;
  const LC *ref = adjoint ? ref_adj_ : ref_fwd_;
  const LR den = adjoint ? nrm_f_ : nrm_fhat_;
  LR worst = 0.0L;
  NFFT_INT j;

  p = NF(plan_ng_guru)(1, &Nl, NULL, &nl, M_, m, NFFT_WINDOW_KAISER_BESSEL, x_,
                       (void *)scratch_hat_, (void *)scratch_f_, FFTW_ESTIMATE,
                       NFFT_ESTIMATE | NFFT_NO_DIRECT);
  if (!p)
    return -1.0L;
  NF(precompute)(p);

  if (adjoint)
    NF(execute_adjoint_on)(p, (void *)out_, (void *)in);
  else
    NF(execute_on)(p, (void *)in, (void *)out_);

  for (j = 0; j < out_len; j++)
  {
    LC a = (LR)CRE(out_[j]) + I * (LR)CIM(out_[j]);
    LR e = cabsl(a - ref[j]);
    if (e > worst)
      worst = e;
  }

  if (secs)
  {
    int r;
    double t0;
    if (adjoint)
      NF(execute_adjoint_on)(p, (void *)out_, (void *)in);
    else
      NF(execute_on)(p, (void *)in, (void *)out_);
    t0 = now_seconds();
    for (r = 0; r < reps; r++)
    {
      if (adjoint)
        NF(execute_adjoint_on)(p, (void *)out_, (void *)in);
      else
        NF(execute_on)(p, (void *)in, (void *)out_);
    }
    *secs = (now_seconds() - t0) / (double)reps;
  }

  NF(plan_ng_destroy)(p);
  return den == 0.0L ? 0.0L : worst / den;
}

typedef struct
{
  NFFT_INT N;
  int shape; /* 0 = N-dominated, 1 = balanced, 2 = M-dominated */
} testcase;

static const char *shape_name(int s)
{
  return s == 0 ? "N-dominated" : (s == 1 ? "balanced" : "M-dominated");
}

int main(int argc, char **argv)
{
  /* Spread over t = next_power_of_2(N)/N, which is what decides whether the
   * dyadic ladder has a rung below the legacy grid at all. The older set was
   * {256, 251, 250, 243, 255, 512}, every one of which has t <= 1.06, so rung
   * 0 was illegal in all of it and the two tuners were indistinguishable.
   *
   *   N    256  512  255  251  250  243  200  320  160  600  1100  140
   *   t   1.00 1.00 1.00 1.02 1.02 1.05 1.28 1.60 1.60 1.71  1.86 1.83
   */
  static const NFFT_INT Ns[] = {256, 512, 255, 251, 250, 243,
                                200, 320, 160, 600, 1100, 140};
  static const double goals_f[] = {1e-2, 1e-4};
  static const double goals_d[] = {1e-4, 1e-8, 1e-11};
  static const double goals_l[] = {1e-8, 1e-14, 1e-20};
#if SWEEP_PREC == 0
  const double *goals = goals_f;
  const int n_goals = 2;
#elif SWEEP_PREC == 1
  const double *goals = goals_d;
  const int n_goals = 3;
#else
  const double *goals = goals_l;
  const int n_goals = 3;
#endif
  const int reps = (argc > 1) ? atoi(argv[1]) : 50;
  int iN, ish, ig, dir;

  printf("prec,N,M,shape,goal,dir,n_dya,m_dya,err_dya,t_dya,"
         "n_leg,m_leg,err_leg,t_leg,dya_met,leg_met\n");

  for (iN = 0; iN < (int)(sizeof(Ns) / sizeof(Ns[0])); iN++)
    for (ish = 0; ish < 3; ish++)
    {
      const NFFT_INT N = Ns[iN];
      const NFFT_INT M = ish == 0 ? (N / 4 > 1 ? N / 4 : 1)
                                  : (ish == 1 ? N : 4 * N);
      LR *xl;
      LC *fhat_l, *f_l;
      NFFT_INT j, k;

      N_ = N;
      M_ = M;
      xl = (LR *)malloc((size_t)M * sizeof(LR));
      fhat_l = (LC *)malloc((size_t)N * sizeof(LC));
      f_l = (LC *)malloc((size_t)M * sizeof(LC));
      ref_fwd_ = (LC *)malloc((size_t)M * sizeof(LC));
      ref_adj_ = (LC *)malloc((size_t)N * sizeof(LC));
      x_ = (R *)NF(malloc)((size_t)M * sizeof(R));
      in_fwd_ = (C *)NF(malloc)((size_t)N * sizeof(C));
      in_adj_ = (C *)NF(malloc)((size_t)M * sizeof(C));
      out_ = (C *)NF(malloc)((size_t)(N > M ? N : M) * sizeof(C));
      scratch_hat_ = (C *)NF(malloc)((size_t)N * sizeof(C));
      scratch_f_ = (C *)NF(malloc)((size_t)M * sizeof(C));
      memset(scratch_hat_, 0, (size_t)N * sizeof(C));
      memset(scratch_f_, 0, (size_t)M * sizeof(C));

      rng_seed(0x101d0000ULL + (unsigned long long)(N * 8 + ish));
      for (j = 0; j < M; j++)
      {
        xl[j] = (LR)(rng_uniform() - 0.5);
        x_[j] = (R)xl[j];
      }
      /* Real and imaginary parts on [0, 1), as gsweep.c draws them. Centred
       * data measures a forward error up to 2.6x smaller, so drawing it here
       * would time the tuner on an easier distribution than it was fitted
       * for. */
      nrm_fhat_ = 0.0L;
      for (k = 0; k < N; k++)
      {
        R re = (R)rng_uniform(), im = (R)rng_uniform();
        in_fwd_[k] = re + I * im;
        fhat_l[k] = (LR)re + I * (LR)im;
        nrm_fhat_ += cabsl(fhat_l[k]);
      }
      nrm_f_ = 0.0L;
      for (j = 0; j < M; j++)
      {
        R re = (R)rng_uniform(), im = (R)rng_uniform();
        in_adj_[j] = re + I * im;
        f_l[j] = (LR)re + I * (LR)im;
        nrm_f_ += cabsl(f_l[j]);
      }

      {
        const LR k0 = (LR)(-N / 2);
        for (j = 0; j < M; j++)
        {
          LC v = 0.0L;
          for (k = 0; k < N; k++)
            v += fhat_l[k] * twiddle(k0 + (LR)k, xl[j], -1);
          ref_fwd_[j] = v;
        }
        for (k = 0; k < N; k++)
        {
          LC v = 0.0L;
          for (j = 0; j < M; j++)
            v += f_l[j] * twiddle(k0 + (LR)k, xl[j], +1);
          ref_adj_[k] = v;
        }
      }

      for (ig = 0; ig < n_goals; ig++)
        for (dir = 0; dir <= 1; dir++)
        {
          const R goal = (R)goals[ig];
          NFFT_INT n_dya = 0, n_leg;
          int m_dya = 0, m_leg = 0;
          R att = 0;
          LR err_dya, err_leg = -1.0L;
          double t_dya = 0.0, t_leg = 0.0;
          int m;

          if (NF(tune_plan_dyadic)(N, M, dir, goal, &n_dya, &m_dya, &att) < 0)
            continue;
          err_dya = run(n_dya, m_dya, dir, &t_dya, reps);

          /* Legacy geometry, with an oracle search for its cut-off. */
          n_leg = 2 * next_pow2(N);
          for (m = 1; m <= MAX_M && m <= (int)(n_leg / 2 - 1); m++)
          {
            LR e = run(n_leg, m, dir, NULL, 0);
            if (e >= 0.0L && e <= (LR)goal)
            {
              m_leg = m;
              err_leg = e;
              break;
            }
          }
          if (m_leg == 0)
          {
            /* Never reached the goal: report its best cut-off. */
            LR best = -1.0L;
            int best_m = 1;
            for (m = 1; m <= MAX_M && m <= (int)(n_leg / 2 - 1); m++)
            {
              LR e = run(n_leg, m, dir, NULL, 0);
              if (e >= 0.0L && (best < 0.0L || e < best))
              {
                best = e;
                best_m = m;
              }
            }
            m_leg = best_m;
            err_leg = best;
          }
          run(n_leg, m_leg, dir, &t_leg, reps);

          printf("%s,%d,%d,%s,%.3e,%s,%d,%d,%.6Le,%.9f,%d,%d,%.6Le,%.9f,"
                 "%d,%d\n",
                 PREC_NAME, (int)N, (int)M, shape_name(ish), (double)goal,
                 dir ? "adjoint" : "forward",
                 (int)n_dya, m_dya, err_dya, t_dya,
                 (int)n_leg, m_leg, err_leg, t_leg,
                 err_dya <= (LR)goal ? 1 : 0, err_leg <= (LR)goal ? 1 : 0);
          fflush(stdout);
        }

      free(xl);
      free(fhat_l);
      free(f_l);
      free(ref_fwd_);
      free(ref_adj_);
      NF(free)(x_);
      NF(free)(in_fwd_);
      NF(free)(in_adj_);
      NF(free)(out_);
      NF(free)(scratch_hat_);
      NF(free)(scratch_f_);
    }

  return 0;
}
