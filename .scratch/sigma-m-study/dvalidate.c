/* Does the dyadic tuner's cut-off hold on data the fit never saw?
 *
 * The envelope dominates a sweep that records the worst of five trials per
 * cell. The tuner's promise is about one draw, which is a weaker statement,
 * and the gap between the two is what sizes the safety margin. So: take the
 * pair Y(tune_plan_dyadic) returns, measure it against a direct NDFT over
 * many fresh draws, and count how many exceed the goal.
 *
 * One CSV row per (N, M, goal, direction), argv[1] draws per row.
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nfft3.h"
#include "infft.h"

static R measure(INT N, INT n, INT M, int m, int adjoint, int seed)
{
  R *x;
  C *in, *out, *ref, *fh, *ff;
  Y(plan_ng) * fast;
  Y(plan_ng) * direct;
  const INT in_len = adjoint ? M : N;
  const INT out_len = adjoint ? N : M;
  R num = K(0.0), den = K(0.0);
  INT j;

  x = (R *)Y(malloc)((size_t)M * sizeof(R));
  in = (C *)Y(malloc)((size_t)in_len * sizeof(C));
  out = (C *)Y(malloc)((size_t)out_len * sizeof(C));
  ref = (C *)Y(malloc)((size_t)out_len * sizeof(C));
  fh = (C *)Y(malloc)((size_t)N * sizeof(C));
  ff = (C *)Y(malloc)((size_t)M * sizeof(C));
  memset(fh, 0, (size_t)N * sizeof(C));
  memset(ff, 0, (size_t)M * sizeof(C));

  Y(srand48)(seed);
  Y(vrand_shifted_unit_double)(x, M);
  Y(vrand_unit_complex)(in, in_len);

  direct = Y(plan_ng_guru)(1, &N, 0, &n, M, m, NFFT_WINDOW_KAISER_BESSEL, x,
                           (void *)fh, (void *)ff, 0u,
                           NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
  if (!direct)
    return K(-1.0);
  Y(precompute)(direct);
  if (adjoint)
    Y(execute_adjoint_on)(direct, (void *)ref, (void *)in);
  else
    Y(execute_on)(direct, (void *)in, (void *)ref);
  Y(plan_ng_destroy)(direct);

  fast = Y(plan_ng_guru)(1, &N, 0, &n, M, m, NFFT_WINDOW_KAISER_BESSEL, x,
                         (void *)fh, (void *)ff, 0u,
                         NFFT_ESTIMATE | NFFT_NO_DIRECT);
  if (!fast)
    return K(-1.0);
  Y(precompute)(fast);
  if (adjoint)
    Y(execute_adjoint_on)(fast, (void *)out, (void *)in);
  else
    Y(execute_on)(fast, (void *)in, (void *)out);
  Y(plan_ng_destroy)(fast);

  for (j = 0; j < out_len; j++) {
    const R e = CABS(out[j] - ref[j]);
    if (e > num)
      num = e;
  }
  for (j = 0; j < in_len; j++)
    den += CABS(in[j]);

  Y(free)(x);
  Y(free)(in);
  Y(free)(out);
  Y(free)(ref);
  Y(free)(fh);
  Y(free)(ff);
  return den > K(0.0) ? num / den : K(-1.0);
}

int main(int argc, char **argv)
{
  /* N spanning t = next_power_of_2(N)/N over (1, 2]. */
  static const INT Ns[] = {64, 100, 130, 140, 160, 200, 251, 256,
                           300, 320, 500, 600};
  static const double goals[] = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
  static const int mfac[] = {-4, 1, 4}; /* M = N/4, N, 4N */
  const int draws = (argc > 1) ? atoi(argv[1]) : 24;
  size_t iN, ig, iM;
  int dir, s;

  printf("N,M,goal,dir,n,m,rc,worst,worst_over_goal,misses,draws\n");

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (iM = 0; iM < sizeof(mfac) / sizeof(mfac[0]); iM++)
      for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++)
        for (dir = 0; dir <= 1; dir++) {
          const INT N = Ns[iN];
          const INT M = mfac[iM] < 0 ? (N / (-mfac[iM]) > 1 ? N / (-mfac[iM])
                                                            : 1)
                                     : N * mfac[iM];
          const R goal = (R)goals[ig];
          INT n = 0;
          int m = 0, misses = 0;
          R att = K(0.0), worst = K(0.0);
          const int rc =
               Y(tune_plan_dyadic)(N, M, dir, goal, &n, &m, &att);

          if (rc != 1)
            continue; /* capped goals are a different question */

          for (s = 0; s < draws; s++) {
            const R e = measure(N, n, M, m, dir, 5000 + s);
            if (e < K(0.0))
              continue;
            if (e > worst)
              worst = e;
            if (e > goal)
              misses++;
          }

          printf("%d,%d,%.3e,%s,%d,%d,%d,%.6e,%.4f,%d,%d\n", (int)N, (int)M,
                 (double)goal, dir ? "adjoint" : "forward", (int)n, m, rc,
                 (double)worst, (double)(worst / goal), misses, draws);
          fflush(stdout);
        }

  return 0;
}
