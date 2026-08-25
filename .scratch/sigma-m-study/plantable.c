/* Show what Y(tune_plan) picks over a grid of bandwidths and goals. */

#include <stdio.h>

#include "nfft3.h"

#if SWEEP_PREC == 0
typedef float R;
#define NF(name) nfftf_##name
#define PREC_NAME "float"
#elif SWEEP_PREC == 1
typedef double R;
#define NF(name) nfft_##name
#define PREC_NAME "double"
#else
typedef long double R;
#define NF(name) nfftl_##name
#define PREC_NAME "long double"
#endif

/* The pair depends on the node count; this table shows one shape. */
#define MNODES(N) (N)
static const double GOALS[] = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-16};
#define NG ((int)(sizeof(GOALS) / sizeof(GOALS[0])))

int main(void)
{
  static const NFFT_INT Ns[] = {128, 500, 509, 512, 517, 1000};
  int dir, iN, ig;

  for (dir = 0; dir <= 1; dir++)
  {
    printf("\n%s, %s -- tune_plan()\n", PREC_NAME,
           dir ? "adjoint" : "forward");
    printf("%-7s", "N");
    for (ig = 0; ig < NG; ig++)
      printf(" %14.0e", GOALS[ig]);
    printf("\n");
    for (iN = 0; iN < (int)(sizeof(Ns) / sizeof(Ns[0])); iN++)
    {
      printf("%-7d", (int)Ns[iN]);
      for (ig = 0; ig < NG; ig++)
      {
        NFFT_INT n = 0;
        int m = 0;
        R att = 0;
        int rc = NF(tune_plan)(Ns[iN], MNODES(Ns[iN]), dir, (R)GOALS[ig], &n, &m,
                               &att);
        if (rc < 0)
          printf(" %14s", "-");
        else
          printf(" %7d/m=%-2d%s", (int)n, m, rc == 1 ? " " : "*");
      }
      printf("\n");
    }
  }

  printf("\nreachable floor (goal capped there), %s:\n", PREC_NAME);
  for (dir = 0; dir <= 1; dir++)
    for (iN = 0; iN < (int)(sizeof(Ns) / sizeof(Ns[0])); iN++)
    {
      NFFT_INT n = 0;
      int m = 0;
      R att = 0;
      NF(tune_plan)(Ns[iN], MNODES(Ns[iN]), dir, (R)1e-30, &n, &m, &att);
      printf("  %-8s N=%-5d -> n=%-6d sigma=%.3f m=%-3d attained=%.3e\n",
             dir ? "adjoint" : "forward", (int)Ns[iN], (int)n,
             (double)n / (double)Ns[iN], m, (double)att);
    }
  printf("\n* goal below the reachable floor; capped\n");
  return 0;
}
