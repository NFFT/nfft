/* Print what Y(tune) and Y(tune_sigma) return over a grid, for the report. */

#include <stdio.h>
#include <stdlib.h>

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

static const double GOALS[] = {1e-2,  1e-4,  1e-6,  1e-8,
                               1e-10, 1e-12, 1e-14, 1e-16};
#define NGOALS ((int)(sizeof(GOALS) / sizeof(GOALS[0])))

int main(void)
{
  static const double sigmas[] = {1.25, 1.5, 2.0, 3.0, 4.0};
  const NFFT_INT N = 1024;
  int dir, is, ig;

  for (dir = 0; dir <= 1; dir++)
  {
    printf("\n%s, %s -- m from tune(), N=%d\n", PREC_NAME,
           dir ? "adjoint" : "forward", (int)N);
    printf("%-7s", "sigma");
    for (ig = 0; ig < NGOALS; ig++)
      printf(" %8.0e", GOALS[ig]);
    printf("\n");
    for (is = 0; is < (int)(sizeof(sigmas) / sizeof(sigmas[0])); is++)
    {
      NFFT_INT n = (NFFT_INT)(sigmas[is] * (double)N + 0.5);
      if (n % 2)
        n++;
      printf("%-7.2f", (double)n / (double)N);
      for (ig = 0; ig < NGOALS; ig++)
      {
        int m = 0;
        R att = 0;
        int rc = NF(tune)(N, n, 2 * N, dir, (R)GOALS[ig], &m, &att);
        if (rc == 1)
          printf(" %8d", m);
        else if (rc == 0)
          printf(" %7d*", m);
        else
          printf("        -");
      }
      printf("\n");
    }
  }

  printf("\n%s -- sigma from tune_sigma(), N=%d\n", PREC_NAME, (int)N);
  printf("%-10s", "goal");
  printf(" %14s %14s\n", "forward", "adjoint");
  for (ig = 0; ig < NGOALS; ig++)
  {
    printf("%-10.0e", GOALS[ig]);
    for (dir = 0; dir <= 1; dir++)
    {
      NFFT_INT n = 0;
      R att = 0;
      int rc = NF(tune_sigma)(N, 2 * N, dir, (R)GOALS[ig], &n, &att);
      int m = 0;
      if (rc == 1)
      {
        NF(tune)(N, n, 2 * N, dir, (R)GOALS[ig], &m, 0);
        printf(" %8.3f m=%-3d", (double)n / (double)N, m);
      }
      else
        printf(" %14s", "out of reach");
    }
    printf("\n");
  }
  printf("\n* goal not attainable; m is the smallest reaching the best "
         "available accuracy\n");
  return 0;
}
