/* Fit error max |poly - phi| / peak vs degree, window at reach m + 1,
 * same sampling as tests/kbpoly.c but 4x denser. Prints err in units of eps. */
#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "nfft3.h"
#include "infft.h"

int main(int argc, char **argv)
{
  const R sig[2] = {K(1.25), K(2.0)};
  const R eps = Y(float_property)(NFFT_EPSILON);
  const INT dmax = argc > 1 ? atoi(argv[1]) : 24;
  INT m, deg;
  unsigned s;

  printf("eps = %.3Le; rows m, cols degree 6..%td, err/eps (max over sigma)\n",
      (long double)eps, (ptrdiff_t)dmax);
  printf("  m");
  for (deg = 6; deg <= dmax; deg++)
    printf(" %6td", (ptrdiff_t)deg);
  printf("\n");

  for (m = 2; m <= 14; m++)
  {
    printf("%3td", (ptrdiff_t)m);
    for (deg = 6; deg <= dmax; deg++)
    {
      R worst = K(0.0);

      for (s = 0; s < 2; s++)
      {
        const INT N = 64, n = 2 * (INT)(sig[s] * (R)N / K(2.0));
        const R b = KPI * (K(2.0) - (R)N / (R)n);
        const R reach = (R)m + K(1.0);
        const R lt = Y(bessel_i0_logtail)(reach * b);
        const R pki = EXP(-reach * b - lt);
        const R peak = Y(kb_phi)(b, lt, pki, reach, K(0.0));
        const INT w = 2 * m + 2;
        R *coef = malloc((size_t)((deg + 1) * KB_POLY_COLS(m)) * sizeof(R));
        R *got = malloc((size_t)w * sizeof(R));
        INT i, l;

        Y(kb_poly_fit)(coef, b, lt, pki, reach, m, deg);

        for (i = 0; i <= 2048; i++)
        {
          const R nx0 = (R)(m - 1) + (R)i / K(1024.0);
          const R nx = (R)(m + 2) * ((R)i / K(1024.0) - K(1.0));
          R err;

          Y(kb_poly_run)(got, coef, m, deg, nx0);
          for (l = 0; l < w; l++)
          {
            err = FABS(got[l] - Y(kb_phi)(b, lt, pki, reach, nx0 - (R)l)) / peak;
            if (err > worst)
              worst = err;
          }
          err = FABS(Y(kb_poly_phi)(coef, m, deg, nx)
              - Y(kb_phi)(b, lt, pki, reach, nx)) / peak;
          if (err > worst)
            worst = err;
        }
        free(got);
        free(coef);
      }
      printf(" %6.3g", (double)(worst / eps));
    }
    printf("\n");
  }
  return 0;
}
