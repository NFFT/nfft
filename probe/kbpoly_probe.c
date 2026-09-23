/* Sweep the polynomial Kaiser-Bessel fit error, max |poly - kb_phi| / peak,
 * over m and degree, for the precision the library was built in. */
#include <stdio.h>
#include "config.h"
#include "nfft3.h"
#include "infft.h"

int main(void)
{
  const R sig[] = {K(1.25), K(2.0)};
  unsigned s;
  INT m, deg, i, l;

  printf("MANT_DIG=%d eps=%.3Le cap=%d\n", (int)MANT_DIG,
      (long double)Y(float_property)(NFFT_EPSILON), (int)KB_POLY_DEG_MAX);
  printf("sigma  m  deg  err/eps  err\n");

  for (s = 0; s < 2; s++)
    for (m = 2; m <= 14; m++)
      for (deg = 8; deg <= 28; deg++)
      {
        const INT N = 64, n = 2 * (INT)(sig[s] * (R)N / K(2.0));
        const R b = KPI * (K(2.0) - (R)N / (R)n);
        const R lt = Y(bessel_i0_logtail)((R)m * b);
        const R pki = EXP(-(R)m * b - lt);
        const R peak = Y(kb_phi)(b, lt, pki, (R)m, K(0.0));
        const INT cols = KB_POLY_COLS(m), w = 2 * m + 2;
        R *coef = (R*) Y(malloc)((size_t)((deg + 1) * cols) * sizeof(R));
        R got[64], worst = K(0.0);

        Y(kb_poly_fit)(coef, b, lt, pki, (R)m, m, deg);

        for (i = 0; i <= 1024; i++)
        {
          const R nx0 = (R)(m - 1) + (R)i / K(512.0);

          Y(kb_poly_run)(got, coef, m, deg, nx0);
          for (l = 0; l < w; l++)
          {
            const R e = FABS(got[l] - Y(kb_phi)(b, lt, pki, (R)m,
                nx0 - (R)l)) / peak;
            if (e > worst)
              worst = e;
          }
        }

        printf("%.2f %2d %3d %8.1f  %.2Le\n", (double)sig[s], (int)m,
            (int)deg, (double)(worst / Y(float_property)(NFFT_EPSILON)),
            (long double)worst);
        Y(free)(coef);
      }

  return 0;
}
