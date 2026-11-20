/* Consumer without <complex.h>: fftw3.h then makes NFFT_C a 2-element real
 * array. Same check as use_nfft3.c, f[0] = N_total for one node at x = 0 and
 * all Fourier coefficients 1. Precision comes from NFFT_PRECISION_*. */
#include <stdio.h>
#include "nfft3.h"
#include "nfft3mp.h"

#ifdef NFFT_CREAL
#error "nfft3mp.h defines the C99 complex helpers without <complex.h>"
#endif

int main(void)
{
  NFFT(plan) p;
  NFFT(init_1d)(&p, 16, 1);

  p.x[0] = NFFT_K(0.0);
  for (int k = 0; k < p.N_total; k++)
  {
    p.f_hat[k][0] = NFFT_K(1.0);
    p.f_hat[k][1] = NFFT_K(0.0);
  }

  NFFT(trafo_direct)(&p);

  double re = (double)p.f[0][0];
  double im = (double)p.f[0][1];
  NFFT(finalize)(&p);

  if (re < 15.999 || re > 16.001 || im < -0.001 || im > 0.001)
  {
    fprintf(stderr, "use_nfft3_nocomplex: expected f[0]=16+0i, got %g%+gi\n", re, im);
    return 1;
  }
  return 0;
}
