/* Permanent per-precision regression fixture. Selects the mangled API at compile
 * time via -DNFFT_PREFIX=<prefix> (nfft_ / nfftf_ / nfftl_), so one source proves
 * any precision package links and runs. Same math as the double consumer
 * (use_nfft3.c): a 1-node NDFT at x=0 with all Fourier coefficients 1 yields
 * f[0] = N_total.
 * The tolerance is widened to 1e-3 (vs the double consumer's 1e-9) so the check
 * also holds in single (float) precision.
 * nfft3.h declares all three manglings, so this compiles regardless of the
 * library's build precision; linking the wrong-precision library fails at link
 * time (unresolved <prefix>init_1d), which is exactly what we want to catch.
 * <complex.h> must be included before nfft3.h so that fftw3.h sees _Complex_I /
 * complex / I and uses "R _Complex" rather than "R[2]" for the complex typedefs;
 * without it, __real__/__imag__ and scalar assignment would both fail. */
#include <complex.h>
#include <stdio.h>

#define NFFT_CAT_(a, b) a##b
#define NFFT_CAT(a, b) NFFT_CAT_(a, b)
#define API(name) NFFT_CAT(NFFT_PREFIX, name)

#include "nfft3.h"

int main(void)
{
  API(plan) p;
  API(init_1d)(&p, 16, 1);          /* N = 16 coeffs, M = 1 node */

  p.x[0] = 0.0;
  for (int k = 0; k < p.N_total; k++)
    p.f_hat[k] = 1.0;               /* real part 1; imag part 0 by C99 real-to-complex conversion */

  API(trafo_direct)(&p);

  double re = (double) __real__ p.f[0];
  double im = (double) __imag__ p.f[0];
  API(finalize)(&p);

  if (re < 15.999 || re > 16.001 || im < -0.001 || im > 0.001) {
    fprintf(stderr, "use_nfft3_precision: expected f[0]=16+0i, got %g%+gi\n", re, im);
    return 1;
  }
  return 0;
}
