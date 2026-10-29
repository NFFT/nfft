/*
 * Copyright (c) 2026 The NFFT Authors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see <https://www.gnu.org/licenses/>.
 */

/* Permanent regression smoke test for a downstream consumer of the NFFT3
 * package (used by both the installed find_package path and the FetchContent
 * path, on Linux and Windows CI). Exercises the public double-precision API
 * end to end with a direct transform (no precomputation) and asserts a known
 * value, so a wrong-precision or otherwise mis-linked library fails loudly.
 *
 * Math: for the 1D NDFT with all Fourier coefficients set to 1, the transform
 * evaluated at the node x = 0 is f(0) = sum_k f_hat[k] = N_total. So with one
 * node at x = 0 we must get f[0] == N_total (real), 0 (imag), exactly up to
 * floating-point round-off. Returns 0 on success, 1 on failure. */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "nfft3.h"

int main(void)
{
  nfft_plan p;
  nfft_init_1d(&p, 16, 1); /* N = 16 Fourier coeffs, M = 1 node */

  p.x[0] = 0.0; /* single node at the origin */
  for (int k = 0; k < p.N_total; k++)
    p.f_hat[k] = 1.0 + 0.0 * I; /* all coefficients = 1 */

  nfft_trafo_direct(&p); /* slow/direct transform: no precompute */

  double re = creal(p.f[0]);
  double im = cimag(p.f[0]);
  nfft_finalize(&p);

  if (fabs(re - (double)16) > 1e-9 || fabs(im) > 1e-9)
  {
    fprintf(stderr, "use_nfft3: expected f[0]=16+0i, got %g%+gi\n", re, im);
    return 1;
  }
  return 0;
}
