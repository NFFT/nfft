/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* RNFFT, the NFFT of real-valued functions: trafo, adjoint and the c2r
 * solver. */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "nfft3.h"
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

static void simple_test_rnfft_1d(void)
{
  RNFFT(plan) p;
  const char *error_str;

  RNFFT(init_1d)(&p, 14, 19);
  NFFT(vrand_shifted_unit_double)(p.x, p.M_total);
  if (p.flags & PRE_ONE_PSI)
    RNFFT(precompute_one_psi)(&p);

  error_str = RNFFT(check)(&p);
  if (error_str != 0)
  {
    printf("Error in rnfft module: %s\n", error_str);
    RNFFT(finalize)(&p);
    return;
  }

  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);
  NFFT(vpr_complex)(p.f_hat, p.N_total,
    "given Fourier coefficients, vector f_hat");

  RNFFT(trafo_direct)(&p);
  NFFT(vpr_double)(p.f, p.M_total, "ndft, vector f");

  RNFFT(trafo)(&p);
  NFFT(vpr_double)(p.f, p.M_total, "rnfft, vector f");

  RNFFT(adjoint_direct)(&p);
  NFFT(vpr_complex)(p.f_hat, p.N_total, "adjoint ndft, vector f_hat");

  RNFFT(adjoint)(&p);
  NFFT(vpr_complex)(p.f_hat, p.N_total, "adjoint rnfft, vector f_hat");

  RNFFT(finalize)(&p);
}

static void simple_test_rnfft_2d(void)
{
  RNFFT(plan) p;
  NFFT_C *ref;
  NFFT_R *fref, err = NFFT_K(0.0), norm = NFFT_K(0.0);
  int N[2] = {32, 14}, n[2] = {64, 32}, M = 32 * 14, k;
  double t0, t1;
  const char *error_str;

  RNFFT(init_guru)(&p, 2, N, M, n, 7, PRE_PHI_HUT | PRE_PSI | MALLOC_F_HAT
    | MALLOC_X | MALLOC_F | FFTW_INIT | NFFT_SORT_NODES,
    FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);
  if (p.flags & PRE_ONE_PSI)
    RNFFT(precompute_one_psi)(&p);

  error_str = RNFFT(check)(&p);
  if (error_str != 0)
  {
    printf("Error in rnfft module: %s\n", error_str);
    RNFFT(finalize)(&p);
    return;
  }

  ref = NFFT(malloc)((size_t)p.N_total * sizeof(NFFT_C));
  fref = NFFT(malloc)((size_t)p.M_total * sizeof(NFFT_R));
  NFFT(vrand_unit_complex)(ref, p.N_total);

  for (k = 0; k < p.N_total; k++)
    p.f_hat[k] = ref[k];
  t0 = NFFT(clock_gettime_seconds)();
  RNFFT(trafo_direct)(&p);
  t1 = NFFT(clock_gettime_seconds)();
  for (k = 0; k < p.M_total; k++)
    fref[k] = p.f[k];
  printf("ndft took %.2E seconds.\n", t1 - t0);

  for (k = 0; k < p.N_total; k++)
    p.f_hat[k] = ref[k];
  t0 = NFFT(clock_gettime_seconds)();
  RNFFT(trafo)(&p);
  t1 = NFFT(clock_gettime_seconds)();
  for (k = 0; k < p.M_total; k++)
  {
    err = NFFT_M(fmax)(err, NFFT_M(fabs)(p.f[k] - fref[k]));
    norm = NFFT_M(fmax)(norm, NFFT_M(fabs)(fref[k]));
  }
  printf("rnfft took %.2E seconds, relative error %.2E.\n", t1 - t0,
    (double)(err / norm));

  NFFT(free)(fref);
  NFFT(free)(ref);
  RNFFT(finalize)(&p);
}

/* Reconstruct coefficients from real samples with the c2r solver. */
static void simple_test_solver_rnfft_1d(int N, int M, int iter)
{
  RNFFT(plan) p;
  SOLVER(plan_c2r) ip;
  int k, l;

  RNFFT(init_1d)(&p, N, M);
  NFFT(vrand_shifted_unit_double)(p.x, p.M_total);
  if (p.flags & PRE_ONE_PSI)
    RNFFT(precompute_one_psi)(&p);

  SOLVER(init_c2r)(&ip, (NFFT(mv_plan_c2r)*)&p);

  NFFT(vrand_shifted_unit_double)(ip.y, p.M_total);
  NFFT(vpr_double)(ip.y, p.M_total, "given data, vector y");

  for (k = 0; k < p.N_total; k++)
    ip.f_hat_iter[k] = NFFT_K(0.0);

  SOLVER(before_loop_c2r)(&ip);
  printf(" Initial residual r=%" NFFT__FES__ "\n",
    NFFT_M(sqrt)(ip.dot_r_iter));

  for (l = 0; l < iter; l++)
  {
    SOLVER(loop_one_step_c2r)(&ip);
    printf(" Iteration %2d, residual r=%" NFFT__FES__ "\n", l,
      NFFT_M(sqrt)(ip.dot_r_iter));
  }
  NFFT(vpr_complex)(ip.f_hat_iter, p.N_total,
    "approximate solution, vector f_hat_iter");

  SOLVER(finalize_c2r)(&ip);
  RNFFT(finalize)(&p);
}

int main(void)
{
  printf("1) computing a one dimensional ndft, rnfft and an adjoint rnfft\n\n");
  simple_test_rnfft_1d();

  printf("\n2) computing a two dimensional ndft and rnfft\n\n");
  simple_test_rnfft_2d();

  printf("\n3) solving a one dimensional inverse rnfft\n\n");
  simple_test_solver_rnfft_1d(8, 6, 9);

  return EXIT_SUCCESS;
}
