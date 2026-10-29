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

/* Tests of the NSFFT: the fast transforms on the hyperbolic cross against the
 * direct ones, in 2D (own plain grids and FFTW plans of the NFFT blocks) and
 * 3D, for several levels J, with the relative error measure of the NSFFT
 * example. Only failures are printed. */

#include <stdio.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "nsfft.h"

#ifdef CHECK_NSFFT

#define NSFFT_TEST_M 6
#define NSFFT_TEST_SAMPLES 500
#define NSFFT_TEST_TOL 1.0E-5

static int nsfft_check_one(const int d, const int J)
{
  nsfft_plan p;
  C *ref_tr, *ref_ad;
  R err_tr, err_ad;
  int ok;

  nsfft_init(&p, d, J, NSFFT_TEST_SAMPLES, NSFFT_TEST_M, NSDFT);
  ref_tr = (C*)NFFT(malloc)((size_t)(p.M_total) * sizeof(C));
  ref_ad = (C*)NFFT(malloc)((size_t)(p.N_total) * sizeof(C));

  NFFT(srand48)(20261029L + 16 * d + J);
  nsfft_init_random_nodes_coeffs(&p);

  nsfft_trafo_direct(&p);
  memcpy(ref_tr, p.f, (size_t)(p.M_total) * sizeof(C));
  nsfft_trafo(&p);
  err_tr = NFFT(error_l_infty_1_complex)(ref_tr, p.f, p.M_total, p.f_hat,
    p.N_total);

  NFFT(vrand_unit_complex)(p.f, p.M_total);
  nsfft_adjoint_direct(&p);
  memcpy(ref_ad, p.f_hat, (size_t)(p.N_total) * sizeof(C));
  nsfft_adjoint(&p);
  err_ad = NFFT(error_l_infty_1_complex)(ref_ad, p.f_hat, p.N_total, p.f,
    p.M_total);

  ok = err_tr <= NSFFT_TEST_TOL && err_ad <= NSFFT_TEST_TOL;
  if (!ok)
    printf("nsfft d = %d, J = %d: trafo %.3e, adjoint %.3e (%.1e) -> FAIL\n",
      d, J, (double)err_tr, (double)err_ad, NSFFT_TEST_TOL);

  NFFT(free)(ref_tr);
  NFFT(free)(ref_ad);
  nsfft_finalize(&p);
  return ok;
}

void NFFT(check_nsfft)(void)
{
  int J, ok = 1;

  /* J = 1 is too small for the fast transform */
  for (J = 2; J <= 6; J++)
    ok &= nsfft_check_one(2, J);
  for (J = 2; J <= 5; J++)
    ok &= nsfft_check_one(3, J);

  CU_ASSERT(ok);
}

#endif
