/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: simple_test.c 3198 2009-05-27 14:16:50Z keiner $ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"
#include <time.h>

int main(void)
{
  nfft_plan p;
  const int N = 10000000;
  const int M = 100000;
//  const int N = 20000;
//  const int M = 20000;
  time_t t0, t1;

  fftw_init_threads();

  fftw_plan_with_nthreads(2);

  /* one dimensional plan */
  nfft_init_1d(&p,N,M);

  /* pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x,p.M_total);

  /* precompute psi, that is, the entries of the matrix B */
  if(p.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);

  /* pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat,p.N_total);

  /** direct transformation */
  t0 = time(0);
  nfft_trafo(&p);
  t1 = time(0);
  fprintf(stderr,"elapsed time: %d seconds\n",t1-t0);
  fflush(stderr);
//  nfft_vpr_complex(p.f,p.M_total,"ndft, vector f");

  /* approximate transformation */
//  nfft_trafo(&p);
//  nfft_vpr_complex(p.f,p.M_total,"nfft, vector f");

  /** finalize plan */
  nfft_finalize(&p);

  fftw_cleanup_threads();

  return EXIT_SUCCESS;
}
