/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <omp.h>

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"
//#include <time.h>

int main(void)
{
  nfft_plan p;
  const int N = 1000000;
  const int M = 1000000;
  ticks t0, t1;
  double t;

  printf("nthreads = %d\n", nfft_get_num_threads());

  /* init */
  fftw_init_threads();
  nfft_init_1d(&p,N,M);

  /* pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x,p.M_total);

  /* precompute psi, that is, the entries of the matrix B */
  t0 = getticks();
  if(p.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);
  t1 = getticks();
  t = nfft_elapsed_seconds(t1,t0);
  fprintf(stderr,"precompute elapsed time: %.3f seconds\n",t);

  /* pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat,p.N_total);

  /* transformation */
  t0 = getticks();
  nfft_trafo(&p);
  t1 = getticks();
  t = nfft_elapsed_seconds(t1,t0);
  fprintf(stderr,"compute    elapsed time: %.3f seconds\n",t);
  fflush(stderr);
//  nfft_vpr_complex(p.f,p.M_total,"ndft, vector f");

  /* cleanup */
  nfft_finalize(&p);
  fftw_cleanup_threads();

  return EXIT_SUCCESS;
}
