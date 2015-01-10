/*
 * Copyright (c) 2002, 2015 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "nfft3.h"
#include "infft.h"
//#include <time.h>

int main(void)
{
  NFFT(plan) p;
  const int N = 1000000;
  const int M = 1000000;
  ticks t0, t1;
  double t;

  printf("nthreads = %d\n", NFFT(get_num_threads)());

  /* init */
  FFTW(init_threads)();
  NFFT(init_1d)(&p,N,M);

  /* pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x,p.M_total);

  /* precompute psi, that is, the entries of the matrix B */
  t0 = getticks();
  if(p.flags & PRE_ONE_PSI)
      NFFT(precompute_one_psi)(&p);
  t1 = getticks();
  t = NFFT(elapsed_seconds(t1,t0));
  fprintf(stderr,"precompute elapsed time: %.3f seconds\n",t);

  /* pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat,p.N_total);

  /* transformation */
  t0 = getticks();
  NFFT(trafo)(&p);
  t1 = getticks();
  t = NFFT(elapsed_seconds)(t1,t0);
  fprintf(stderr,"compute    elapsed time: %.3f seconds\n",t);
  fflush(stderr);
//  NFFT(vpr_complex)(p.f,p.M_total,"ndft, vector f");

  /* cleanup */
  NFFT(finalize)(&p);
  FFTW(cleanup_threads)();

  return EXIT_SUCCESS;
}
