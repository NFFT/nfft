/*
 * Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
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

/* $Id$ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

void simple_test_nfst_1d(void)
{
  int j,k;
  nfst_plan p;

  int N=14;
  int M=19;

  /** init an one dimensional plan */
  nfst_init_1d(&p,N,M);

  /** init pseudo random nodes */
  for(j = 0; j < p.d*p.M_total; j++)
    p.x[j] = 0.5 * ((double)rand()) / RAND_MAX;

  /** precompute psi, the entries of the matrix B */
  if( p.nfst_flags & PRE_PSI)
    nfst_precompute_psi( &p);

  /** init pseudo random Fourier coefficients and show them */
  for(k = 0; k < p.N_total; k++)
    p.f_hat[k] = (double)rand() / RAND_MAX;

  nfft_vpr_double(p.f_hat,p.N_total,"given Fourier coefficients, vector f_hat");

  /** direct trafo and show the result */
  ndst_trafo(&p);
  nfft_vpr_double(p.f,p.M_total,"ndst, vector f");

  /** approx. trafo and show the result */
  nfst_trafo(&p);
  nfft_vpr_double(p.f,p.M_total,"nfst, vector f");

  /** approx. adjoint and show the result */
  ndst_adjoint(&p);
  nfft_vpr_double(p.f_hat,p.N_total,"adjoint ndst, vector f_hat");

  /** approx. adjoint and show the result */
  nfst_adjoint(&p);
  nfft_vpr_double(p.f_hat,p.N_total,"adjoint nfst, vector f_hat");

  /** finalise the one dimensional plan */
  nfst_finalize(&p);
}

int main(void)
{
  system("clear");
  printf("computing one dimensional ndst, nfst and adjoint ndst, nfst\n\n");
  simple_test_nfst_1d();
  printf("\n\n");

  return 1;
}
