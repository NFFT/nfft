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

/* $Id$ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "nfft3.h"
#include "infft.h"

static void simple_test_nfct_1d(void)
{
  int j,k;
  nfct_plan p;

  int N=14;
  int M=19;

  /** init an one dimensional plan */
  nfct_init_1d(&p,N,M);

  /** init pseudo random nodes */
  for(j = 0; j < p.d*p.M_total; j++)
    p.x[j] = 0.5 * ((double)rand()) / RAND_MAX;

  /** precompute psi, the entries of the matrix B */
  if( p.flags & PRE_PSI)
    nfct_precompute_psi( &p);

  /** init pseudo random Fourier coefficients and show them */
  for(k = 0; k < p.N_total; k++)
    p.f_hat[k] = (double)rand() / RAND_MAX;

  nfft_vpr_double(p.f_hat,p.N_total,"given Fourier coefficients, vector f_hat");

  /** direct trafo and show the result */
  nfct_trafo_direct(&p);
  nfft_vpr_double(p.f,p.M_total,"ndct, vector f");

  /** approx. trafo and show the result */
  nfct_trafo(&p);
  nfft_vpr_double(p.f,p.M_total,"nfct, vector f");

  /** approx. adjoint and show the result */
  nfct_adjoint_direct(&p);
  nfft_vpr_double(p.f_hat,p.N_total,"adjoint ndct, vector f_hat");

  /** approx. adjoint and show the result */
  nfct_adjoint(&p);
  nfft_vpr_double(p.f_hat,p.N_total,"adjoint nfct, vector f_hat");

  /** finalise the one dimensional plan */
  nfct_finalize(&p);
}

int main(void)
{
  system("clear");
  printf("computing one dimensional ndct, nfct and adjoint ndct, nfct\n\n");
  simple_test_nfct_1d();
  printf("\n\n");

  return 1;
}
