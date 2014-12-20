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

#undef X
#define X(name) NFST(name)

static void simple_test_nfst_1d(void)
{
  X(plan) p;

  int N = 14;
  int M = 19;

  /** init an one dimensional plan */
  X(init_1d)(&p,N,M);

  /** init pseudo random nodes */
  Y(vrand_real)(p.x, p.M_total, K(0.0), K(0.5));

  /** precompute psi, the entries of the matrix B */
  if( p.flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&p);

  /** init pseudo random Fourier coefficients and show them */
  Y(vrand_real)(p.f_hat, p.N_total, K(0.0), K(1.0));
  Y(vpr_double)(p.f_hat,p.N_total,"given Fourier coefficients, vector f_hat");

  /** direct trafo and show the result */
  X(trafo_direct)(&p);
  Y(vpr_double)(p.f,p.M_total,"ndst, vector f");

  /** approx. trafo and show the result */
  X(trafo)(&p);
  Y(vpr_double)(p.f,p.M_total,"nfst, vector f");

  /** approx. adjoint and show the result */
  X(adjoint_direct)(&p);
  Y(vpr_double)(p.f_hat,p.N_total,"adjoint ndst, vector f_hat");

  /** approx. adjoint and show the result */
  X(adjoint)(&p);
  Y(vpr_double)(p.f_hat,p.N_total,"adjoint nfst, vector f_hat");

  /** finalise the one dimensional plan */
  X(finalize)(&p);
}

int main(void)
{
  printf("Computing one dimensional ndct, nfct, adjoint ndct, and adjoint nfct...\n\n");
  simple_test_nfst_1d();
  printf("\n\n");

  return EXIT_SUCCESS;
}
