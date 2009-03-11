/*
 * Copyright (c) 2002, 2012 Jens Keiner, Daniel Potts, Stefan Kunis
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
#include <complex.h>

#include "util.h"
#include "nfft3.h"

void simple_test_nsfft(int d, int J, int M)
{
  int K=12;
  nsfft_plan p;

  nsfft_init(&p, d, J, M, 6, NSDFT);

  nsfft_init_random_nodes_coeffs(&p);

  nfft_vpr_complex(p.f_hat, K, "frequencies, vector f_hat (first few entries)");

  /** direct trafo and show the result */
  nsdft_trafo(&p);
  nfft_vpr_complex(p.f, K, "nsdft, vector f (first few entries)");

  /** approx. trafo and show the result */
  nsfft_trafo(&p);
  nfft_vpr_complex(p.f, K, "nsfft, vector f (first few entries)");

  /** direct adjoint and show the result */
  nsdft_adjoint(&p);
  nfft_vpr_complex(p.f_hat, K, "adjoint nsdft, vector f_hat, (first few entries)");

  /** approx. adjoint and show the result */
  nsfft_adjoint(&p);
  nfft_vpr_complex(p.f_hat, K, "adjoint nsfft, vector f_hat, (first few entries)");

  /** finalise the one dimensional plan */
  nsfft_finalize(&p);
}

int main(int argc,char **argv)
{
  int d, J, M;

  system("clear");
  printf("1) computing a two dimensional nsdft, nsfft and adjoints\n\n");
  d=2;
  J=5;
  M=(J+4)*nfft_int_2_pow(J+1);
  simple_test_nsfft(d,J,M);
  getc(stdin);

  system("clear");
  printf("2) computing a three dimensional nsdft, nsfft and adjoints\n\n");
  d=3;
  J=5;
  M=6*nfft_int_2_pow(J)*(nfft_int_2_pow((J+1)/2+1)-1)+nfft_int_2_pow(3*(J/2+1));
  simple_test_nsfft(d,J,M);

  return 1;
}
