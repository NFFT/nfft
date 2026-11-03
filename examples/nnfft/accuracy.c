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

#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3.h"
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

void accuracy(int d)
{
  int m,t;
  NNFFT(plan) my_plan;
  NFFT_C *slow;

  int N[d],n[d];
  int M_total,N_total;
  M_total=10000;N_total=1;

  slow=(NFFT_C*)NFFT(malloc)(M_total*sizeof(NFFT_C));

  for(t=0; t<d; t++)
    {
      N[t]=(1<<(12/d));
      n[t]=2*N[t];
      N_total*=N[t];
    }

  /** init a plan */
  for(m=0; m<10; m++)
    {
      NNFFT(init_guru)(&my_plan, d, N_total, M_total, N, n, m,
		      PRE_PSI| PRE_PHI_HUT|
		      MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F);


      /** init pseudo random nodes */
      NFFT(vrand_shifted_unit_double)(my_plan.x, d*my_plan.M_total);
      NFFT(vrand_shifted_unit_double)(my_plan.v, d*my_plan.N_total);

      /** precompute psi, the entries of the matrix B */
      if(my_plan.nnfft_flags & PRE_PSI)
      	NNFFT(precompute_psi)(&my_plan);

      if(my_plan.nnfft_flags & PRE_LIN_PSI)
        NNFFT(precompute_lin_psi)(&my_plan);

      if(my_plan.nnfft_flags & PRE_FULL_PSI)
        NNFFT(precompute_full_psi)(&my_plan);

      /** precompute psi, the entries of the matrix D */
      if(my_plan.nnfft_flags & PRE_PHI_HUT)
        NNFFT(precompute_phi_hut)(&my_plan);

      /** init pseudo random Fourier coefficients */
      NFFT(vrand_unit_complex)(my_plan.f_hat, my_plan.N_total);

      /** direct trafo and show the result */
      NNFFT(trafo_direct)(&my_plan);

      NFFT_CSWAP(my_plan.f,slow);

      /** approx. trafo and show the result */
      NNFFT(trafo)(&my_plan);

      printf(NFFT__FR__ ", " NFFT__FR__ "\n",
	     NFFT(error_l_infty_complex)(slow, my_plan.f, M_total),
	     NFFT(error_l_infty_1_complex)(slow, my_plan.f, M_total, my_plan.f_hat,
				     my_plan.N_total));

      /** finalise the one dimensional plan */
      NNFFT(finalize)(&my_plan);
    }
}

int main(void)
{
  int d;
  for(d=1; d<4; d++)
    accuracy(d);

  return 1;
}
