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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3.h"
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

static void accuracy_nsfft(int d, int J, int M, int m)
{
  NSFFT(plan) p;
  NFFT_C *swap_sndft_trafo, *swap_sndft_adjoint;

  NSFFT(init)(&p, d, J, M, m, NSDFT);

  swap_sndft_trafo=(NFFT_C*) NFFT(malloc)(p.M_total*
						 sizeof(NFFT_C));
  swap_sndft_adjoint=(NFFT_C*) NFFT(malloc)(p.N_total*
						   sizeof(NFFT_C));

  NSFFT(init_random_nodes_coeffs)(&p);

  /** direct trafo */
  NSFFT(trafo_direct)(&p);

  NFFT_CSWAP(swap_sndft_trafo,p.f);

  /** approx. trafo */
  NSFFT(trafo)(&p);

  printf("%5d\t %+.5" NFFT__FES__ "\t",J,
         NFFT(error_l_infty_1_complex)(swap_sndft_trafo, p.f, p.M_total,
                                 p.f_hat, p.N_total));
  fflush(stdout);

  NFFT(vrand_unit_complex)(p.f, p.M_total);

  /** direct adjoint */
  NSFFT(adjoint_direct)(&p);

  NFFT_CSWAP(swap_sndft_adjoint,p.f_hat);

  /** approx. adjoint */
  NSFFT(adjoint)(&p);

  printf("%+.5" NFFT__FES__ "\n",
         NFFT(error_l_infty_1_complex)(swap_sndft_adjoint, p.f_hat,
                                 p.N_total,
                                 p.f, p.M_total));
  fflush(stdout);

  NFFT(free)(swap_sndft_adjoint);
  NFFT(free)(swap_sndft_trafo);

  /** finalise the one dimensional plan */
  NSFFT(finalize)(&p);
}

static void time_nsfft(int d, int J, int M, unsigned test_nsdft, unsigned test_nfft)
{
  int r, N[d], n[d];
  double t, t_nsdft, t_nfft, t_nsfft;
  double t0, t1;

  NSFFT(plan) p;
  NFFT(plan) np;

  for(r=0;r<d;r++)
  {
    N[r]= NFFT(exp2i)(J+2);
    n[r]=(3*N[r])/2;
    /*n[r]=2*N[r];*/
  }

  /** init */
  NSFFT(init)(&p, d, J, M, 4, NSDFT);
  NSFFT(init_random_nodes_coeffs)(&p);

  /* transforms */
  if(test_nsdft)
  {
    t_nsdft=0;
    r=0;
    while(t_nsdft<0.1)
    {
      r++;
      t0 = NFFT(clock_gettime_seconds)();
      NSFFT(trafo_direct)(&p);
      t1 = NFFT(clock_gettime_seconds)();
      t = (t1 - t0);
      t_nsdft+=t;
    }
    t_nsdft/=r;
  }
  else
    t_nsdft=nan("");

  if(test_nfft)
  {
    NFFT(init_guru)(&np,d,N,M,n,6, FG_PSI| MALLOC_F_HAT| MALLOC_F| FFTW_INIT,
      FFTW_MEASURE);
    np.x=p.act_nfft_plan->x;
    NFFT(precompute_one_psi)(&np);
    NSFFT(cp)(&p, &np);

    t_nfft=0;
    r=0;
    while(t_nfft<0.1)
    {
      r++;
      t0 = NFFT(clock_gettime_seconds)();
      NFFT(trafo)(&np);
      t1 = NFFT(clock_gettime_seconds)();
      t = (t1 - t0);
      t_nfft+=t;
    }
    t_nfft/=r;

    NFFT(finalize)(&np);
  }
  else
  {
    t_nfft=nan("");
  }

  t_nsfft=0;
  r=0;
  while(t_nsfft<0.1)
    {
      r++;
      t0 = NFFT(clock_gettime_seconds)();
      NSFFT(trafo)(&p);
      t1 = NFFT(clock_gettime_seconds)();
      t = (t1 - t0);
      t_nsfft+=t;
    }
  t_nsfft/=r;

  printf("%d\t%.2e\t%.2e\t%.2e\n", J, t_nsdft, t_nfft, t_nsfft);

  fflush(stdout);

  /** finalise */
  NSFFT(finalize)(&p);
}


int main(int argc,char **argv)
{
  int d, J, M;

  if(argc<=2)
  {
    fprintf(stderr,"nsfft_test type d [first last trials]\n");
    return -1;
  }

  d=atoi(argv[2]);
  fprintf(stderr,"Testing the nfft on the hyperbolic cross (nsfft).\n");

  if(atoi(argv[1])==1)
  {
    fprintf(stderr,"Testing the accuracy of the nsfft vs. nsdft\n");
    fprintf(stderr,"Columns: d, E_{1,\\infty}(trafo) E_{1,\\infty}(adjoint)\n\n");
    for(J=1; J<10; J++)
      accuracy_nsfft(d, J, 1000, 6);
  }

  if(atoi(argv[1])==2)
  {
    fprintf(stderr,"Testing the computation time of the nsdft, nfft, and nsfft\n");
    fprintf(stderr,"Columns: d, J, M, t_nsdft, t_nfft, t_nsfft\n\n");
    for(J=atoi(argv[3]); J<=atoi(argv[4]); J++)
    {
      if(d==2)
	M=(J+4)*NFFT(exp2i)(J+1);
      else
	M=6*NFFT(exp2i)(J)*(NFFT(exp2i)((J+1)/2+1)-1)+NFFT(exp2i)(3*(J/2+1));

      if(d*(J+2)<=24)
	time_nsfft(d, J, M, 1, 1);
      else
	if(d*(J+2)<=24)
	  time_nsfft(d, J, M, 0, 1);
	else
	  time_nsfft(d, J, M, 0, 0);
    }
  }

  return 1;
}
