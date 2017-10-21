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

/** Swap two vectors. */
#define CSWAP(x,y) {double _Complex * NFFT_SWAP_temp__; \
  NFFT_SWAP_temp__=(x); (x)=(y); (y)=NFFT_SWAP_temp__;}

static void accuracy_nsfft(int d, int J, int M, int m)
{
  nsfft_plan p;
  double _Complex *swap_sndft_trafo, *swap_sndft_adjoint;

  nsfft_init(&p, d, J, M, m, NSDFT);

  swap_sndft_trafo=(double _Complex*) nfft_malloc(p.M_total*
						 sizeof(double _Complex));
  swap_sndft_adjoint=(double _Complex*) nfft_malloc(p.N_total*
						   sizeof(double _Complex));

  nsfft_init_random_nodes_coeffs(&p);

  /** direct trafo */
  nsfft_trafo_direct(&p);

  CSWAP(swap_sndft_trafo,p.f);

  /** approx. trafo */
  nsfft_trafo(&p);

  printf("%5d\t %+.5E\t",J,
         nfft_error_l_infty_1_complex(swap_sndft_trafo, p.f, p.M_total,
                                 p.f_hat, p.N_total));
  fflush(stdout);

  nfft_vrand_unit_complex(p.f, p.M_total);

  /** direct adjoint */
  nsfft_adjoint_direct(&p);

  CSWAP(swap_sndft_adjoint,p.f_hat);

  /** approx. adjoint */
  nsfft_adjoint(&p);

  printf("%+.5E\n",
         nfft_error_l_infty_1_complex(swap_sndft_adjoint, p.f_hat,
                                 p.N_total,
                                 p.f, p.M_total));
  fflush(stdout);

  nfft_free(swap_sndft_adjoint);
  nfft_free(swap_sndft_trafo);

  /** finalise the one dimensional plan */
  nsfft_finalize(&p);
}

static void time_nsfft(int d, int J, int M, unsigned test_nsdft, unsigned test_nfft)
{
  int r, N[d], n[d];
  double t, t_nsdft, t_nfft, t_nsfft;
  double t0, t1;

  nsfft_plan p;
  nfft_plan np;

  for(r=0;r<d;r++)
  {
    N[r]= nfft_exp2i(J+2);
    n[r]=(3*N[r])/2;
    /*n[r]=2*N[r];*/
  }

  /** init */
  nsfft_init(&p, d, J, M, 4, NSDFT);
  nsfft_init_random_nodes_coeffs(&p);

  /* transforms */
  if(test_nsdft)
  {
    t_nsdft=0;
    r=0;
    while(t_nsdft<0.1)
    {
      r++;
      t0 = nfft_clock_gettime_seconds();
      nsfft_trafo_direct(&p);
      t1 = nfft_clock_gettime_seconds();
      t = (t1 - t0);
      t_nsdft+=t;
    }
    t_nsdft/=r;
  }
  else
    t_nsdft=nan("");

  if(test_nfft)
  {
    nfft_init_guru(&np,d,N,M,n,6, FG_PSI| MALLOC_F_HAT| MALLOC_F| FFTW_INIT,
      FFTW_MEASURE);
    np.x=p.act_nfft_plan->x;
    if(np.flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&np);
    nsfft_cp(&p, &np);

    t_nfft=0;
    r=0;
    while(t_nfft<0.1)
    {
      r++;
      t0 = nfft_clock_gettime_seconds();
      nfft_trafo(&np);
      t1 = nfft_clock_gettime_seconds();
      t = (t1 - t0);
      t_nfft+=t;
    }
    t_nfft/=r;

    nfft_finalize(&np);
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
      t0 = nfft_clock_gettime_seconds();
      nsfft_trafo(&p);
      t1 = nfft_clock_gettime_seconds();
      t = (t1 - t0);
      t_nsfft+=t;
    }
  t_nsfft/=r;

  printf("%d\t%.2e\t%.2e\t%.2e\n", J, t_nsdft, t_nfft, t_nsfft);

  fflush(stdout);

  /** finalise */
  nsfft_finalize(&p);
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
	M=(J+4)*nfft_exp2i(J+1);
      else
	M=6*nfft_exp2i(J)*(nfft_exp2i((J+1)/2+1)-1)+nfft_exp2i(3*(J/2+1));

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
