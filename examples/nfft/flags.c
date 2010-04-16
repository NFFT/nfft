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

/* $Id$ */

/*! \file flags.c
 *
 * \brief Testing the nfft.
 *
 * \author Stefan Kunis
 *
 * References: Time and Memory Requirements of the Nonequispaced FFT
 */
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"

#ifdef GAUSSIAN
  unsigned test_fg=1;
#else
  unsigned test_fg=0;
#endif

#ifdef MEASURE_TIME_FFTW
  unsigned test_fftw=1;
#else
  unsigned test_fftw=0;
#endif

#ifdef MEASURE_TIME
  unsigned test=1;
#else
  unsigned test=0;
#endif

void flags_cp(nfft_plan *dst, nfft_plan *src)
{
  dst->x=src->x;
  dst->f_hat=src->f_hat;
  dst->f=src->f;
  dst->g1=src->g1;
  dst->g2=src->g2;
  dst->my_fftw_plan1=src->my_fftw_plan1;
  dst->my_fftw_plan2=src->my_fftw_plan2;
}

void time_accuracy(int d, int N, int M, int n, int m, unsigned test_ndft,
                   unsigned test_pre_full_psi)
{
  int r, NN[d], nn[d];
  double t_ndft, t, e;
  double _Complex *swapndft;
  ticks t0, t1;

  nfft_plan p;
  nfft_plan p_pre_phi_hut;
  nfft_plan p_fg_psi;
  nfft_plan p_pre_lin_psi;
  nfft_plan p_pre_fg_psi;
  nfft_plan p_pre_psi;
  nfft_plan p_pre_full_psi;

  printf("%d\t%d\t", d, N);

  for(r=0; r<d; r++)
    {
      NN[r]=N;
      nn[r]=n;
    }

  /** output vector ndft */
  if(test_ndft)
    swapndft=(double _Complex*)nfft_malloc(M*sizeof(double _Complex));

  nfft_init_guru(&p, d, NN, M, nn, m,
                 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  nfft_init_guru(&p_pre_phi_hut, d, NN, M, nn, m, PRE_PHI_HUT,0);
  flags_cp(&p_pre_phi_hut, &p);
  nfft_precompute_one_psi(&p_pre_phi_hut);

  if(test_fg)
    {
      nfft_init_guru(&p_fg_psi, d, NN, M, nn, m, FG_PSI,0);
      flags_cp(&p_fg_psi, &p);
      nfft_precompute_one_psi(&p_fg_psi);
    }

  nfft_init_guru(&p_pre_lin_psi, d, NN, M, nn, m, PRE_LIN_PSI,0);
  flags_cp(&p_pre_lin_psi, &p);
  nfft_precompute_one_psi(&p_pre_lin_psi);

  if(test_fg)
    {
      nfft_init_guru(&p_pre_fg_psi, d, NN, M, nn, m, PRE_FG_PSI,0);
      flags_cp(&p_pre_fg_psi, &p);
      nfft_precompute_one_psi(&p_pre_fg_psi);
    }

  nfft_init_guru(&p_pre_psi, d, NN, M, nn, m, PRE_PSI,0);
  flags_cp(&p_pre_psi, &p);
  nfft_precompute_one_psi(&p_pre_psi);

  if(test_pre_full_psi)
    {
      nfft_init_guru(&p_pre_full_psi, d, NN, M, nn, m, PRE_FULL_PSI,0);
      flags_cp(&p_pre_full_psi, &p);
      nfft_precompute_one_psi(&p_pre_full_psi);
    }

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** NDFT */
  if(test_ndft)
    {
      NFFT_SWAP_complex(p.f,swapndft);

      t_ndft=0;
      r=0;
      while(t_ndft<0.01)
        {
          r++;
          t0 = getticks();
          nfft_direct_trafo(&p);
          t1 = getticks();
          t = nfft_elapsed_seconds(t1,t0);
          t_ndft+=t;
        }
      t_ndft/=r;

      NFFT_SWAP_complex(p.f,swapndft);
    }
  else
    t_ndft=nan("");

  /** NFFTs */
  nfft_trafo(&p);
  nfft_trafo(&p_pre_phi_hut);
  if(test_fg)
    nfft_trafo(&p_fg_psi);
  else
    p_fg_psi.MEASURE_TIME_t[2]=nan("");
  nfft_trafo(&p_pre_lin_psi);
  if(test_fg)
    nfft_trafo(&p_pre_fg_psi);
  else
    p_pre_fg_psi.MEASURE_TIME_t[2]=nan("");
  nfft_trafo(&p_pre_psi);
  if(test_pre_full_psi)
    nfft_trafo(&p_pre_full_psi);
  else
    p_pre_full_psi.MEASURE_TIME_t[2]=nan("");

  if(test_fftw==0)
    p.MEASURE_TIME_t[1]=nan("");

  if(test_ndft)
    e=nfft_error_l_2_complex(swapndft, p.f, p.M_total);
  else
    e=nan("");

  printf("%.2e\t%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
         t_ndft,
         m,
         e,
         p.MEASURE_TIME_t[0],
         p_pre_phi_hut.MEASURE_TIME_t[0],

         p.MEASURE_TIME_t[1],

         p.MEASURE_TIME_t[2],
         p_fg_psi.MEASURE_TIME_t[2],
         p_pre_lin_psi.MEASURE_TIME_t[2],
         p_pre_fg_psi.MEASURE_TIME_t[2],
         p_pre_psi.MEASURE_TIME_t[2],
         p_pre_full_psi.MEASURE_TIME_t[2]);

  fflush(stdout);

  /** finalise */
  if(test_pre_full_psi)
    nfft_finalize(&p_pre_full_psi);
  nfft_finalize(&p_pre_psi);
  if(test_fg)
    nfft_finalize(&p_pre_fg_psi);
  nfft_finalize(&p_pre_lin_psi);
  if(test_fg)
    nfft_finalize(&p_fg_psi);
  nfft_finalize(&p_pre_phi_hut);
  nfft_finalize(&p);

  if(test_ndft)
    nfft_free(swapndft);
}

void accuracy_pre_lin_psi(int d, int N, int M, int n, int m, int K)
{
  int r, NN[d], nn[d];
  double e;
  double _Complex *swapndft;

  nfft_plan p;

  for(r=0; r<d; r++)
    {
      NN[r]=N;
      nn[r]=n;
    }

  /** output vector ndft */
  swapndft=(double _Complex*)nfft_malloc(M*sizeof(double _Complex));

  nfft_init_guru(&p, d, NN, M, nn, m,
                 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                 PRE_PHI_HUT| PRE_LIN_PSI|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** realloc psi */
  nfft_free(p.psi);
  p.K=K;
  p.psi=(double*) nfft_malloc((p.K+1)*p.d*sizeof(double));

  /** precomputation can be done before the nodes are initialised */
  nfft_precompute_one_psi(&p);

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** compute exact result */
  NFFT_SWAP_complex(p.f,swapndft);
  nfft_direct_trafo(&p);
  NFFT_SWAP_complex(p.f,swapndft);

  /** NFFT */
  nfft_trafo(&p);
  e=nfft_error_l_2_complex(swapndft, p.f, p.M_total);

  //  printf("%d\t%d\t%d\t%d\t%.2e\n",d,N,m,K,e);
  printf("$%.1e$&\t",e);

  fflush(stdout);

  /** finalise */
  nfft_finalize(&p);
  nfft_free(swapndft);
}


int main(int argc,char **argv)
{
  int l,m,d,trial,N;

  if(argc<=2)
    {
      fprintf(stderr,"flags type first last trials d m\n");
      return -1;
    }

  if((test==0)&&(atoi(argv[1])<2))
    {
      fprintf(stderr,"MEASURE_TIME in nfft3util.h not set\n");
      return -1;
    }

  fprintf(stderr,"Testing different precomputation schemes for the nfft.\n");
  fprintf(stderr,"Columns: d, N=M, t_ndft, e_nfft, t_D, t_pre_phi_hut, ");
  fprintf(stderr,"t_fftw, t_B, t_fg_psi, t_pre_lin_psi, t_pre_fg_psi, ");
  fprintf(stderr,"t_pre_psi, t_pre_full_psi\n\n");

  d=atoi(argv[5]);
  m=atoi(argv[6]);

  /* time vs. N=M */
  if(atoi(argv[1])==0)
    for(l=atoi(argv[2]); l<=atoi(argv[3]); l++)
      for(trial=0; trial<atoi(argv[4]); trial++)
	{
	  if(l<=20)
	    time_accuracy(d, (1U<< l), (1U<< (d*l)), (1U<< (l+1)), m, 0, 0);
	  else
	    time_accuracy(d, (1U<< l), (1U<< (d*l)), (1U<< (l+1)), m, 0, 0);
	}

  d=atoi(argv[5]);
  N=atoi(argv[6]);

  /* accuracy vs. time */
  if(atoi(argv[1])==1)
    for(m=atoi(argv[2]); m<=atoi(argv[3]); m++)
      for(trial=0; trial<atoi(argv[4]); trial++)
        time_accuracy(d, N, (int)pow(N,d), 2*N, m, 1, 1);

  d=atoi(argv[5]);
  N=atoi(argv[6]);
  m=atoi(argv[7]);

  /* accuracy vs. K for linear interpolation, assumes (m+1)|K */
  if(atoi(argv[1])==2)
    {
      printf("$\\log_2(K/(m+1))$&\t");
      for(l=atoi(argv[2]); l<atoi(argv[3]); l++)
	printf("$%d$&\t",l);

      printf("$%d$\\\\\n",atoi(argv[3]));

      printf("$\\tilde E_2$&\t");
      for(l=atoi(argv[2]); l<=atoi(argv[3]); l++)
	accuracy_pre_lin_psi(d, N, (int)pow(N,d), 2*N, m, (m+1)*(1U<< l));

      printf("\n");
    }

  return 1;
}
