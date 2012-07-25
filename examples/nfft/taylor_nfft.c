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

/*! \file taylor_nfft.c
 *
 * \brief Testing the nfft againt a Taylor expansion based version.
 *
 * \author Stefan Kunis
 *
 * References: Time and memory requirements of the Nonequispaced FFT
 *
 */
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

typedef struct
{
  nfft_plan p;                          /**< used for fftw and data          */

  int *idx0;                            /**< index of next neighbour of x_j
                                             on the oversampled regular grid */
  double *deltax0;                      /**< distance to the grid point      */
} taylor_plan;

/**
 * Initialisation of a transform plan.
 *
 * \arg ths The pointer to a taylor plan
 * \arg N The multi bandwidth
 * \arg M The number of nodes
 * \arg n The fft length
 * \arg m The order of the Taylor expansion
 *
 * \author Stefan Kunis
 */
static void taylor_init(taylor_plan *ths, int N, int M, int n, int m)
{
  /* Note: no nfft precomputation! */
  nfft_init_guru((nfft_plan*)ths, 1, &N, M, &n, m,
                 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_PRESERVE_INPUT);

  ths->idx0=(int*)nfft_malloc(M*sizeof(int));
  ths->deltax0=(double*)nfft_malloc(M*sizeof(double));
}

/**
 * Precomputation of weights and indices in Taylor expansion.
 *
 * \arg ths The pointer to a taylor plan
 *
 * \author Stefan Kunis
 */
static void taylor_precompute(taylor_plan *ths)
{
  int j;

  nfft_plan* cths=(nfft_plan*)ths;

  for(j=0;j<cths->M_total;j++)
    {
      ths->idx0[j] = ((int)round((cths->x[j]+0.5)*cths->n[0]) +
                                 cths->n[0]/2)%cths->n[0];
      ths->deltax0[j] = cths->x[j] - (round((cths->x[j]+0.5)*cths->n[0]) /
                                      cths->n[0] - 0.5);
    }
}

/**
 * Destroys a transform plan.
 *
 * \arg ths The pointer to a taylor plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
static void taylor_finalize(taylor_plan *ths)
{
  nfft_free(ths->deltax0);
  nfft_free(ths->idx0);

  nfft_finalize((nfft_plan*)ths);
}

/**
 * Executes a Taylor-NFFT, see equation (1.1) in [Guide], computes fast and
 * approximate by means of a Taylor expansion
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 *
 * \arg ths The pointer to a taylor plan
 *
 * \author Stefan Kunis
 */
static void taylor_trafo(taylor_plan *ths)
{
  int j,k,l,ll;
  double _Complex *f, *f_hat, *g1;
  double *deltax;
  int *idx;

  nfft_plan *cths=(nfft_plan*)ths;

  for(j=0, f=cths->f; j<cths->M_total; j++)
    *f++ = 0;

  for(k=0; k<cths->n_total; k++)
    cths->g1[k]=0;

  for(k=-cths->N_total/2, g1=cths->g1+cths->n_total-cths->N_total/2,
      f_hat=cths->f_hat; k<0; k++)
    (*g1++)=cpow( - 2*KPI*_Complex_I*k,cths->m)* (*f_hat++);

  cths->g1[0]=cths->f_hat[cths->N_total/2];

  for(k=1, g1=cths->g1+1, f_hat=cths->f_hat+cths->N_total/2+1;
      k<cths->N_total/2; k++)
    (*g1++)=cpow( - 2*KPI*_Complex_I*k,cths->m)* (*f_hat++);

  for(l=cths->m-1; l>=0; l--)
    {
      for(k=-cths->N_total/2, g1=cths->g1+cths->n_total-cths->N_total/2;
          k<0; k++)
        (*g1++) /= (-2*KPI*_Complex_I*k);

      for(k=1, g1=cths->g1+1; k<cths->N_total/2; k++)
        (*g1++) /= (-2*KPI*_Complex_I*k);

      fftw_execute(cths->my_fftw_plan1);

      ll=(l==0?1:l);
      for(j=0, f=cths->f, deltax=ths->deltax0, idx=ths->idx0; j<cths->M_total;
          j++, f++)
	(*f) = ((*f) * (*deltax++) + cths->g2[*idx++]) /ll;
    }
}

/**
 * Compares NDFT, NFFT, and Taylor-NFFT
 *
 * \arg N The bandwidth
 * \arg N The number of nodes
 * \arg n The FFT-size for the NFFT
 * \arg m The cut-off for window function
 * \arg n_taylor The FFT-size for the Taylor-NFFT
 * \arg m_taylor The order of the Taylor approximation
 * \arg test_accuracy Flag for NDFT computation
 *
 * \author Stefan Kunis
 */
static void taylor_time_accuracy(int N, int M, int n, int m, int n_taylor,
                          int m_taylor, unsigned test_accuracy)
{
  int r;
  double t_ndft, t_nfft, t_taylor, t;
  double _Complex *swapndft = NULL;
  ticks t0, t1;

  taylor_plan tp;
  nfft_plan np;

  printf("%d\t%d\t",N, M);

  taylor_init(&tp,N,M,n_taylor,m_taylor);

  nfft_init_guru(&np, 1, &N, M, &n, m,
                 PRE_PHI_HUT| PRE_FG_PSI|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** share nodes, input, and output vectors */
  np.x=tp.p.x;
  np.f_hat=tp.p.f_hat;
  np.f=tp.p.f;

  /** output vector ndft */
  if(test_accuracy)
    swapndft=(double _Complex*)nfft_malloc(M*sizeof(double _Complex));

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(np.x, np.M_total);

  /** nfft precomputation */
  taylor_precompute(&tp);

  /** nfft precomputation */
  if(np.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&np);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(np.f_hat, np.N_total);

  /** NDFT */
  if(test_accuracy)
    {
      CSWAP(np.f,swapndft);

      t_ndft=0;
      r=0;
      while(t_ndft<0.01)
        {
          r++;
          t0 = getticks();
          nfft_trafo_direct(&np);
          t1 = getticks();
t = nfft_elapsed_seconds(t1,t0);
          t_ndft+=t;
        }
      t_ndft/=r;

      CSWAP(np.f,swapndft);
      printf("%.2e\t",t_ndft);
    }
  else
    printf("nan\t\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.01)
    {
      r++;
      t0 = getticks();
      nfft_trafo(&np);
      t1 = getticks();
t = nfft_elapsed_seconds(t1,t0);
      t_nfft+=t;
    }
  t_nfft/=r;

  printf("%.2f\t%d\t%.2e\t",((double)n)/N, m, t_nfft);

  if(test_accuracy)
    printf("%.2e\t",X(error_l_infty_complex)(swapndft, np.f, np.M_total));
  else
    printf("nan\t\t");

  /** TAYLOR NFFT */
  t_taylor=0;
  r=0;
  while(t_taylor<0.01)
    {
      r++;
      t0 = getticks();
      taylor_trafo(&tp);
      t1 = getticks();
t = nfft_elapsed_seconds(t1,t0);
      t_taylor+=t;
    }
  t_taylor/=r;


  printf("%.2f\t%d\t%.2e\t",((double)n_taylor)/N,m_taylor,t_taylor);

  if(test_accuracy)
    printf("%.2e\n",X(error_l_infty_complex)(swapndft, np.f, np.M_total));
  else
    printf("nan\t\n");

  fflush(stdout);

  /** finalise */
  if(test_accuracy)
    nfft_free(swapndft);

  nfft_finalize(&np);
  taylor_finalize(&tp);
}

int main(int argc,char **argv)
{
  int l,m,trial,N;

  if(argc<=2)
    {
      fprintf(stderr,"taylor_nfft type first last trials sigma_nfft sigma_taylor.\n");
      return -1;
    }

  fprintf(stderr,"Testing the Nfft & a Taylor expansion based version.\n\n");
  fprintf(stderr,"Columns: N, M, t_ndft, sigma_nfft, m_nfft, t_nfft, e_nfft");
  fprintf(stderr,", sigma_taylor, m_taylor, t_taylor, e_taylor\n");

  /* time vs. N=M */
  if(atoi(argv[1])==0)
    {
      fprintf(stderr,"Fixed target accuracy, timings.\n\n");
      for(l=atoi(argv[2]); l<=atoi(argv[3]); l++)
        for(trial=0; trial<atoi(argv[4]); trial++)
          if(l<=10)
            taylor_time_accuracy((1U<< l), (1U<< l), (int)(atof(argv[5])*
                                 (1U<< l)), 6, (int)(atof(argv[6])*(1U<< l)),
                                 6, 1);
          else
            taylor_time_accuracy((1U<< l), (1U<< l), (int)(atof(argv[5])*
                                 (1U<< l)), 6, (int)(atof(argv[6])*(1U<< l)),
                                 6, 0);
    }

  /* error vs. m */
  if(atoi(argv[1])==1)
    {
      N=atoi(argv[7]);
      fprintf(stderr,"Fixed N=M=%d, error vs. m.\n\n",N);
      for(m=atoi(argv[2]); m<=atoi(argv[3]); m++)
        for(trial=0; trial<atoi(argv[4]); trial++)
          taylor_time_accuracy(N,N, (int)(atof(argv[5])*N), m,
                                    (int)(atof(argv[6])*N), m, 1);
    }

  return 1;
}
