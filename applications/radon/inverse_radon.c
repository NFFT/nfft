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

/**
 * \file inverse_radon.c
 * \brief NFFT-based discrete inverse Radon transform.
 *
 * Computes the inverse of the discrete Radon transform
 * \f[
 *    R_{\theta_t} f\left(\frac{s}{R}\right)
 *    = \sum_{r \in I_R} w_r \; \sum_{k \in I_N^2} f_{k}
 *        \mathrm{e}^{-2\pi\mathrm{I} k \; (\frac{r}{R}\theta_t)}
 *        \, \mathrm{e}^{2\pi\mathrm{i} r s / R}
 *    \qquad(t \in I_T, s \in I_R).
 * \f]
 * given at the points \f$\frac{r}{R}\theta_t\f$ of the polar or linogram grid
 * and where \f$w_r\f$ are the weights of the Dirichlet- or Fejer-kernel
 * by 1D-FFTs and the 2D-iNFFT.
 * \author Markus Fenn
 * \date 2005
 */
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"

/** define weights of kernel function for discrete Radon transform */
/*#define KERNEL(r) 1.0 */
#define KERNEL(r) (1.0-fabs((double)(r))/((double)S/2))

/** generates the points x with weights w
 *  for the polar grid with T angles and R offsets
 */
static int polar_grid(int T, int S, double *x, double *w)
{
  int t, r;
  double W=(double)T*(((double)S/2.0)*((double)S/2.0)+1.0/4.0);

  for(t=-T/2; t<T/2; t++)
  {
    for(r=-S/2; r<S/2; r++)
    {
      x[2*((t+T/2)*S+(r+S/2))+0] = (double)r/S*cos(KPI*t/T);
      x[2*((t+T/2)*S+(r+S/2))+1] = (double)r/S*sin(KPI*t/T);
      if (r==0)
        w[(t+T/2)*S+(r+S/2)] = 1.0/4.0/W;
      else
        w[(t+T/2)*S+(r+S/2)] = fabs((double)r)/W;
    }
  }

  return 0;
}

/** generates the points x with weights w
 *  for the linogram grid with T slopes and R offsets
 */
static int linogram_grid(int T, int S, double *x, double *w)
{
  int t, r;
  double W=(double)T*(((double)S/2.0)*((double)S/2.0)+1.0/4.0);

  for(t=-T/2; t<T/2; t++)
  {
    for(r=-S/2; r<S/2; r++)
    {
      if(t<0)
      {
        x[2*((t+T/2)*S+(r+S/2))+0] = (double)r/S;
        x[2*((t+T/2)*S+(r+S/2))+1] = (double)4*(t+T/4)/T*r/S;
      }
      else
      {
        x[2*((t+T/2)*S+(r+S/2))+0] = -(double)4*(t-T/4)/T*r/S;
        x[2*((t+T/2)*S+(r+S/2))+1] = (double)r/S;
      }
      if (r==0)
        w[(t+T/2)*S+(r+S/2)] = 1.0/4.0/W;
      else
        w[(t+T/2)*S+(r+S/2)] = fabs((double)r)/W;
    }
  }

  return 0;
}

/** computes the inverse discrete Radon transform of Rf
 *  on the grid given by gridfcn() with T angles and R offsets
 *  by a NFFT-based CG-type algorithm
 */
int Inverse_Radon_trafo(int (*gridfcn)(), int T, int S, double *Rf, int NN, double *f, int max_i)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */
  solver_plan_complex my_infft_plan;             /**< plan for the inverse nfft        */

  fftw_complex *fft;                    /**< variable for the fftw-1Ds        */
  fftw_plan my_fftw_plan;               /**< plan for the fftw-1Ds            */

  int t,r;                              /**< index for directions and offsets */
  double *x, *w;                        /**< knots and associated weights     */
  int l;                                /**< index for iterations             */

  int N[2],n[2];
  int M=T*S;

  N[0]=NN; n[0]=2*N[0];
  N[1]=NN; n[1]=2*N[1];

  fft = (fftw_complex *)nfft_malloc(S*sizeof(fftw_complex));
  my_fftw_plan = fftw_plan_dft_1d(S,fft,fft,FFTW_FORWARD,FFTW_MEASURE);

  x = (double *)nfft_malloc(2*T*S*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)nfft_malloc(T*S*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, 4,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init two dimensional infft plan */
  solver_init_advanced_complex(&my_infft_plan,(nfft_mv_plan_complex*)(&my_nfft_plan), CGNR | PRECOMPUTE_WEIGHT);

  /** init nodes and weights of grid*/
  gridfcn(T,S,x,w);
  for(j=0;j<my_nfft_plan.M_total;j++)
  {
    my_nfft_plan.x[2*j+0] = x[2*j+0];
    my_nfft_plan.x[2*j+1] = x[2*j+1];
    if (j%S)
      my_infft_plan.w[j]    = w[j];
    else
      my_infft_plan.w[j]    = 0.0;
  }

  /** precompute psi, the entries of the matrix B */
  if(my_nfft_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_nfft_plan);

  if(my_nfft_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_nfft_plan);

  if(my_nfft_plan.nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&my_nfft_plan);

  /** compute 1D-ffts and init given samples and weights */
  for(t=0; t<T; t++)
  {
/*    for(r=0; r<R/2; r++)
       fft[r] = cexp(I*KPI*r)*Rf[t*R+(r+R/2)];
      for(r=0; r<R/2; r++)
       fft[r+R/2] = cexp(I*KPI*r)*Rf[t*R+r];
 */

    for(r=0; r<S; r++)
      fft[r] = Rf[t*S+r] + _Complex_I*0.0;

    nfft_fftshift_complex(fft, 1, &S);
    fftw_execute(my_fftw_plan);
    nfft_fftshift_complex(fft, 1, &S);

    my_infft_plan.y[t*S] = 0.0;
    for(r=-S/2+1; r<S/2; r++)
      my_infft_plan.y[t*S+(r+S/2)] = fft[r+S/2]/KERNEL(r);
  }

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_nfft_plan.N_total;k++)
    my_infft_plan.f_hat_iter[k] = 0.0 + _Complex_I*0.0;

  /** solve the system */
  solver_before_loop_complex(&my_infft_plan);

  if (max_i<1)
  {
    l=1;
    for(k=0;k<my_nfft_plan.N_total;k++)
      my_infft_plan.f_hat_iter[k] = my_infft_plan.p_hat_iter[k];
  }
  else
  {
    for(l=1;l<=max_i;l++)
    {
      solver_loop_one_step_complex(&my_infft_plan);
      /*if (sqrt(my_infft_plan.dot_r_iter)<=1e-12) break;*/
    }
  }
  /*printf("after %d iteration(s): weighted 2-norm of original residual vector = %g\n",l-1,sqrt(my_infft_plan.dot_r_iter));*/

  /** copy result */
  for(k=0;k<my_nfft_plan.N_total;k++)
    f[k] = creal(my_infft_plan.f_hat_iter[k]);

  /** finalise the plans and free the variables */
  fftw_destroy_plan(my_fftw_plan);
  nfft_free(fft);
  solver_finalize_complex(&my_infft_plan);
  nfft_finalize(&my_nfft_plan);
  nfft_free(x);
  nfft_free(w);
  return 0;
}

/** simple test program for the inverse discrete Radon transform
 */
int main(int argc,char **argv)
{
  int (*gridfcn)();                     /**< grid generating function        */
  int T, S;                             /**< number of directions/offsets    */
  FILE *fp;
  int N;                                /**< image size                      */
  double *Rf, *iRf;
  int max_i;                            /**< number of iterations            */

  if( argc!=6 )
  {
    printf("inverse_radon gridfcn N T R max_i\n");
    printf("\n");
    printf("gridfcn    \"polar\" or \"linogram\" \n");
    printf("N          image size NxN            \n");
    printf("T          number of slopes          \n");
    printf("R          number of offsets         \n");
    printf("max_i      number of iterations      \n");
    exit(-1);
  }

  if (strcmp(argv[1],"polar") == 0)
    gridfcn = polar_grid;
  else
    gridfcn = linogram_grid;

  N = atoi(argv[2]);
  T = atoi(argv[3]);
  S = atoi(argv[4]);
  /*printf("N=%d, %s grid with T=%d, R=%d. \n",N,argv[1],T,R);*/
  max_i = atoi(argv[5]);

  Rf  = (double *)nfft_malloc(T*S*(sizeof(double)));
  iRf = (double *)nfft_malloc(N*N*(sizeof(double)));

  /** load data */
  fp=fopen("sinogram_data.bin","rb");
  if (fp==NULL)
    return(-1);
  fread(Rf,sizeof(double),T*S,fp);
  fclose(fp);

  /** inverse Radon transform */
  Inverse_Radon_trafo(gridfcn,T,S,Rf,N,iRf,max_i);

  /** write result */
  fp=fopen("output_data.bin","wb+");
  if (fp==NULL)
    return(-1);
  fwrite(iRf,sizeof(double),N*N,fp);
  fclose(fp);

  /** free the variables */
  nfft_free(Rf);
  nfft_free(iRf);

  return EXIT_SUCCESS;
}
