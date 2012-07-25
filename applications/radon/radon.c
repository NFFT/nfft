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
 * \file radon.c
 * \brief NFFT-based discrete Radon transform.
 *
 * Computes the discrete Radon transform
 * \f[
 *    R_{\theta_t} f\left(\frac{s}{R}\right)
 *    = \sum_{r \in I_R} w_r \; \sum_{k \in I_N^2} f_{k}
 *        \mathrm{e}^{-2\pi\mathrm{I} k \; (\frac{r}{R}\theta_t)}
 *        \, \mathrm{e}^{2\pi\mathrm{i} r s / R}
 *    \qquad(t \in I_T, s \in I_R).
 * \f]
 * by taking the 2D-NFFT of \f$f_k\f$ (\f$k \in I_N^2\f$)
 * at the points \f$\frac{r}{R}\theta_t\f$ of the polar or linogram grid
 * followed by 1D-iFFTs for every direction \f$t \in T\f$,
 * where \f$w_r\f$ are the weights of the Dirichlet- or Fejer-kernel.
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

/** computes the NFFT-based discrete Radon transform of f
 *  on the grid given by gridfcn() with T angles and R offsets
 */
int Radon_trafo(int (*gridfcn)(), int T, int S, double *f, int NN, double *Rf)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */

  fftw_complex *fft;                    /**< variable for the fftw-1Ds        */
  fftw_plan my_fftw_plan;               /**< plan for the fftw-1Ds            */

  int t,r;                              /**< index for directions and offsets */
  double *x, *w;                        /**< knots and associated weights     */

  int N[2],n[2];
  int M=T*S;

  N[0]=NN; n[0]=2*N[0];
  N[1]=NN; n[1]=2*N[1];

  fft = (fftw_complex *)nfft_malloc(S*sizeof(fftw_complex));
  my_fftw_plan = fftw_plan_dft_1d(S,fft,fft,FFTW_BACKWARD,FFTW_MEASURE);

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


  /** init nodes from grid*/
  gridfcn(T,S,x,w);
  for(j=0;j<my_nfft_plan.M_total;j++)
  {
    my_nfft_plan.x[2*j+0] = x[2*j+0];
    my_nfft_plan.x[2*j+1] = x[2*j+1];
  }

  /** precompute psi, the entries of the matrix B */
  if(my_nfft_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_nfft_plan);

  if(my_nfft_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_nfft_plan);

  if(my_nfft_plan.nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&my_nfft_plan);

  /** init Fourier coefficients from given image */
  for(k=0;k<my_nfft_plan.N_total;k++)
    my_nfft_plan.f_hat[k] = f[k] + _Complex_I*0.0;

  /** NFFT-2D */
  nfft_trafo(&my_nfft_plan);

  /** FFTW-1Ds */
  for(t=0; t<T; t++)
  {
    fft[0]=0.0;
    for(r=-S/2+1; r<S/2; r++)
      fft[r+S/2] = KERNEL(r)*my_nfft_plan.f[t*S+(r+S/2)];

    nfft_fftshift_complex(fft, 1, &S);
    fftw_execute(my_fftw_plan);
    nfft_fftshift_complex(fft, 1, &S);

    for(r=0; r<S; r++)
      Rf[t*S+r] = creal(fft[r])/S;

/*    for(r=0; r<R/2; r++)
      Rf[t*R+(r+R/2)] = creal(cexp(-I*KPI*r)*fft[r]);
    for(r=0; r<R/2; r++)
      Rf[t*R+r] = creal(cexp(-I*KPI*r)*fft[r+R/2]);
 */
  }

  /** finalise the plans and free the variables */
  fftw_destroy_plan(my_fftw_plan);
  nfft_free(fft);
  nfft_finalize(&my_nfft_plan);
  nfft_free(x);
  nfft_free(w);
  return 0;
}

/** simple test program for the discrete Radon transform
 */
int main(int argc,char **argv)
{
  int (*gridfcn)();                     /**< grid generating function        */
  int T, S;                             /**< number of directions/offsets    */
  FILE *fp;
  int N;                                /**< image size                      */
  double *f, *Rf;

  if( argc!=5 )
  {
    printf("radon gridfcn N T R\n");
    printf("\n");
    printf("gridfcn    \"polar\" or \"linogram\" \n");
    printf("N          image size NxN            \n");
    printf("T          number of slopes          \n");
    printf("R          number of offsets         \n");
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

  f   = (double *)nfft_malloc(N*N*(sizeof(double)));
  Rf  = (double *)nfft_malloc(T*S*(sizeof(double)));

  /** load data */
  fp=fopen("input_data.bin","rb");
  if (fp==NULL)
    return(-1);
  fread(f,sizeof(double),N*N,fp);
  fclose(fp);

  /** Radon transform */
  Radon_trafo(gridfcn,T,S,f,N,Rf);

  /** write result */
  fp=fopen("sinogram_data.bin","wb+");
  if (fp==NULL)
    return(-1);
  fwrite(Rf,sizeof(double),T*S,fp);
  fclose(fp);

  /** free the variables */
  nfft_free(f);
  nfft_free(Rf);

  return EXIT_SUCCESS;
}
