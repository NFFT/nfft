/**
 * \file polarFFT/polar_fft_test.c
 * \brief NFFT-based polar FFT and inverse.
 *
 * Computes the NFFT-based polar FFT and its inverse for various parameters.
 * \author Markus Fenn
 * \date 2006
 */
#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "util.h"
#include "nfft3.h"

/** 
 * \defgroup applications_polarFFT_polar polar_fft_test
 * \ingroup applications_polarFFT
 * \{
 */

/** Generates the points \f$x_{t,j}\f$ with weights \f$w_{t,j}\f$
 *  for the polar grid with \f$T\f$ angles and \f$R\f$ offsets.
 *
 *  The nodes of the polar grid lie on concentric circles around the origin.
 *  They are given for \f$(j,t)^{\top}\in I_R\times I_T\f$ by
 *  a signed radius \f$r_j := \frac{j}{R} \in [-\frac{1}{2},\frac{1}{2})\f$ and
 *  an angle \f$\theta_t := \frac{\pi t}{T} \in [-\frac{\pi}{2},\frac{\pi}{2})\f$
 *  as 
 *  \f[
 *    x_{t,j} := r_j\left(\cos\theta_t, \sin\theta_t\right)^{\top}\,.
 *  \f]
 *  The total number of nodes is \f$M=TR\f$,
 *  whereas the origin is included multiple times.
 *
 *  Weights are introduced to compensate for local sampling density variations.
 *  For every point in the sampling set, we associate a small surrounding area.
 *  In case of the polar grid, we choose small ring segments.
 *  The area of such a ring segment around \f$x_{t,j}\f$ (\f$j \ne 0\f$) is
 *  \f[
 *    w_{t,j}
 *    = \frac{\pi}{2TR^2}\left(\left(|j|+\frac{1}{2}\right)^2-
 *      \left(|j|-\frac{1}{2}\right)^2\right)
 *    = \frac{\pi |j| }{TR^2}\, .
 *  \f]
 *  The area of the small circle of radius \f$\frac{1}{2R}\f$ around the origin is
 *  \f$\frac{\pi}{4R^2}\f$.
 *  Divided by the multiplicity of the origin in the sampling set, we get
 *  \f$w_{t,0} := \frac{\pi}{4TR^2}\f$.
 *  Thus, the sum of all weights is \f$\frac{\pi}{4}(1+\frac{1}{R^2})\f$ and
 *  we divide by this value for normalization.
 */
int polar_grid(int T, int R, double *x, double *w)
{
  int t, r;
  double W=(double)T*(((double)R/2.0)*((double)R/2.0)+1.0/4.0);

  for(t=-T/2; t<T/2; t++)
  {
    for(r=-R/2; r<R/2; r++)
    {
      x[2*((t+T/2)*R+(r+R/2))+0] = (double)r/R*cos(PI*t/T);
      x[2*((t+T/2)*R+(r+R/2))+1] = (double)r/R*sin(PI*t/T);
      if (r==0)
        w[(t+T/2)*R+(r+R/2)] = 1.0/4.0/W;
      else
        w[(t+T/2)*R+(r+R/2)] = fabs((double)r)/W;
    }
  }

  return T*R;                           /** return the number of knots        */
}

/** discrete polar FFT */
int polar_dft(fftw_complex *f_hat, int NN, fftw_complex *f, int T, int R, int m)
{
  int j,k;                              /**< index for nodes and frequencies  */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */

  double *x, *w;                        /**< knots and associated weights     */

  int N[2],n[2];
  int M=T*R;                            /**< number of knots                  */

  N[0]=NN; n[0]=2*N[0];                 /**< oversampling factor sigma=2      */
  N[1]=NN; n[1]=2*N[1];                 /**< oversampling factor sigma=2      */

  x = (double *)malloc(2*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, m,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init nodes from polar grid*/
  polar_grid(T,R,x,w);
  for(j=0;j<my_nfft_plan.M_total;j++)
  {
    my_nfft_plan.x[2*j+0] = x[2*j+0];
    my_nfft_plan.x[2*j+1] = x[2*j+1];
  }

  /** init Fourier coefficients from given image */
  for(k=0;k<my_nfft_plan.N_total;k++)
    my_nfft_plan.f_hat[k] = f_hat[k];

  /** NDFT-2D */
  ndft_trafo(&my_nfft_plan);

  /** copy result */
  for(j=0;j<my_nfft_plan.M_total;j++)
    f[j] = my_nfft_plan.f[j];

  /** finalise the plans and free the variables */
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);

  return EXIT_SUCCESS;
}

/** NFFT-based polar FFT */
int polar_fft(fftw_complex *f_hat, int NN, fftw_complex *f, int T, int R, int m)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */

  double *x, *w;                        /**< knots and associated weights     */

  int N[2],n[2];
  int M=T*R;                            /**< number of knots                  */

  N[0]=NN; n[0]=2*N[0];                 /**< oversampling factor sigma=2      */
  N[1]=NN; n[1]=2*N[1];                 /**< oversampling factor sigma=2      */

  x = (double *)malloc(2*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, m,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init nodes from polar grid*/
  polar_grid(T,R,x,w);
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
    my_nfft_plan.f_hat[k] = f_hat[k];

  /** NFFT-2D */
  nfft_trafo(&my_nfft_plan);

  /** copy result */
  for(j=0;j<my_nfft_plan.M_total;j++)
    f[j] = my_nfft_plan.f[j];

  /** finalise the plans and free the variables */
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);

  return EXIT_SUCCESS;
}

/** inverse NFFT-based polar FFT */
int inverse_polar_fft(fftw_complex *f, int T, int R, fftw_complex *f_hat, int NN, int max_i, int m)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */
  infft_plan my_infft_plan;             /**< plan for the inverse nfft        */

  double *x, *w;                        /**< knots and associated weights     */
  int l;                                /**< index for iterations             */

  int N[2],n[2];
  int M=T*R;                            /**< number of knots                  */

  N[0]=NN; n[0]=2*N[0];                 /**< oversampling factor sigma=2      */
  N[1]=NN; n[1]=2*N[1];                 /**< oversampling factor sigma=2      */

  x = (double *)malloc(2*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, m,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init two dimensional infft plan */
  infft_init_advanced(&my_infft_plan,&my_nfft_plan, CGNR | PRECOMPUTE_WEIGHT );

  /** init nodes, given samples and weights */
  polar_grid(T,R,x,w);
  for(j=0;j<my_nfft_plan.M_total;j++)
  {
    my_nfft_plan.x[2*j+0] = x[2*j+0];
    my_nfft_plan.x[2*j+1] = x[2*j+1];
    my_infft_plan.y[j]    = f[j];
    my_infft_plan.w[j]    = w[j];
  }

  /** precompute psi, the entries of the matrix B */
  if(my_nfft_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_nfft_plan);

  if(my_nfft_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_nfft_plan);

  if(my_nfft_plan.nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&my_nfft_plan);

  /** initialise damping factors */
  if(my_infft_plan.flags & PRECOMPUTE_DAMP)
    for(j=0;j<my_infft_plan.mv->N[0];j++)
      for(k=0;k<my_infft_plan.mv->N[1];k++)
  {
    my_infft_plan.w_hat[j*my_infft_plan.mv->N[1]+k]=
        (sqrt(pow(j-my_infft_plan.mv->N[0]/2,2)+pow(k-my_infft_plan.mv->N[1]/2,2))>(my_infft_plan.mv->N[0]/2)?0:1);
  }

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_nfft_plan.N_total;k++)
    my_infft_plan.f_hat_iter[k] = 0.0 + _Complex_I*0.0;

  /** solve the system */
  infft_before_loop(&my_infft_plan);

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
      infft_loop_one_step(&my_infft_plan);
    }
  }

  /** copy result */
  for(k=0;k<my_nfft_plan.N_total;k++)
    f_hat[k] = my_infft_plan.f_hat_iter[k];

  /** finalise the plans and free the variables */
  infft_finalize(&my_infft_plan);
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);

  return EXIT_SUCCESS;
}

/** test program for various parameters */
int main(int argc,char **argv)
{
  int N;                                /**< mpolar FFT size NxN              */
  int T, R;                             /**< number of directions/offsets     */
  int M;                                /**< number of knots of mpolar grid   */
  double *x, *w;                        /**< knots and associated weights     */
  fftw_complex *f_hat, *f, *f_direct, *f_tilde;
  int k;
  int max_i;                            /**< number of iterations             */
  int m;
  double temp1, temp2, E_max=0.0;
  FILE *fp1, *fp2;
  char filename[30];

  if( argc!=4 )
  {
    printf("polar_fft_test N T R \n");
    printf("\n");
    printf("N          polar FFT of size NxN     \n");
    printf("T          number of slopes          \n");
    printf("R          number of offsets         \n");
    exit(-1);
  }

  N = atoi(argv[1]);
  T = atoi(argv[2]);
  R = atoi(argv[3]);
  printf("N=%d, polar grid with T=%d, R=%d => ",N,T,R);

  x = (double *)malloc(2*1.25*T*R*(sizeof(double)));
  w = (double *)malloc(1.25*T*R*(sizeof(double)));

  f_hat    = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N*N);
  f        = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*T*R);
  f_direct = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*T*R);
  f_tilde  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N*N);

  /** generate knots of mpolar grid */
  M=polar_grid(T,R,x,w); printf("M=%d.\n",M);

  /** load data */
  fp1=fopen("input_data_r.dat","r");
  fp2=fopen("input_data_i.dat","r");
  if (fp1==NULL)
    return(-1);
  for(k=0;k<N*N;k++)
  {
    fscanf(fp1,"%le ",&temp1);
    fscanf(fp2,"%le ",&temp2);
    f_hat[k]=temp1+ _Complex_I*temp2;
  }
  fclose(fp1);
  fclose(fp2);

  /** direct polar FFT */
    polar_dft(f_hat,N,f_direct,T,R,m);
  //  polar_fft(f_hat,N,f_direct,T,R,12);

  /** Test of the polar FFT with different m */
  printf("\nTest of the polar FFT: \n");
  fp1=fopen("polar_fft_error.dat","w+");
  for (m=1; m<=12; m++)
  {
    /** fast polar FFT */
    polar_fft(f_hat,N,f,T,R,m);

    /** compute error of fast polar FFT */
    E_max=nfft_error_l_infty_complex(f_direct,f,M);
    printf("m=%2d: E_max = %e\n",m,E_max);
    fprintf(fp1,"%e\n",E_max);
  }
  fclose(fp1);

  /** Test of the inverse polar FFT for different m in dependece of the iteration number*/
  for (m=3; m<=9; m+=3)
  {
    printf("\nTest of the inverse polar FFT for m=%d: \n",m);
    sprintf(filename,"polar_ifft_error%d.dat",m);
    fp1=fopen(filename,"w+");
    for (max_i=0; max_i<=100; max_i+=10)
    {
      /** inverse polar FFT */
      inverse_polar_fft(f_direct,T,R,f_tilde,N,max_i,m);

      /** compute maximum relative error */
      /* E_max=0.0;
      for(k=0;k<N*N;k++)
      {
        temp = cabs((f_hat[k]-f_tilde[k])/f_hat[k]);
        if (temp>E_max) E_max=temp;
      }
      */
       E_max=nfft_error_l_infty_complex(f_hat,f_tilde,N*N);
      printf("%3d iterations: E_max = %e\n",max_i,E_max);
      fprintf(fp1,"%e\n",E_max);
    }
    fclose(fp1);
  }

  /** free the variables */
  free(x);
  free(w);
  free(f_hat);
  free(f);
  free(f_direct);
  free(f_tilde);

  return 0;
}
/* \} */
