/**
 * \file radon.c
 * \brief NFFT-based discrete Radon transform and inverse.
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
 * followed by a 1D-iFFTs for every direction \f$t \in T\f$,
 * where \f$w_r\f$ are the weights of the Dirichlet- or Fejer-kernel.
 * The inverse transform is computed by 1D-FFTs and the 2D-iNFFT.
 * \author Markus Fenn
 * \date 2005
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

/** generates the points x with weights w
 *  for the linogram grid with T slopes and R offsets
 */
int linogram_grid(int T, int R, double *x, double *w)
{
  int t, r;
  double W=(double)T*(((double)R/2.0)*((double)R/2.0)+1.0/4.0);

  for(t=-T/2; t<T/2; t++)
  {
    for(r=-R/2; r<R/2; r++)
    {
      if(t<0)
      {
        x[2*((t+T/2)*R+(r+R/2))+0] = (double)r/R;
        x[2*((t+T/2)*R+(r+R/2))+1] = (double)4*(t+T/4)/T*r/R;
      }
      else
      {
        x[2*((t+T/2)*R+(r+R/2))+0] = -(double)4*(t-T/4)/T*r/R;
        x[2*((t+T/2)*R+(r+R/2))+1] = (double)r/R;
      }
      if (r==0)
        w[(t+T/2)*R+(r+R/2)] = 1.0/4.0/W;
      else
        w[(t+T/2)*R+(r+R/2)] = fabs((double)r)/W;
    }
  }

  return 0;
}

int Radon_trafo(int (*gridfcn)(), int T, int R, double *f, int NN, double *Rf)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */

  fftw_complex *fft;                    /**< variable for the fftw-1Ds        */
  fftw_plan my_fftw_plan;               /**< plan for the fftw-1Ds            */

  int t,r;                              /**< index for directions and offsets */
  double *x, *w;                        /**< knots and associated weights     */
  double W;

  int N[2],n[2];
  int M=T*R;

  N[0]=NN; n[0]=2*N[0];
  N[1]=NN; n[1]=2*N[1];

  fft = (fftw_complex *)fftw_malloc(R*sizeof(fftw_complex));
  my_fftw_plan = fftw_plan_dft_1d(R,fft,fft,FFTW_BACKWARD,FFTW_MEASURE);

  x = (double *)malloc(2*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, 4,
                 PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                 FFTW_MEASURE| FFTW_DESTROY_INPUT);


  /** init nodes from grid*/
  gridfcn(T,R,x,w);
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
    my_nfft_plan.f_hat[k] = f[k] + I*0.0;

  /** NFFT-2D */
  nfft_trafo(&my_nfft_plan);

  /** FFTW-1Ds */
  for(t=0; t<T; t++)
  {
    fft[0]=0.0;
    for(r=-R/2+1; r<R/2; r++)
      fft[r+R/2] = (1.0-fabs((double)r)/((double)R/2))*my_nfft_plan.f[t*R+(r+R/2)];
      //fft[r+R/2] = my_nfft_plan.f[t*R+(r+R/2)];

    fftshift_complex(fft, 1, &R);
    fftw_execute(my_fftw_plan);
    fftshift_complex(fft, 1, &R);

    for(r=0; r<R; r++)
      Rf[t*R+r] = creal(fft[r])/R;

//    for(r=0; r<R/2; r++)
//      Rf[t*R+(r+R/2)] = creal(cexp(-I*PI*r)*fft[r]);
//    for(r=0; r<R/2; r++)
//      Rf[t*R+r] = creal(cexp(-I*PI*r)*fft[r+R/2]);
  }

  /** finalise the plans and free the variables */
  fftw_destroy_plan(my_fftw_plan);
  fftw_free(fft);
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);
}

int Inverse_Radon_trafo(int (*gridfcn)(), int T, int R, double *Rf, int NN, double *f, int max_i)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */
  infft_plan my_infft_plan;             /**< plan for the inverse nfft        */

  fftw_complex *fft;                    /**< variable for the fftw-1Ds        */
  fftw_plan my_fftw_plan;               /**< plan for the fftw-1Ds            */

  int t,r;                              /**< index for directions and offsets */
  double *x, *w;                        /**< knots and associated weights     */
  double W;
  int l;                                /**< index for iterations             */

  int N[2],n[2];
  int M=T*R;

  N[0]=NN; n[0]=2*N[0];
  N[1]=NN; n[1]=2*N[1];

  fft = (fftw_complex *)fftw_malloc(R*sizeof(fftw_complex));
  my_fftw_plan = fftw_plan_dft_1d(R,fft,fft,FFTW_FORWARD,FFTW_MEASURE);

  x = (double *)malloc(2*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, 4,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init two dimensional infft plan */
  infft_init_advanced(&my_infft_plan,&my_nfft_plan, CGNR | PRECOMPUTE_WEIGHT);

  /** init nodes and weights of grid*/
  gridfcn(T,R,x,w);
  for(j=0;j<my_nfft_plan.M_total;j++)
  {
    my_nfft_plan.x[2*j+0] = x[2*j+0];
    my_nfft_plan.x[2*j+1] = x[2*j+1];
    if (j%R)
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
//     for(r=0; r<R/2; r++)
//       fft[r] = cexp(I*PI*r)*Rf[t*R+(r+R/2)];
//     for(r=0; r<R/2; r++)
//       fft[r+R/2] = cexp(I*PI*r)*Rf[t*R+r];

    for(r=0; r<R; r++)
      fft[r] = Rf[t*R+r] + I*0.0;

    fftshift_complex(fft, 1, &R);
    fftw_execute(my_fftw_plan);
    fftshift_complex(fft, 1, &R);

    my_infft_plan.y[t*R] = 0.0;
    for(r=-R/2+1; r<R/2; r++)
      my_infft_plan.y[t*R+(r+R/2)] = fft[r+R/2]/(1.0-fabs((double)r)/((double)R/2));
      //my_infft_plan.y[t*R+(r+R/2)] = fft[r+R/2];
  }

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_nfft_plan.N_total;k++)
    my_infft_plan.f_hat_iter[k] = 0.0 + I*0.0;

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
      //if (sqrt(my_infft_plan.dot_r_iter)<=1e-12) break;
    }
  }
  printf("after %d iteration(s): weighted 2-norm of original residual vector = %g\n",
         l-1,sqrt(my_infft_plan.dot_r_iter));

  /* copy result */
  for(k=0;k<my_nfft_plan.N_total;k++)
    f[k] = creal(my_infft_plan.f_hat_iter[k]);

  /** finalise the plans and free the variables */
  fftw_destroy_plan(my_fftw_plan);
  fftw_free(fft);
  infft_finalize(&my_infft_plan);
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);
}

int main(int argc,char **argv)
{
  int T, R;                             /**< number of directions/offsets    */
  FILE *fp;
  int N;                                /**< image size                      */
  double *f, *Rf, *iRf;
  int k;
  int max_i;                            /**< number of iterations            */

  if( argc!=5 )
  {
    printf("Radon_trafo N T R max_i\n");
    printf("\n");
    /*printf("gridfcn    \"polar\" or \"linogram\" \n");*/
    printf("N          image size NxN            \n");
    printf("T          number of slopes          \n");
    printf("R          number of offsets         \n");
    printf("max_i      number of iterations      \n");
    exit(-1);
  }

  N = atoi(argv[1]);
  T = atoi(argv[2]);
  R = atoi(argv[3]);
  printf("N=%d, T=%d, R=%d. \n",N,T,R);
  max_i = atoi(argv[4]);

  f   = (double *)malloc(N*N*(sizeof(double)));
  Rf  = (double *)malloc(T*R*(sizeof(double)));
  iRf = (double *)malloc(N*N*(sizeof(double)));

  /** load data */
  fp=fopen("input_data.bin","rb");
  if (fp==NULL)
    return(-1);
  fread(f,sizeof(double),N*N,fp);
  fclose(fp);

  /** Radon transform */
  Radon_trafo(linogram_grid,T,R,f,N,Rf);

  /** write result */
  fp=fopen("sinogram_data.bin","wb+");
  if (fp==NULL)
    return(-1);
  fwrite(Rf,sizeof(double),T*R,fp);
  fclose(fp);

  /** inverse Radon transform */
  Inverse_Radon_trafo(linogram_grid,T,R,Rf,N,iRf,max_i);

  /** write result */
  fp=fopen("output_data.bin","wb+");
  if (fp==NULL)
    return(-1);
  fwrite(iRf,sizeof(double),N*N,fp);
  fclose(fp);

  /** free the variables */
  free(f);
  free(Rf);
  free(iRf);

  return 0;
}
