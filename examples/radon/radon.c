/*gcc -o Radon_trafo Radon_trafo.c -lnfft3 -lfftw3 -lm -I/home/mfenn/NFFT3_develop/lib/trunk/include -L/home/mfenn/NFFT3_develop/lib/trunk/.libs -L/usr/local/lib*/
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
        x[2*((t+T/2)*R+(r+R/2))+0] = (double)-4*(t-T/4)/T*r/R;
        x[2*((t+T/2)*R+(r+R/2))+1] = (double)r/R;
      }
      if (r==0)
        w[(t+T/2)*R+(r+R/2)] = 1.0/4.0;
      else
        w[(t+T/2)*R+(r+R/2)] = fabs((double)r);
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
                 PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT,
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
    for(r=-R/2; r<R/2; r++)
      fft[r+R/2] = (1.0-fabs((double)r)/((double)R/2))*my_nfft_plan.f[t*R+(r+R/2)];
    fft[0]=0.0;

    fftw_execute(my_fftw_plan);

    for(r=0; r<R/2; r++)
      Rf[t*R+(r+R/2)] = creal(cexp(-I*PI*r)*fft[r]);
    for(r=0; r<R/2; r++)
      Rf[t*R+r] = creal(cexp(-I*PI*r)*fft[r+R/2]);
  }

  /** finalise the plans and free the variables */
  fftw_destroy_plan(my_fftw_plan);
  fftw_free(fft);
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);
}

int main(int argc,char **argv)
{
  int T, R;                             /**< number of directions/offsets    */
  FILE* fp;
  int N;                                /**< image size                      */
  double *f, *Rf;
  int k;

  if( argc!=4 )
  {
    printf("Radon_trafo N T R \n");
    printf("\n");
    /*printf("gridfcn    \"polar\" or \"linogram\" \n");*/
    printf("N          image size NxN            \n");
    printf("T          number of slopes          \n");
    printf("R          number of offsets         \n");
    exit(-1);
  }

  N = atoi(argv[1]);
  T = atoi(argv[2]);
  R = atoi(argv[3]);
  printf("N=%d, T=%d, R=%d. \n",N,T,R);

  f  = (double *)malloc(N*N*(sizeof(double)));
  Rf = (double *)malloc(T*R*(sizeof(double)));

  /** load data */
  fp=fopen("input_data.dat","r");
  if (fp==NULL)
    return(-1);
  for(k=0;k<N*N;k++)
    fscanf(fp,"%le ",&f[k]);
  fclose(fp);

  Radon_trafo(linogram_grid,T,R,f,N,Rf);

  /** write result */
  fp=fopen("output_data.dat","w+");
  if (fp==NULL)
    return(-1);

  for(k=0;k<T*R;k++)
  {
    if ((k>0)&(!(k%R)))
      fprintf(fp,"\n");
    fprintf(fp," % 15.7e",Rf[k]);
  }
  fclose(fp);

  return 0;
}
