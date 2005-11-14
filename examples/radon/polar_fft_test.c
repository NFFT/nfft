/*gcc -o polar_fft_test polar_fft_test.c -lnfft3 -lfftw3 -lm -I/home/mfenn/NFFT3_develop/lib/trunk/include -L/home/mfenn/NFFT3_develop/lib/trunk/.libs -L/usr/local/lib*/
#include <math.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

/** generates the points x with weights w
 *  for the polar grid with T angles and R offsets
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

  return 0;
}

int polar_fft(fftw_complex *f_hat, int NN, fftw_complex *f, int T, int R)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */

  double *x, *w;                        /**< knots and associated weights     */

  int N[2],n[2];
  int M=T*R;

  N[0]=NN; n[0]=2*N[0];
  N[1]=NN; n[1]=2*N[1];

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

int inverse_polar_fft(fftw_complex *f, int T, int R, fftw_complex *f_hat, int NN, int max_i)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */
  infft_plan my_infft_plan;             /**< plan for the inverse nfft        */

  double *x, *w;                        /**< knots and associated weights     */
  int l;                                /**< index for iterations             */

  int N[2],n[2];
  int M=T*R;

  N[0]=NN; n[0]=2*N[0];
  N[1]=NN; n[1]=2*N[1];

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

int main(int argc,char **argv)
{
  int N;                                /**< polar FFT size NxN              */
  int T, R;                             /**< number of directions/offsets    */
  fftw_complex *f_hat, *f, *f_hat2;
  int k;
  int max_i;                            /**< number of iterations            */
  double temp, E_max=0.0;
  FILE *fp;

  if( argc!=5 )
  {
    printf("polar_fft_test N T R max_i\n");
    printf("\n");
    printf("N          polar FFT of size NxN     \n");
    printf("T          number of slopes          \n");
    printf("R          number of offsets         \n");
    printf("max_i      number of iterations      \n");
    exit(-1);
  }

  N = atoi(argv[1]);
  T = atoi(argv[2]);
  R = atoi(argv[3]);
  printf("N=%d, polar grid with T=%d, R=%d. \n",N,T,R);
  max_i = atoi(argv[4]);

  f_hat  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N*N);
  f      = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*T*R);
  f_hat2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N*N);

  /** load data */
  fp=fopen("input_data.dat","r");
  if (fp==NULL)
    return(-1);
  for(k=0;k<N*N;k++)
  {
    fscanf(fp,"%le ",&temp);
    f_hat[k]=temp;
  }
  fclose(fp);

  // /** init pseudo random Fourier coefficients and show them */
  //for(k=0;k<N*N;k++)
  //  f_hat[k] = ((double)rand())/RAND_MAX + I*((double)rand())/RAND_MAX;

  /** polar FFT */
  polar_fft(f_hat,N,f,T,R);

  /** inverse polar FFT */
  inverse_polar_fft(f,T,R,f_hat2,N,max_i);

  /** compute maximum absolute error */
  for(k=0;k<N*N;k++)
  {
    temp = cabs(f_hat[k]-f_hat2[k]);
    if (temp>E_max) E_max=temp;
  }
  printf("E_max = %g\n",E_max);

  /** write result */
  fp=fopen("polar_fft_data.dat","w+");
  if (fp==NULL)
    return(-1);
  for(k=0;k<N*N;k++)
  {
    if ((k>0)&(!(k%N)))
      fprintf(fp,"\n");
    fprintf(fp," % 15.7e ",creal(f_hat2[k]));
  }
  fclose(fp);

    /** free the variables */
  free(f_hat);
  free(f);
  free(f_hat2);

  return 0;
}
