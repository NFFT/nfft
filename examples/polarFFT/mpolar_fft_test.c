/**
 * \file mpolar_fft_test.c
 * \brief NFFT-based polar FFT and inverse on  modified polar grid.
 *
 * Computes the NFFT-based polar FFT and its inverse
 * on a modified polar grid for various parameters.
 * \author Markus Fenn
 * \date 2006
 */
#include <math.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

double GLOBAL_elapsed_time;

/** generates the points x with weights w
 *  for the mpolar grid with T angles and R offsets
 */
int mpolar_grid(int T, int R, double *x, double *w)
{
  int t, r;
  double W;
  int R2=2*ceil(sqrt(2)*R/2);
  double xx, yy;
  int M=0;

  for(t=-T/2; t<T/2; t++)
  {
    for(r=-R2/2; r<R2/2; r++)
    {
      xx = (double)r/R*cos(PI*t/T);
      yy = (double)r/R*sin(PI*t/T);

      if ( (-0.5-1.0/(double)R<=xx) & (xx<=0.5+1.0/(double)R & (-0.5-1.0/(double)R)<=yy) & (yy<=0.5+1.0/(double)R) )
      {
        x[2*M+0] = xx;
        x[2*M+1] = yy;

        if (r==0)
          w[M] = 1.0/4.0;
        else
          w[M] = fabs((double)r);

        M++; /** count the knots */
      }
    }
  }

   /** normalize the weights */
   W=0.0;
   for (t=0; t<M; t++)
      W+=w[t];

   for (t=0; t<M; t++)
    w[t]/=W;

  return M;                             /** return the number of knots        */
}

/** discrete mpolar FFT */
int mpolar_dft(fftw_complex *f_hat, int NN, fftw_complex *f, int T, int R, int m)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */

  double *x, *w;                        /**< knots and associated weights     */

  int N[2],n[2];
  int M;                                /**< number of knots                  */

  N[0]=NN; n[0]=2*N[0];                 /**< oversampling factor sigma=2      */
  N[1]=NN; n[1]=2*N[1];                 /**< oversampling factor sigma=2      */

  x = (double *)malloc(2*1.25*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(1.25*T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  M=mpolar_grid(T,R,x,w);
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, m,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init nodes from mpolar grid*/
  for(j=0;j<my_nfft_plan.M_total;j++)
  {
    my_nfft_plan.x[2*j+0] = x[2*j+0];
    my_nfft_plan.x[2*j+1] = x[2*j+1];
  }

  /** init Fourier coefficients from given image */
  for(k=0;k<my_nfft_plan.N_total;k++)
    my_nfft_plan.f_hat[k] = f_hat[k];

  GLOBAL_elapsed_time=second();

  /** NDFT-2D */
  ndft_trafo(&my_nfft_plan);

  GLOBAL_elapsed_time=second()-GLOBAL_elapsed_time;

  /** copy result */
  for(j=0;j<my_nfft_plan.M_total;j++)
    f[j] = my_nfft_plan.f[j];

  /** finalise the plans and free the variables */
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);

  return EXIT_SUCCESS;
}

/** NFFT-based mpolar FFT */
int mpolar_fft(fftw_complex *f_hat, int NN, fftw_complex *f, int T, int R, int m)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */

  double *x, *w;                        /**< knots and associated weights     */

  int N[2],n[2];
  int M;                                /**< number of knots                  */

  N[0]=NN; n[0]=2*N[0];                 /**< oversampling factor sigma=2      */
  N[1]=NN; n[1]=2*N[1];                 /**< oversampling factor sigma=2      */

  x = (double *)malloc(2*1.25*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(1.25*T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  M=mpolar_grid(T,R,x,w);
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, m,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init nodes from mpolar grid*/

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

  GLOBAL_elapsed_time=second();

  /** NFFT-2D */
  nfft_trafo(&my_nfft_plan);

  GLOBAL_elapsed_time=second()-GLOBAL_elapsed_time;

  /** copy result */
  for(j=0;j<my_nfft_plan.M_total;j++)
    f[j] = my_nfft_plan.f[j];

  /** finalise the plans and free the variables */
  nfft_finalize(&my_nfft_plan);
  free(x);
  free(w);

  return EXIT_SUCCESS;
}

/** inverse NFFT-based mpolar FFT */
int inverse_mpolar_fft(fftw_complex *f, int T, int R, fftw_complex *f_hat, int NN, int max_i, int m)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_nfft_plan;               /**< plan for the nfft-2D             */
  infft_plan my_infft_plan;             /**< plan for the inverse nfft        */

  double *x, *w;                        /**< knots and associated weights     */
  int l;                                /**< index for iterations             */

  int N[2],n[2];
  int M;                                /**< number of knots                  */

  N[0]=NN; n[0]=2*N[0];                 /**< oversampling factor sigma=2      */
  N[1]=NN; n[1]=2*N[1];                 /**< oversampling factor sigma=2      */

  x = (double *)malloc(2*1.25*T*R*(sizeof(double)));
  if (x==NULL)
    return -1;

  w = (double *)malloc(1.25*T*R*(sizeof(double)));
  if (w==NULL)
    return -1;

  /** init two dimensional NFFT plan */
  M=mpolar_grid(T,R,x,w);
  nfft_init_guru(&my_nfft_plan, 2, N, M, n, m,
                  PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init two dimensional infft plan */
    infft_init_advanced(&my_infft_plan,&my_nfft_plan, CGNR | PRECOMPUTE_WEIGHT );

  /** init nodes, given samples and weights */
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
      my_infft_plan.f_hat_iter[k] = 0.0 + I*0.0; 

  GLOBAL_elapsed_time=second();

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

  GLOBAL_elapsed_time=second()-GLOBAL_elapsed_time;

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

/** Comparison of the FFTW, mpolar FFT, and inverse mpolar FFT */
int comparison_fft(FILE *fp, int N, int T, int R)
{
  fftw_plan my_fftw_plan;
  fftw_complex *f_hat,*f;
  int m,k;
  double t_fft, t_dft_mpolar;

  f_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N*N);
  f     = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(T*R/4)*5);

  my_fftw_plan = fftw_plan_dft_2d(N,N,f_hat,f,FFTW_BACKWARD,FFTW_MEASURE);

  for(k=0; k<N*N; k++)
    f_hat[k] = drand48() + I* drand48();
  
  GLOBAL_elapsed_time=second();
  for(m=0;m<65536/N;m++)
    {
      fftw_execute(my_fftw_plan);
      /* touch */
      f_hat[2]=2*f_hat[0];
    }
  GLOBAL_elapsed_time=second()-GLOBAL_elapsed_time;
  t_fft=N*GLOBAL_elapsed_time/65536;

  if(N<256)
    {
      mpolar_dft(f_hat,N,f,T,R,m);
      t_dft_mpolar=GLOBAL_elapsed_time;
    }
      
  for (m=3; m<=9; m+=3)
    {
      if((m==3)&&(N<256))
        fprintf(fp,"%d\t&\t&\t%1.1e&\t%1.1e&\t%d\t",N,t_fft,t_dft_mpolar,m);
      else
        if(m==3)
	  fprintf(fp,"%d\t&\t&\t%1.1e&\t       &\t%d\t",N,t_fft,m);
	else
	  fprintf(fp,"  \t&\t&\t       &\t       &\t%\t",m);  

      printf("N=%d\tt_fft=%1.1e\tt_dft_mpolar=%1.1e\tm=%d\t",N,t_fft,t_dft_mpolar,m);              
   
      mpolar_fft(f_hat,N,f,T,R,m);
      fprintf(fp,"%1.1e&\t",GLOBAL_elapsed_time);
      printf("t_mpolar=%1.1e\t",GLOBAL_elapsed_time);
      inverse_mpolar_fft(f,T,R,f_hat,N,2*m,m);
      if(m==9)
	fprintf(fp,"%1.1e\\\\\\hline\n",GLOBAL_elapsed_time);
      else
	fprintf(fp,"%1.1e\\\\\n",GLOBAL_elapsed_time);
      printf("t_impolar=%1.1e\n",GLOBAL_elapsed_time);
    }

  fflush(fp);

  fftw_free(f);
  fftw_free(f_hat);  

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
  double temp, E_max=0.0;
  FILE *fp;
  char filename[30];
  int logN;

  if( argc!=4 )
  {
    printf("mpolar_fft_test N T R \n");
    printf("\n");
    printf("N          mpolar FFT of size NxN    \n");
    printf("T          number of slopes          \n");
    printf("R          number of offsets         \n");

    /** Hence, comparison of the FFTW, mpolar FFT, and inverse mpolar FFT */
    printf("\nHence, comparison FFTW, linogram FFT and inverse linogram FFT\n");
    fp=fopen("mpolar_comparison_fft.dat","w");
    if (fp==NULL)
	return(-1);
    for (logN=4; logN<=8; logN++)
	comparison_fft(fp,(1U<< logN), 3*(1U<< logN), 3*(1U<< (logN-1)));
    fclose(fp);

    exit(-1);
  }

  N = atoi(argv[1]);
  T = atoi(argv[2]);
  R = atoi(argv[3]);
  printf("N=%d, modified polar grid with T=%d, R=%d => ",N,T,R);

  x = (double *)malloc(2*1.25*T*R*(sizeof(double)));
  w = (double *)malloc(1.25*T*R*(sizeof(double)));

  f_hat    = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N*N);
  f        = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*1.25*T*R);  /* 4/pi*log(1+sqrt(2)) = 1.122... < 1.25 */
  f_direct = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*1.25*T*R);
  f_tilde  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N*N);

  /** generate knots of mpolar grid */
  M=mpolar_grid(T,R,x,w); printf("M=%d.\n",M);

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


  /** direct mpolar FFT */
      mpolar_dft(f_hat,N,f_direct,T,R,m);
  //  mpolar_fft(f_hat,N,f_direct,T,R,12);

  /** Test of the mpolar FFT with different m */
  printf("\nTest of the mpolar FFT: \n");
  fp=fopen("mpolar_fft_error.dat","w+");
  for (m=1; m<=12; m++)
  {
    /** fast mpolar FFT */
    mpolar_fft(f_hat,N,f,T,R,m);

    /** compute error of fast mpolar FFT */
    E_max=error_l_infty_complex(f_direct,f,M);
    printf("m=%2d: E_max = %e\n",m,E_max);
    fprintf(fp,"%e\n",E_max);
  }
  fclose(fp);

  /** Test of the inverse mpolar FFT for different m in dependece of the iteration number*/
  for (m=3; m<=9; m+=3)
  {
    printf("\nTest of the inverse mpolar FFT for m=%d: \n",m);
    sprintf(filename,"mpolar_ifft_error%d.dat",m);
    fp=fopen(filename,"w+");
    for (max_i=0; max_i<=20; max_i+=2)
    {
      /** inverse mpolar FFT */
      inverse_mpolar_fft(f_direct,T,R,f_tilde,N,max_i,m);

      /** compute maximum relativ error */
      E_max=error_l_infty_complex(f_hat,f_tilde,N*N);
      printf("%3d iterations: E_max = %e\n",max_i,E_max);
      fprintf(fp,"%e\n",E_max);
    }
    fclose(fp);
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
