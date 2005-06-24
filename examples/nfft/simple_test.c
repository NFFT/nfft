#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

void simple_test_nfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies  */
  nfft_plan my_plan;                    /**< plan for the nfft               */

  int N[1],n[1];
  N[0]=14; n[0]=32;

  /** init an one dimensional plan */
 nfft_init_guru(&my_plan, 1, N, 19, n, 6,
                 PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
    }

  /** precompute psi, the entries of the matrix B */
    if(my_plan.nfft_flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&my_plan);

  if(my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);
      
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,my_plan.N_total,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** approx. adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,my_plan.N_total,"adjoint ndft, vector f_hat");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,my_plan.N_total,"adjoint nfft, vector f_hat");

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}

void simple_test_nfft_2d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */

  int N[2],n[2];
  N[0]=12; n[0]=32;
  N[1]=14; n[1]=64;

  /** init an one dimensional plan */
  //nfft_init_2d(&my_plan,12,12,19);
  nfft_init_guru(&my_plan, 2, N, 19, n, 6,
		 PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT, FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
      my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
    }

  /** precompute psi, the entries of the matrix B */
    if(my_plan.nfft_flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&my_plan);

  if(my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);
      
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,12,
	"given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** direct adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint ndft, vector f_hat (first 12 entries)");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint nfft, vector f_hat (first 12 entries)"); 

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}

void simple_test_nfft_2d_huge()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */

  int N[2],n[2];
  N[0]=2050; n[0]=4096;
  N[1]=2050; n[1]=4096;

  /** init an one dimensional plan */

  nfft_init_guru(&my_plan, 2, N, 10, n, 
		 6, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT| MALLOC_X| MALLOC_F |
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
      my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
    }

  /** precompute psi, the entries of the matrix B */
    if(my_plan.nfft_flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&my_plan);

  if(my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);
      
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,12,
	"given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** direct adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint ndft, vector f_hat (first 12 entries)");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint nfft, vector f_hat (first 12 entries)"); 

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}



void measure_time_nfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N;
  double t;

  for(my_N=16; my_N<=1024*1024; my_N*=2)
    {
      nfft_init_1d(&my_plan,my_N,my_N);

      for(j=0;j<my_plan.M_total;j++)
	my_plan.x[j]=((double)rand())/RAND_MAX-0.5;

      if(my_plan.nfft_flags & PRE_PSI)
	nfft_precompute_psi(&my_plan);

      for(k=0;k<my_plan.N_total;k++)
	my_plan.f_hat[k]=((double)rand())/RAND_MAX + I*((double)rand())/RAND_MAX;

      if(my_N<=8192)
	{
	  t=second();
	  ndft_trafo(&my_plan);
	  t=second()-t;
	  printf("t_ndft=%e,\t",t);
	}
      else
	printf("t_ndft=nan\t");

      t=second();
      nfft_trafo(&my_plan);
      t=second()-t;
      printf("t_nfft=%e\t",t);
      
      printf("(N=M_total=%d)\n",my_N);

      nfft_finalize(&my_plan);  
    }
} 

void simple_test_infft_1d()
{
  int j,k,l;                            /**< index for nodes, freqencies,iter*/
  nfft_plan my_plan;                    /**< plan for the nfft               */
  infft_plan my_iplan;                  /**< plan for the inverse nfft       */
 
  /** initialise an one dimensional plan */
  nfft_init_1d(&my_plan, 8, 4);

  /** initialise my_iplan */
  infft_init(&my_iplan,&my_plan);

  /** init pseudo random nodes */
  vrand_shifted_unit_double(my_plan.x,my_plan.M_total);
  
  /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  /** init pseudo random samples and show them */
  vrand_unit_complex(my_iplan.y,my_plan.M_total);

  vpr_complex(my_iplan.y,my_plan.M_total,"given data, vector y");

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_plan.N_total;k++)
      my_iplan.f_hat_iter[k]=0; 

  vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
	"approximate solution, vector f_hat_iter");

  /** solve the system */
  infft_before_loop(&my_iplan);
  for(l=0;l<3;l++)
    {
      printf("iteration l=%d\n",l);
      infft_loop_one_step(&my_iplan);
      vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
	    "approximate solution, vector f_hat_iter");
      
      SWAP_complex(my_iplan.f_hat_iter,my_plan.f_hat);
      nfft_trafo(&my_plan);
      vpr_complex(my_plan.f,my_plan.M_total,"fitting the data, vector f");
      SWAP_complex(my_iplan.f_hat_iter,my_plan.f_hat);
    }
  
  infft_finalize(&my_iplan);  
  nfft_finalize(&my_plan);  
}

void simple_test_taylor_nfft_1d()
{
  int j,k,l;                            /**< index for nodes and freqencies  */
  nfft_plan my_plan, taylor_plan;       /**< plan for the nfft               */

  int N,M,m,n_taylor,facl;
  int* idx0;
  double* deltax0;
  int sigma;

  N=4096;
  M=4096;

  /** Taylor nfft */
  sigma=8;
  m=4;

  n_taylor=sigma*N;

  idx0=(int*)fftw_malloc(M*sizeof(int));
  deltax0=(double*)fftw_malloc(M*sizeof(double));

  /** init an one dimensional plan */
  nfft_init_1d(&my_plan, N, M);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
      idx0[j]=((int)round((my_plan.x[j]+0.5)*n_taylor) + n_taylor/2)%n_taylor;
      deltax0[j]=my_plan.x[j] - (round((my_plan.x[j]+0.5)*n_taylor) / n_taylor - 0.5);
    }

  /** precompute psi, the entries of the matrix B */
    if(my_plan.nfft_flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&my_plan);

  if(my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);
      
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,10,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,10,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,10,"nfft, vector f");  

  /** Taylor nfft */

  /** init an one dimensional plan */
  nfft_init_guru(&taylor_plan, 1, &N, M, &n_taylor, 0,
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_PRESERVE_INPUT);

  taylor_plan.x=my_plan.x;
  taylor_plan.f_hat=my_plan.f_hat;
  taylor_plan.f=my_plan.f;

  for(j=0; j<M; j++)
    taylor_plan.f[j]=0;

  for(k=0; k<n_taylor; k++)
    taylor_plan.g1[k]=0;

  for(l=0; l<m; l++)
    {
      facl=((l==0)?1:l*facl);

      for(k=-N/2; k<0; k++)
        taylor_plan.g1[n_taylor+k]=cpow( - 2*PI*I*k,l)*taylor_plan.f_hat[k+N/2];
      
      taylor_plan.g1[0]=taylor_plan.f_hat[N/2];

      for(k=1; k<N/2; k++)
        taylor_plan.g1[k]=cpow( - 2*PI*I*k,l)*taylor_plan.f_hat[k+N/2];

      fftw_execute(taylor_plan.my_fftw_plan1);

      for(j=0; j<M; j++)
        taylor_plan.f[j] += taylor_plan.g2[idx0[j]]*pow(deltax0[j],l)/facl;

      vpr_complex(taylor_plan.f,10,"taylor, vector f");  
    } 

  /** finalise the one dimensional plan */
  nfft_finalize(&taylor_plan);
  nfft_finalize(&my_plan);
  fftw_free(deltax0);
  fftw_free(idx0);
}

void nfft_taylor_nfft(nfft_plan *ths, double *deltax0, int *idx0)
{
  int j,k,l,ll;

  complex *f, *f_hat, *g1;
  double *deltax;
  int *idx;

  for(j=0, f=ths->f; j<ths->M_total; j++)
    *f++ = 0;

  for(k=-ths->N_total/2, g1=ths->g1+ths->n_total-ths->N_total/2, f_hat=ths->f_hat; k<0; k++)
    (*g1++)=cpow( - 2*PI*I*k,ths->m)* (*f_hat++);
   
  ths->g1[0]=ths->f_hat[ths->N_total/2];

  for(k=1, g1=ths->g1+1, f_hat=ths->f_hat+ths->N_total/2+1; k<ths->N_total/2; k++)
    (*g1++)=cpow( - 2*PI*I*k,ths->m)* (*f_hat++);

  for(l=ths->m-1; l>=0; l--)
    {
      for(k=-ths->N_total/2, g1=ths->g1+ths->n_total-ths->N_total/2; k<0; k++)
        (*g1++) /= (-2*PI*I*k);

      for(k=1, g1=ths->g1+1; k<ths->N_total/2; k++)
        (*g1++) /= (-2*PI*I*k);

      fftw_execute(ths->my_fftw_plan1);

      ll=(l==0?1:l);
      for(j=0, f=ths->f, deltax=deltax0, idx=idx0; j<ths->M_total; j++)
        (*f++) = ((*f) * (*deltax++) + ths->g2[*idx++]) /ll;
    }
}

void time_accuracy_taylor_nfft_1d(int N, int M, int n, int m, int n_taylor, int m_taylor, unsigned mode)
{
  int j,k,l,ll,r;                            /**< index for nodes and freqencies  */
  nfft_plan my_plan, taylor_plan;       /**< plan for the nfft               */

  int facl;
  int *idx0;
  double *deltax0;
  double t_ndft, t_nfft, t_taylor, t;
  complex *swapndft;
  double sigma,sigma_taylor;

  sigma_taylor=((double)n_taylor)/N;
  sigma=((double)n)/N;

  idx0=(int*)fftw_malloc(M*sizeof(int));
  deltax0=(double*)fftw_malloc(M*sizeof(double));
  swapndft=(complex*)fftw_malloc(M*sizeof(complex));

printf("%d\t%d\t%.1f\t%d\t%.1f\t%d\t",N, M, sigma, m, sigma_taylor, m_taylor);

  /** init an one dimensional plan PRE_PHI_HUT| PRE_FULL_PSI| */
  nfft_init_guru(&my_plan, 1, &N, M, &n, m,
                 PRE_PHI_HUT| PRE_FULL_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
      idx0[j]=((int)round((my_plan.x[j]+0.5)*n_taylor) + n_taylor/2)%n_taylor;
      deltax0[j]=my_plan.x[j] - (round((my_plan.x[j]+0.5)*n_taylor) / n_taylor - 0.5);
    }

  /** precompute psi, the entries of the matrix B */
    if(my_plan.nfft_flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&my_plan);

  if(my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);
      
  if(my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  /** direct trafo */
  if(mode)
    {
      SWAP_complex(my_plan.f,swapndft);
      
      t_ndft=0;
      r=0; 
      while(t_ndft<0.01)
        {
          r++;
          t=second();
          ndft_trafo(&my_plan);
          t=second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;

      SWAP_complex(my_plan.f,swapndft);
    }
  else 
    {
      t_ndft=0;
    }

  /** approx. trafo */
  t_nfft=0;
  r=0;
  while(t_nfft<0.01)
    {
      r++;
      t=second();
      nfft_trafo(&my_plan);
      t=second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;

  if(mode)
    printf("%.2e\t",error_l_infty_1_complex(swapndft, my_plan.f, my_plan.M_total,
                    my_plan.f_hat, my_plan.N_total));
  else
    printf("--------\t");

  /** Taylor nfft */

  /** init an one dimensional plan */
  nfft_init_guru(&taylor_plan, 1, &N, M, &n_taylor, m_taylor,
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_PRESERVE_INPUT);

  taylor_plan.x=my_plan.x;
  taylor_plan.f_hat=my_plan.f_hat;
  taylor_plan.f=my_plan.f;

  for(k=0; k<n_taylor; k++)
    taylor_plan.g1[k]=0;

  t_taylor=0;
  r=0;
  while(t_taylor<0.01)
    {
      r++;
      t=second();
      nfft_taylor_nfft(&taylor_plan, deltax0, idx0);
      t=second()-t;
      t_taylor+=t;
    }
  t_taylor/=r;

  if(mode)
    printf("%.2e\t",error_l_infty_1_complex(swapndft, my_plan.f, my_plan.M_total,
                    my_plan.f_hat, my_plan.N_total));
  else
    printf("--------\t");
  

  printf("%.2e\t%.2e\t%.2e\n",t_ndft,t_nfft,t_taylor);

  /** finalise the one dimensional plan */
  nfft_finalize(&taylor_plan);
  nfft_finalize(&my_plan);
  fftw_free(deltax0);
  fftw_free(idx0);
}

int main()
{
  int l,m;

  time_accuracy_taylor_nfft_1d((1U<< 4), 100, (1U<< 5), 3, (1U<< 7), 4, 1); exit(-1);

printf("polynomial degree N,\nnumber of nodes M,\nnfft oversampling factor sigma,\nnfft truncation parameter m,\n");
printf("taylor nfft oversampling factor D,\ntaylor nfft truncation parameter L,\nerrors e=|F-F_approx|_infty/|f|_1, and\ntimes t in sec.\n\n"); 

printf("N\tM\tsigma\tm\tD\tL\te_nfft\t\te_taylornfft\tt_ndft\t\tt_nfft\t\tt_taylornfft\n");

//  for(m=1;m<14;m++)
//    time_accuracy_taylor_nfft_1d((1U<< 12), 10000, (1U<< 13), m, (1U<< 16), m, 1);

  for(l=4;l<20;l++)
    if(l<13)
      time_accuracy_taylor_nfft_1d((1U<< l), (1U<< l), (1U<< (l+1)), 3, (1U<< (l+3)), 4, 1);
    else
      time_accuracy_taylor_nfft_1d((1U<< l), (1U<< l), (1U<< (l+1)), 3, (1U<< (l+3)), 4, 0);

  exit(-1);
//  simple_test_taylor_nfft_1d();  exit(-1);
/*  simple_test_nfft_2d_huge();
  exit(-1);
*/  
  system("clear");
  printf("1) computing an one dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_1d();

  getc(stdin);
/*
  system("clear"); 
  printf("2) computing a two dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_2d();

  getc(stdin);
*/
  system("clear"); 
  printf("3) computing an one dimensional infft\n\n");
  simple_test_infft_1d();
/*
  getc(stdin);
    
  system("clear"); 
  printf("3) computing times for one dimensional nfft\n\n");
  measure_time_nfft_1d();
*/
  return 1;
}
