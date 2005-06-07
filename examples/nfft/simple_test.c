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
  N[0]=16; n[0]=32;

  /** init an one dimensional plan */
  /*nfft_init_1d(&my_plan,12,19);*/
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

  vpr_c(my_plan.f_hat,my_plan.N_total,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** approx. adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,my_plan.N_total,"adjoint ndft, vector f_hat");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,my_plan.N_total,"adjoint nfft, vector f_hat");

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
		 PRE_PHI_HUT| PRE_LIN_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
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

  vpr_c(my_plan.f_hat,12,
	"given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** direct adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,12,"adjoint ndft, vector f_hat (first 12 entries)");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,12,"adjoint nfft, vector f_hat (first 12 entries)"); 

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

  vpr_c(my_plan.f_hat,12,
	"given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** direct adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,12,"adjoint ndft, vector f_hat (first 12 entries)");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,12,"adjoint nfft, vector f_hat (first 12 entries)"); 

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}



void measure_time_nfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N;
  double t;

  for(my_N=16; my_N<=1024*1024; my_N*=2)simple_test_nfft_2d();
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


int main()
{
  simple_test_nfft_2d_huge();
  exit(-1);
    
  printf("1) computing an one dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_1d();

  getc(stdin);

  system("clear"); 
  printf("2) computing a two dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_2d();

  getc(stdin);
  
  system("clear"); 
  printf("3) computing times for one dimensional nfft\n\n");
  measure_time_nfft_1d();

  return 1;
}
