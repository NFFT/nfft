#include "utils.h"
#include "nfft.h"

void simple_test_nfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */

  int N[1],n[1];
  N[0]=12; n[0]=32;

  /** init an one dimensional plan */
  /*nfft_init_1d(&my_plan,12,19);*/
  nfft_init_specific(&my_plan, 1, N, 19, n, 6,
		     PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		     FFT_OUT_OF_PLACE, FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
  
  /** init pseudo random nodes */
  for(j=0;j<my_plan.M;j++)
    {
      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
    }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_PSI)
    {
      nfft_precompute_psi(&my_plan);
      
      if(my_plan.nfft_flags & PRE_FULL_PSI)
	nfft_full_psi(&my_plan,pow(10,-10));
    }

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_L;k++)
    {
      my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
      my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
    }
  vpr_c(my_plan.f_hat,my_plan.N_L,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M,"nfft, vector f");

  /** approx. adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,my_plan.N_L,"adjoint ndft, vector f_hat");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,my_plan.N_L,"adjoint nfft, vector f_hat");

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
  nfft_init_specific(&my_plan, 2, N, 19, n, 6,
		     MALLOC_X| MALLOC_F_HAT| MALLOC_F| PRE_FULL_PSI, 
		     FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M;j++)
    {
      my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
      my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
    }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_PSI)
    {
      nfft_precompute_psi(&my_plan);
      
      if(my_plan.nfft_flags & PRE_FULL_PSI)
	nfft_full_psi(&my_plan,pow(10,-10));
    }

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_L;k++)
    {
      my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
      my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
    }
  vpr_c(my_plan.f_hat,12,
	"given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_c(my_plan.f,my_plan.M,"nfft, vector f");

  /** direct adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,12,"adjoint ndft, vector f_hat (first 12 entries)");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_c(my_plan.f_hat,12,"adjoint nfft, vector f_hat (first 12 entries)"); 

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}

void simple_test_infft_1d()
{
  int j,k,l;                            /**< index for nodes, freqencies, iter*/
  nfft_plan my_plan;                    /**< plan for the nfft                */
  infft_plan my_iplan;                  /**< plan for the inverse nfft        */
 
  /** initialise an one dimensional plan */
  nfft_init_1d(&my_plan, 18, 4);

  /** initialise my_iplan */
  infft_init(&my_iplan,&my_plan);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M;j++)
    my_plan.x[j]=((double)rand())/RAND_MAX-0.5;

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  /** init pseudo random samples (real) and show them */
  for(j=0;j<my_plan.M;j++)
    {
      my_iplan.given_f[j][0]=((double)rand())/RAND_MAX;
      my_iplan.given_f[j][1]=0.0;
    }
  vpr_c(my_iplan.given_f,my_plan.M,"given data, vector given_f");

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_plan.N_L;k++)
    {
      my_iplan.f_hat_iter[k][0]=0.0; 
      my_iplan.f_hat_iter[k][1]=0.0;
    }
  vpr_c(my_iplan.f_hat_iter,my_plan.N_L,
	"approximate solution, vector f_hat_iter");

  /** solve the system */
  infft_before_loop(&my_iplan);
  for(l=0;l<4;l++)
    {
      printf("iteration l=%d\n",l);
      infft_loop_one_step(&my_iplan);
      vpr_c(my_iplan.f_hat_iter,my_plan.N_L,
	    "approximate solution, vector f_hat_iter");
      
      SWAPC(my_iplan.f_hat_iter,my_plan.f_hat);
      nfft_trafo(&my_plan);
      vpr_c(my_plan.f,my_plan.M,"fitting the data, vector f");
      SWAPC(my_iplan.f_hat_iter,my_plan.f_hat);
    }
  
  infft_finalize(&my_iplan);  
  nfft_finalize(&my_plan);  
}

void measure_time_nfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N;
  double t;

  for(my_N=16; my_N<=16384; my_N*=2)
    {
      nfft_init_1d(&my_plan,my_N,my_N);

      for(j=0;j<my_plan.M;j++)
	my_plan.x[j]=((double)rand())/RAND_MAX-0.5;

      if(my_plan.nfft_flags & PRE_PSI)
	nfft_precompute_psi(&my_plan);

      for(k=0;k<my_plan.N_L;k++)
	{
	  my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
	  my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
	}

      t=second();
      ndft_trafo(&my_plan);
      t=second()-t;
      printf("t_ndft=%e,\t",t);

      t=second();
      nfft_trafo(&my_plan);
      t=second()-t;
      printf("t_nfft=%e\t",t);
      
      printf("(N=M=%d)\n",my_N);

      nfft_finalize(&my_plan);  
    }
} 


int main()
{ 
  system("clear"); 
  printf("1) computing an one dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_1d();

  getc(stdin);

  system("clear"); 
  printf("2) computing a two dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_2d();

  getc(stdin);

  system("clear"); 
  printf("3) computing an one dimensional infft\n\n");
  simple_test_infft_1d();

  getc(stdin);
  
  system("clear"); 
  printf("4) computing times for one dimensional nfft\n\n");
  measure_time_nfft_1d();

  return 1;
}
