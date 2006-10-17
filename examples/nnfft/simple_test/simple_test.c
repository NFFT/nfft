#include "util.h"
#include "nfft3.h"
#include "stdlib.h"

void simple_test_nnfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nnfft_plan my_plan;                    /**< plan for the nfft                */

  int N[1];
  N[0]=12;

  /** init an one dimensional plan */
  nnfft_init(&my_plan, 1, 12, 19, N);
  
  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
  }  
  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[j]=((double)rand())/RAND_MAX-0.5;
  }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    nnfft_precompute_psi(&my_plan);
      
  if(my_plan.nnfft_flags & PRE_FULL_PSI)
    nnfft_precompute_full_psi(&my_plan);
    
  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    nnfft_precompute_lin_psi(&my_plan);
  
  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);
    
  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX +I*((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,my_plan.N_total,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  nndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nndft, vector f"); 

  /** approx. trafo and show the result */
  nnfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nnfft, vector f");

  /** finalise the one dimensional plan */
  nnfft_finalize(&my_plan);
}

void simple_test_nnfft_2d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nnfft_plan my_plan;                    /**< plan for the nfft                */

  int N[2];
  N[0]=12;
  N[1]=14;

  /** init an one dimensional plan */
  nnfft_init(&my_plan, 2,12*14,19, N);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
    my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
  }
    
  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[2*j]=((double)rand())/RAND_MAX-0.5;
    my_plan.v[2*j+1]=((double)rand())/RAND_MAX-0.5;
  }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    nnfft_precompute_psi(&my_plan);
      
  if(my_plan.nnfft_flags & PRE_FULL_PSI)
    nnfft_precompute_full_psi(&my_plan);
  
  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    nnfft_precompute_lin_psi(&my_plan);
        
  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);
    
  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I*((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,12,
        "given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  nndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nnfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** finalise the one dimensional plan */
  nnfft_finalize(&my_plan);
}

void simple_test_innfft_1d()
{
  int j,k,l,N=8;                        /**< index for nodes, freqencies, iter*/
  nnfft_plan my_plan;                    /**< plan for the nfft                */
  innfft_plan my_iplan;                  /**< plan for the inverse nfft        */

  /** initialise an one dimensional plan */
  nnfft_init(&my_plan,1 ,8 ,8 ,&N);

  /** initialise my_iplan */
  innfft_init_advanced(&my_iplan,&my_plan,CGNR);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
  
  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
    my_plan.v[j]=((double)rand())/RAND_MAX-0.5;
    
  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    nnfft_precompute_psi(&my_plan);
  
  if(my_plan.nnfft_flags & PRE_FULL_PSI)
      nnfft_precompute_full_psi(&my_plan);
  
  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);
    
  /** init pseudo random samples (real) and show them */
  for(j=0;j<my_plan.M_total;j++)
    my_iplan.y[j] = ((double)rand())/RAND_MAX;

  vpr_complex(my_iplan.y,my_plan.M_total,"given data, vector given_f");

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k] = 0.0;

  vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
        "approximate solution, vector f_hat_iter");

  /** solve the system */
  innfft_before_loop(&my_iplan);
  for(l=0;l<8;l++)
  {
    printf("iteration l=%d\n",l);
    innfft_loop_one_step(&my_iplan);
    vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
          "approximate solution, vector f_hat_iter");
      
    SWAP_complex(my_iplan.f_hat_iter,my_plan.f_hat);
    nnfft_trafo(&my_plan);
    vpr_complex(my_plan.f,my_plan.M_total,"fitting the data, vector f");
    SWAP_complex(my_iplan.f_hat_iter,my_plan.f_hat);
  }
  
  innfft_finalize(&my_iplan);  
  nnfft_finalize(&my_plan);  
}

void measure_time_nnfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nnfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N,N=4;
  double t;

  for(my_N=16; my_N<=16384; my_N*=2)
  {
    nnfft_init(&my_plan,1,my_N,my_N,&N);

    for(j=0;j<my_plan.M_total;j++)
      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;

    for(j=0;j<my_plan.N_total;j++)
      my_plan.v[j]=((double)rand())/RAND_MAX-0.5;        
        
    if(my_plan.nnfft_flags & PRE_PSI)
      nnfft_precompute_psi(&my_plan);
      
    if(my_plan.nnfft_flags & PRE_FULL_PSI)
        nnfft_precompute_full_psi(&my_plan);

    if(my_plan.nnfft_flags & PRE_PHI_HUT)
      nnfft_precompute_phi_hut(&my_plan);
      
    for(k=0;k<my_plan.N_total;k++)
      my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I*((double)rand())/RAND_MAX;

    t=second();
    nndft_trafo(&my_plan);
    t=second()-t;
    printf("t_nndft=%e,\t",t);

    t=second();
    nnfft_trafo(&my_plan);
    t=second()-t;
    printf("t_nnfft=%e\t",t);
      
    printf("(N=M=%d)\n",my_N);

    nnfft_finalize(&my_plan);  
  }
} 


int main()
{ 
  system("clear");
  printf("1) computing an one dimensional nndft, nnfft\n\n");
  simple_test_nnfft_1d();

  getc(stdin);

  system("clear"); 
  printf("2) computing a two dimensional nndft, nfft\n\n");
  simple_test_nnfft_2d();
  
  getc(stdin);
  
  system("clear"); 
  printf("3) computing an one dimensional infft\n\n");
  simple_test_innfft_1d();
  
  getc(stdin);
  
  system("clear"); 
  printf("4) computing times for one dimensional nnfft\n\n");
  measure_time_nnfft_1d();

  return 1;
}
