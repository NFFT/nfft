#include "nfft3.h"


void simple_test_infct_1d( )
{
  int j,k,l;                            /**< index for nodes, freqencies, iter*/
  nfct_plan my_plan;                    /**< plan for the nfct                */
  infct_plan my_iplan;                  /**< plan for the inverse nfct        */

  double *f_hat_C_orig, *f_C_orig;
  
  srand( 27031950);

  /** initialise an one dimensional plan */
  nfct_init_1d( &my_plan, 512, 10000);

  /** initialise my_iplan */
  infct_init( &my_iplan, &my_plan);

  /** init pseudo random nodes */
  for( j = 0; j < my_plan.M; j++)
    my_plan.x[j] = 0.5 * ((double)rand()) / RAND_MAX;

  /** precompute psi, the entries of the matrix B */
  if( my_plan.nfct_flags & PRE_PSI)
    nfct_precompute_psi( &my_plan);

  f_hat_C_orig = (double*)malloc( my_plan.N_L * sizeof( double));
  f_C_orig =  (double*)malloc( my_plan.M * sizeof( double));

  for( k = 0; k < my_plan.N_L; k++) {
    my_plan.f_hat_C[k] = ((double)rand()) / RAND_MAX;
    f_hat_C_orig[k] = my_plan.f_hat_C[k];
  }

  ndct_trafo( &my_plan);


  /** init pseudo random samples (real) and show them */
  for( j = 0; j < my_plan.M; j++) {
    my_iplan.given_f[j] = my_plan.f_C[j];
    f_C_orig[j] = my_plan.f_C[j];
  }
  

  /** initialise some guess f_hat_0 */
  for( k = 0; k < my_plan.N_L; k++)
    my_iplan.f_hat_iter[k] = 1.0;


  /** solve the system */
  infct_before_loop(&my_iplan);

  for( l = 0; l < 30; l++)
    {
      infct_loop_one_step( &my_iplan);
      
      printf("%d   %e\n", l,
        E_2_error( f_hat_C_orig, my_iplan.f_hat_iter, my_plan.N_L));
    }

  infct_finalize( &my_iplan);  
  nfct_finalize( &my_plan);  
  free( f_hat_C_orig);
  free( f_C_orig);
}



int main() { 

  printf("1) computing an one dimensional infct\n\n"); simple_test_infct_1d();

  return EXIT_SUCCESS;
}
