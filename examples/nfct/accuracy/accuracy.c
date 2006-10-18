#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

void accuracy(int d, int M, int m1, int m2)
{
  int j,k,t,m;                         
  nfct_plan my_plan;                   

  double *slow;
  double *my_f_hat;
  double *my_x;
  
  int N[d], n[d];

  srand( 27031950);

  for( t = 0; t < d; t++) {
    N[t] = (1<<(12/d - 1)<<1);
    n[t] = 2 * nfft_next_power_of_2( N[t]) + 1;
    
    printf( "N[%d]=%d\n", t, N[t]);
  }

  nfct_init_guru(&my_plan, d, N, M, n, 6,
     	  MALLOC_X| MALLOC_F_HAT| MALLOC_F|
          FFTW_INIT| FFT_OUT_OF_PLACE, 
          FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  slow = (double*)fftw_malloc( M * sizeof( double));
  my_x = (double*)fftw_malloc( d * M * sizeof( double));
  my_f_hat = (double*)fftw_malloc( my_plan.N_total * sizeof( double));


  /** init pseudo random nodes */
  for( j = 0; j < d*M; j++)
    my_x[j] = 0.5 * ((double)rand()) / RAND_MAX;

  /** init pseudo random Fourier coefficients and show them */
  for( k = 0; k < my_plan.N_total; k++)
    my_f_hat[k] = (double)rand() / RAND_MAX;
 
  for( j = 0; j < d*M; j++)
    my_plan.x[j] = my_x[j];

  /** init pseudo random Fourier coefficients and show them */
  for( k = 0; k < my_plan.N_total; k++)
    my_plan.f_hat[k] = my_f_hat[k];
   
  /** direct trafo and show the result */
  ndct_trafo( &my_plan);

  for( j = 0; j < M; j++)
    slow[j] = my_plan.f[j];

  nfct_finalize( &my_plan);

    
  for( m = m1; m < m2+1; m++) {

    nfct_init_m(&my_plan, d, N, M, m);
      
      /** init pseudo random nodes */
    for( j = 0; j < d * M; j++)
      my_plan.x[j] = my_x[j];
      
    /** precompute psi, the entries of the matrix B */
    if( my_plan.nfct_flags & PRE_PSI)
      nfct_precompute_psi( &my_plan);
      
    /** init pseudo random Fourier coefficients and show them */
    for( k = 0; k < my_plan.N_total; k++)
      my_plan.f_hat[k] = my_f_hat[k];
      
    /** approx. trafo and show the result */
    nfct_trafo( &my_plan);

    printf( "%d  %e, %e\n", m,
      nfft_error_l_2_double( slow, my_plan.f, M),
      nfft_error_l_infty_1_double( slow, my_plan.f, M, my_plan.f_hat, my_plan.N_total));
      
    /** finalize the one dimensional plan */
    nfct_finalize(&my_plan);
  }

  fftw_free( slow);
  fftw_free( my_x);
  fftw_free( my_f_hat);
}

int main(int argc, char **argv) { 

  accuracy( atoi( argv[1]), atoi( argv[2]), atoi( argv[3]), atoi( argv[4]));

  return EXIT_SUCCESS;
}
