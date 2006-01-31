/* Include standard C headers. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* Include NFFT3 library header. */
#include "util.h"

/* Include NFFT 3 utilities headers. */
#include "nfft3.h"

void simple_test_nfsft()
{
  int j;                      /**< Index for nodes                            */
  int k;                      /**< Index for freqency degree                  */
  int n;                      /**< Index for freqency degree                  */
  nfsft_plan plan;            /**< Plan for the nfft                          */
  const int N = 8;            /**< The bandwidth M                            */
  const int M = 8;            /**< The number of nodes M                      */
  const int THRESHOLD = 1000; /**< The threshold for the NFSFT stabilization 
                                   procedure.                                 */

  /* Init a transform plan using the guru interface. All arrays for input and 
   * output variables are allocated by nfsft_init_guru(). Computations are 
   * performed with respect to L^2-normalized spherical harmonics Y_k^n. The 
   * array of spherical Fourier coefficients is preserved during 
   * transformations. The internal NFFT uses a cut-off parameter of 6. 
   */
  nfsft_init_guru(&plan, N, M, NFSFT_MALLOC_X | NFSFT_MALLOC_F | 
    NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_PRESERVE_F_HAT |
    NFSFT_PRESERVE_F | NFSFT_PRESERVE_X | NFSFT_USE_DPT, 6);
    
  /* Init pseudo random nodes. */
  for (j = 0; j < plan.M_total; j++)
  {
    plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
    plan.x[2*j+1]=0.5*((double)rand())/RAND_MAX;
  }

  /* Precompute. */
  nfsft_precompute(N,THRESHOLD,0U);

  /* Init pseudo random Fourier coefficients and display them. */
  for (k = 0; k <= plan.N; k++)
  {
    for (n = -k; n <= k; n++)
    {
      plan.f_hat[NFSFT_INDEX(k,n,&plan)] = 
        ((double)rand())/RAND_MAX - 0.5 +  I*(((double)rand())/RAND_MAX - 0.5);
    }
  }

  vpr_complex(plan.f_hat,plan.N_total,"given Fourier coefficients, vector f_hat"); 

  /* Compute direct transformation and display the result. */
  ndsft_trafo(&plan);
  vpr_complex(plan.f,plan.M_total,"ndsft, vector f"); 

  /* Compute approximate transformation and display the result. */
  nfsft_trafo(&plan);
  vpr_complex(plan.f,plan.M_total,"nfsft, vector f");

  /* Compute direct adjoint transformation and display the result. */
  ndsft_adjoint(&plan);
  vpr_complex(plan.f_hat,plan.N_total,"adjoint ndsft, vector f_hat");

  /* COmpute approximate adjoint transformation and display the result */
  nfsft_adjoint(&plan);
  vpr_complex(plan.f_hat,plan.N_total,"adjoint nfsft, vector f_hat");

  /* Finalise the plan. */
  nfsft_finalize(&plan);
}

int main()
{
  int l,m;

  system("clear");
  printf("1) computing a ndsft, a nfsft, an adjoint ndsft, and an adjoint nfsft\n\n");
  simple_test_nfsft();

  /*getc(stdin);

  system("clear"); 
  printf("2) computing a two dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_2d();*/
/*
  getc(stdin);

  system("clear"); 
  printf("3) computing an one dimensional infft\n\n");
  simple_test_infft_1d();

  getc(stdin);
    
  system("clear"); 
  printf("3) computing times for one dimensional nfft\n\n");
  measure_time_nfft_1d();
*/
  return EXIT_SUCCESS;
}
