/* Include standard C headers. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* Include NFFT3 library header. */
#include "nfft3.h"

/* Include NFFT 3 utilities headers. */
#include "util.h"

void simple_test_nfsft()
{
  int j;                      /**< Index for nodes                                 */
  int k;                      /**< Index for freqency degree                       */
  int n;                      /**< Index for freqency degree                       */
  nfsft_plan plan;            /**< Plan for the nfft                               */
  const int N = 8;            /**< The bandwidth M                                 */
  const int M = 8;            /**< The number of nodes M                           */
  const int THRESHOLD = 1000; /**< The threshold for the NFSFT stabilization
                                   procedure.                                      */

  /* Precompute. */
  nfsft_precompute(N,THRESHOLD,0U,0U);

  /* Init a transform plan using the guru interface. All arrays for input and
   * output variables are allocated by nfsft_init_guru(). Computations are
   * performed with respect to L^2-normalized spherical harmonics Y_k^n. The
   * array of spherical Fourier coefficients is preserved during
   * transformations. The internal NFFT uses a cut-off parameter of 6.
   */
  nfsft_init_guru(&plan, N, M, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
    NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_PRESERVE_F_HAT,
    ((N>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT |
    FFT_OUT_OF_PLACE, 6);

  /* Init pseudo random nodes. */
  for (j = 0; j < plan.M_total; j++)
  {
    plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
    plan.x[2*j+1]=0.5*((double)rand())/RAND_MAX;
  }

  /* Do precomputation for nodes. */
  nfsft_precompute_x(&plan);

  /* Init pseudo random Fourier coefficients. */
  for (k = 0; k <= plan.N; k++)
  {
    for (n = -k; n <= k; n++)
    {
      plan.f_hat[NFSFT_INDEX(k,n,&plan)] =
        ((double)rand())/RAND_MAX - 0.5 +  I*(((double)rand())/RAND_MAX - 0.5);
    }
  }

  /*for (k = 0; k < plan.N_total; k++)
  {
     fprintf(stderr,"f_hat[%d] = %le +I*%le\n",k,creal(plan.f_hat[k]),
       cimag(plan.f_hat[k]));
  }*/

  //nfft_vpr_complex(plan.f_hat,plan.N_total,"given Fourier coefficients, vector f_hat");

  /* Compute direct transformation and display the result. */
  ndsft_trafo(&plan);
  nfft_vpr_complex(plan.f,plan.M_total,"ndsft, vector f");

  /* Compute approximate transformation and display the result. */
  nfsft_trafo(&plan);
  nfft_vpr_complex(plan.f, plan.M_total,"nfsft, vector f");
  /*for (k = 0; k < plan.N_total; k++)
  {
     fprintf(stderr,"f_hat[%d] = %le +I*%le\n",k,creal(plan.f_hat[k]),
       cimag(plan.f_hat[k]));
  }*/

  /* Compute direct adjoint transformation and display the result. */
  nfsft_adjoint(&plan);
  for (k = 0; k <= plan.N; k++)
  {
    for (n = -k; n <= k; n++)
    {
      fprintf(stdout,"f_hat[%d,%d] = %le + I*%le\n",k,n,
        creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]),
        cimag(plan.f_hat[NFSFT_INDEX(k,n,&plan)]));
    }
  }
  //nfft_vpr_complex(plan.f_hat,plan.N_total,"adjoint ndsft, vector f_hat");

  /* COmpute approximate adjoint transformation and display the result */
  ndsft_adjoint(&plan);
  for (k = 0; k <= plan.N; k++)
  {
    for (n = -k; n <= k; n++)
    {
      fprintf(stdout,"f_hat[%d,%d] = %le + I*%le\n",k,n,
        creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]),
        cimag(plan.f_hat[NFSFT_INDEX(k,n,&plan)]));
    }
  }
  //nfft_vpr_complex(plan.f_hat,plan.N_total,"adjoint nfsft, vector f_hat");

  /* Finalise the plan. */
  nfsft_finalize(&plan);
	nfsft_forget();
}

int main()
{
  int l,m;

  system("clear");
  printf("1) computing a ndsft, a nfsft, an adjoint ndsft, and an adjoint nfsft\n\n");
  simple_test_nfsft();

  return EXIT_SUCCESS;
}
