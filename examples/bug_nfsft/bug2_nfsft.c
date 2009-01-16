/* Include standard C headers. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

/* Include NFFT 3 utilities headers. */
#include "util.h"
/* Include NFFT3 library header. */
#include "nfft3.h"
#include "infft.h"

static void simple_test_nfsft(void)
{
  int j;                       /**< Index for nodes                                 */
  int k;                       /**< Index for freqency degree                       */
  int n;                       /**< Index for freqency degree                       */
  nfsft_plan plan;             /**< Plan for the nfft                               */
  const int N = 50;             /**< The bandwidth M                                 */
  const int M = 100;             /**< The number of nodes M                           */
  const double THRESHOLD = 1000.0; /**< The threshold for the NFSFT stabilization
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
    FFT_OUT_OF_PLACE, 9);

  /* Init pseudo random nodes. */
  for (j = 0; j < plan.M_total; j++)
  {
    plan.x[2*j]=RAND-0.5;
    plan.x[2*j+1]=0.5*RAND;
  }

  /* Do precomputation for nodes. */
  nfsft_precompute_x(&plan);

  /* Init pseudo random function values */
  for(j=0;j<M;j++)
  {
    plan.f[j]=RAND;
  }

  /* Compute approximate adjoint transformation */
  nfsft_adjoint(&plan);

  double * diff_f_hat_real = malloc((N+1)*(N+1)*sizeof(double));
  double * diff_f_hat_imag = malloc((N+1)*(N+1)*sizeof(double));

  j=0;
   for (k = 0; k <= N; k++)
   {
     for (n = -k; n <= k; n++)
     {
	  diff_f_hat_real[j] = creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]);
	  diff_f_hat_imag[j] = creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]);
          j++;
     }
   }
  
  /* Compute direct adjoint transformation */
  ndsft_adjoint(&plan);
  
  /* Compute the difference of both */

  j=0;
  for (k = 0; k <= N; k++)
  {
    for (n = -k; n <= k; n++)
    {
	  diff_f_hat_real[j] = diff_f_hat_real[j] - creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]);
	  diff_f_hat_imag[j] = diff_f_hat_imag[j] - creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]);
       /* Show the coefficients, which differ too much */
       if(fabs(diff_f_hat_real[j])+fabs(diff_f_hat_imag[j]) > 1e-5)
         fprintf(stdout,"diff_f_hat[%d,%d] = %le + I*%le\n",k,n, diff_f_hat_real[j],diff_f_hat_imag[j]);
       j++;
    }
  }

  /* Finalise the plan. */
  nfsft_finalize(&plan);
	nfsft_forget();
}

int main(void)
{
  system("clear");
  printf("1) computing a ndsft, a nfsft, an adjoint ndsft, and an adjoint nfsft\n\n");
  simple_test_nfsft();

  return EXIT_SUCCESS;
}
