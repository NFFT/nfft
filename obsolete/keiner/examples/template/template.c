/**
 * Program template for using the NFSFT library.
 */

/* Include NFSFT header */
#include<nfsft.h>
#include<complex.h>
#include<stdlib.h>
#include<stdio.h>

/**
 * The main program
 */
int main()
{
  /** The bandwidth */
  const int M = 1<<3;
  /** Next greater power of two with respect to M */
  const int N = 1<<ngpt(M);
  /** 
   * For the Fourier coefficients. 
   * Layout: f_hat[M+n][k] = a_k^n for n=-M,...,M; k=0,...,M
   */
  complex **f_hat;
  /** The number of nodes */
  const int D = 1<<2;
  /** 
   * The nodes on S^2 
   * Layout: angles[2d] = phi_d, angles[2d+1] = theta_d, for d=0,...,D-1
   */
  double *angles;
  /** Function values */
  complex *f;
  
  /* NFSFT transform plans */
  nfsft_plan plan, plan_adjoint;
  
  /** Loop counters */
  int n,k,d;
  
  /* Allocate memory */
  f_hat = (complex**) malloc((2*M+1)*sizeof(complex*));
  for (n = -M; n <= M; n++)
  {
    /* Length must be a power of two + 1 >= M */
    f_hat[n+M] = (complex*) malloc((N+1)*sizeof(complex));
  }  
  angles = (double*) malloc(2*D*sizeof(double));
  f = (complex*) malloc(D*sizeof(complex));
  
  /* Generate random Fourier coeffcients. */
  for (n=-M; n <= M; n++)
  {
    for (k=abs(n); k <= M; k++)
    {
      f_hat[n+M][k] = drand48() + I * drand48();
    }  
  }  
  
  /* Generate random nodes */
  for (d = 0; d < D; d++)
  {
    /* phi \in [-1/2,1/2) */
    angles[2*d] = drand48()-0.5;
    /* theta \in [0,1/2] */
    angles[2*d+1] = 0.5*drand48();
  }  
  
  /* Compute NFSFT */
  plan = nfsft_init(D, M, angles, f_hat, f, 0U);
  nfsft_trafo(plan);
  nfsft_finalize(plan);
  
  /* Print result */
  for (d = 0; d < D; d++)
  {
    fprintf(stdout,"f(%.16f,%.16f) = %.16f + %.16fi\n",angles[2*d],
            angles[2*d+1],creal(f[d]),cimag(f[d]));
  }  
  
  /* Compute adjoint NFSFT */
  plan_adjoint = nfsft_init(D, M, angles, f_hat, f, 0U);
  nfsft_adjoint(plan_adjoint);
  nfsft_finalize(plan_adjoint);
  
  /* Print result */
  for (k= 0; k <= M; k++)
  {
    for (n=-k; n <= k; n++)
    {
      fprintf(stdout,"a_%d^%d = %.16f + %.16fi\n",k,n,creal(f_hat[n+M][k]),
              cimag(f_hat[n+M][k]));
    }  
  }   
  
  /* Free memory */
  free(f);
  free(angles);
  for (n = -M; n <= M; n++)
  {
    free(f_hat[n+M]);
  }  
  free(f_hat);
  
  return EXIT_SUCCESS;
}