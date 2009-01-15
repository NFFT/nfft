#include "nfft3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415926535897932384626433832795029L

void nodesTetraeder(double ** theta, double ** phi, int * M)
{
  *M = 4;
  
  double * theta_temp = malloc(4*sizeof(double));
  double * phi_temp = malloc(4*sizeof(double));
  
  theta_temp[0]=acos(sqrt(1./3));
  phi_temp[0]=0;
  theta_temp[1]=acos(sqrt(1./3));
  phi_temp[1]=PI;
  theta_temp[2]=acos(-sqrt(1./3));
  phi_temp[2]=PI/2;
  theta_temp[3]=acos(-sqrt(1./3));
  phi_temp[3]=3./2*PI;
  
  *theta = theta_temp;
  *phi = phi_temp;
}


int main(int argc, char *argv[])
{
  int M; //number of nodes
  int N = 19;//polynomial degree
  const double THRESHOLD = 1000.0;

  nfsft_plan plan;
  
  double * theta = NULL; //array for theta-coords
  double * phi = NULL;  //array for phi-coords

  nodesTetraeder(&theta, &phi, &M);

fprintf(stdout,"M=%d, N=%d, phi=%p, theta=%p\n",M,N,(void *) phi,(void *) theta);
//int opop; for(opop=0;opop<M;opop++) fprintf(stdout,"theta_temp = %f\n",phi[opop]);

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

fprintf(stdout,"M=%d, N=%d, phi=%p, theta=%p\n",M,N,(void *) phi,(void *) theta);
//for(opop=0;opop<M;opop++) fprintf(stdout,"theta_temp = %f\n",phi[opop]);

nfsft_finalize(&plan);
nfsft_forget();
  
free(theta);free(phi);

  return EXIT_SUCCESS;
}
