#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

void simple_test_nfct_1d()
{
  int j,k;
  nfct_plan p;

  int N=14;
  int M=19;

  /** init an one dimensional plan */
  nfct_init_1d(&p,N,M);

  /** init pseudo random nodes */
  for(j = 0; j < p.d*p.M_total; j++)
    p.x[j] = 0.5 * ((double)rand()) / RAND_MAX;

  /** precompute psi, the entries of the matrix B */
  if( p.nfct_flags & PRE_PSI)
    nfct_precompute_psi( &p);

  /** init pseudo random Fourier coefficients and show them */
  for(k = 0; k < p.N_total; k++)
    p.f_hat[k] = (double)rand() / RAND_MAX;

  nfft_vpr_double(p.f_hat,p.N_total,"given Fourier coefficients, vector f_hat");

  /** direct trafo and show the result */
  ndct_trafo(&p);
  nfft_vpr_double(p.f,p.M_total,"ndct, vector f"); 

  /** approx. trafo and show the result */
  nfct_trafo(&p);
  nfft_vpr_double(p.f,p.M_total,"nfct, vector f");
  
  /** approx. adjoint and show the result */
  ndct_adjoint(&p);
  nfft_vpr_double(p.f_hat,p.N_total,"adjoint ndct, vector f_hat");

  /** approx. adjoint and show the result */
  nfct_adjoint(&p);
  nfft_vpr_double(p.f_hat,p.N_total,"adjoint nfct, vector f_hat");

  /** finalise the one dimensional plan */
  nfct_finalize(&p);
}

int main()
{
  system("clear");
  printf("computing an one dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfct_1d();
  printf("\n\n");

  return 1;
}
