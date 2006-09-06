#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

void simple_test_nfft_1d()
{
  int N[1],n[1];
  nfft_plan p;

  N[0]=14; n[0]=32;

  /** init an one dimensional plan */
  nfft_init_guru(&p, 1, N, 19, n, 12,
                 PRE_PHI_HUT| PRE_LIN_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  vrand_shifted_unit_double(p.x,p.M_total);

  /** precompute psi, the entries of the matrix B */
  if(p.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);

  /** init pseudo random Fourier coefficients and show them */
  vrand_unit_complex(p.f_hat,p.N_total);
  vpr_complex(p.f_hat,p.N_total,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  ndft_trafo(&p);
  vpr_complex(p.f,p.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&p);
  vpr_complex(p.f,p.M_total,"nfft, vector f");
  
  /** approx. adjoint and show the result */
  ndft_adjoint(&p);
  vpr_complex(p.f_hat,p.N_total,"adjoint ndft, vector f_hat");

  /** approx. adjoint and show the result */
  nfft_adjoint(&p);
  vpr_complex(p.f_hat,p.N_total,"adjoint nfft, vector f_hat");

  /** finalise the one dimensional plan */
  nfft_finalize(&p);
}

void simple_test_nfft_2d()
{
  int K,N[2],n[2];
  double t;

  nfft_plan p;
   
  N[0]=70; n[0]=128;
  N[1]=50; n[1]=128;
  K=12;

  t=second();
  /** init a two dimensional plan */
  nfft_init_guru(&p, 2, N, N[0]*N[1], n, 4,
		 PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT| MALLOC_X| MALLOC_F |
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  vrand_shifted_unit_double(p.x,p.d*p.M_total);

  /** precompute psi, the entries of the matrix B */
  if(p.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&p);

  /** init pseudo random Fourier coefficients and show them */
  vrand_unit_complex(p.f_hat,p.N_total);

  t=second()-t;
  vpr_complex(p.f_hat,K,
              "given Fourier coefficients, vector f_hat (first few entries)");
  printf(" ... initialisation took %e seconds.\n",t);

  /** direct trafo and show the result */
  t=second();
  ndft_trafo(&p);
  t=second()-t;
  vpr_complex(p.f,K,"ndft, vector f (first few entries)");
  printf(" took %e seconds.\n",t);

  /** approx. trafo and show the result */
  t=second();
  nfft_trafo(&p);
  t=second()-t;
  vpr_complex(p.f,K,"nfft, vector f (first few entries)");
  printf(" took %e seconds.\n",t);

  /** direct adjoint and show the result */
  t=second();
  ndft_adjoint(&p);
  t=second()-t;
  vpr_complex(p.f_hat,K,"adjoint ndft, vector f_hat (first few entries)");
  printf(" took %e seconds.\n",t);

  /** approx. adjoint and show the result */
  t=second();
  nfft_adjoint(&p);
  t=second()-t;
  vpr_complex(p.f_hat,K,"adjoint nfft, vector f_hat (first few entries)"); 
  printf(" took %e seconds.\n",t);

  /** finalise the two dimensional plan */
  nfft_finalize(&p);
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

  return 1;
}
