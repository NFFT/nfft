#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

void simple_test_nsfft(int d, int J, int M)
{
  int K=12;
  nsfft_plan p;

  nsfft_init(&p, d, J, M, 6, NSDFT);

  nsfft_init_random_nodes_coeffs(&p);

  nfft_vpr_complex(p.f_hat, K, "frequencies, vector f_hat (first few entries)");

  /** direct trafo and show the result */
  nsdft_trafo(&p);
  nfft_vpr_complex(p.f, K, "nsdft, vector f (first few entries)"); 

  /** approx. trafo and show the result */
  nsfft_trafo(&p);
  nfft_vpr_complex(p.f, K, "nsfft, vector f (first few entries)");

  /** direct adjoint and show the result */
  nsdft_adjoint(&p);
  nfft_vpr_complex(p.f_hat, K, "adjoint nsdft, vector f_hat, (first few entries)");

  /** approx. adjoint and show the result */
  nsfft_adjoint(&p);
  nfft_vpr_complex(p.f_hat, K, "adjoint nsfft, vector f_hat, (first few entries)");

  /** finalise the one dimensional plan */
  nsfft_finalize(&p);
}

int main(int argc,char **argv)
{
  int d, J, M;

  system("clear");
  printf("1) computing a two dimensional nsdft, nsfft and adjoints\n\n");
  d=2;
  J=5;
  M=(J+4)*nfft_int_2_pow(J+1);
  simple_test_nsfft(d,J,M);
  getc(stdin);

  system("clear");
  printf("2) computing a three dimensional nsdft, nsfft and adjoints\n\n");
  d=3;
  J=5;
  M=6*nfft_int_2_pow(J)*(nfft_int_2_pow((J+1)/2+1)-1)+nfft_int_2_pow(3*(J/2+1));
  simple_test_nsfft(d,J,M);
    
  return 1;
}
