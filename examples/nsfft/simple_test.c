#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

void simple_test_nsfft(int d, int J, int M)
{
  int j,k;                              /**< index for nodes and freqencies  */
  nsfft_plan my_plan;                   /**< plan for the nfft               */
  complex *swap_sndft;

  nsfft_init(&my_plan, d, J, M, 6, SNDFT);

  swap_sndft=(complex*) fftw_malloc(M*sizeof(complex));

  nsfft_init_random_nodes_coeffs(&my_plan);

  vpr_complex(my_plan.f_hat, 16,"frequencies, vector f_hat, 0,...,16");

  /** direct trafo and show the result */
  nsdft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nsdft, vector f"); 

  /** approx. trafo and show the result */
  nsfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nsfft, vector f");

  /** finalise the one dimensional plan */
  nsfft_finalize(&my_plan);
}

int main()
{ 
  system("clear");
  printf("1) computing a two dimensional nfft, nsdft, nsfft\n\n");
  simple_test_nsfft(2,5,8);

  getc(stdin);

  system("clear");
  printf("1) computing a three dimensional nfft, nsdft, nsfft\n\n");
  simple_test_nsfft(3,10,8);

  return 1;
}
