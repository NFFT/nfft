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

  vpr_complex(my_plan.f_hat, 8,"frequencies, vector f_hat, 0,...,7");

  /** direct trafo and show the result */
  nsdft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nsdft, vector f"); 

  /** approx. trafo and show the result */
  nsfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nsfft, vector f");

  /** direct adjoint and show the result */
  nsdft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat, 8,"adjoint nsdft, vector f_hat, 0,...,7");

  /** approx. adjoint and show the result */
  nsfft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat, 8,"adjoint nsfft, vector f_hat, 0,...,7");

  /** finalise the one dimensional plan */
  nsfft_finalize(&my_plan);
}

void accuracy_nsfft(int d, int J, int M)
{
  int j,k;                              /**< index for nodes and freqencies  */
  nsfft_plan my_plan;                   /**< plan for the nfft               */
  complex *swap_sndft_trafo, *swap_sndft_adjoint;

  nsfft_init(&my_plan, d, J, M, 6, SNDFT);

  swap_sndft_trafo=(complex*) fftw_malloc(my_plan.M_total*sizeof(complex));
  swap_sndft_adjoint=(complex*) fftw_malloc(my_plan.N_total*sizeof(complex));

  nsfft_init_random_nodes_coeffs(&my_plan);

  /** direct trafo */
  nsdft_trafo(&my_plan);
  
  SWAP_complex(swap_sndft_trafo,my_plan.f);

  /** approx. trafo */
  nsfft_trafo(&my_plan);
  
  printf("%5d\t %+.5E\t",J, 
         error_l_infty_1_complex(swap_sndft_trafo, my_plan.f, my_plan.M_total,
                                 my_plan.f_hat, my_plan.N_total));
  fflush(stdout);

  vrand_unit_complex(my_plan.f, my_plan.M_total);

  /** direct adjoint */
  nsdft_adjoint(&my_plan);
  
  SWAP_complex(swap_sndft_adjoint,my_plan.f_hat);

  /** approx. adjoint */
  nsfft_adjoint(&my_plan);
  
  printf("%+.5E\n", 
         error_l_infty_1_complex(swap_sndft_adjoint, my_plan.f_hat,
                                 my_plan.N_total,
                                 my_plan.f, my_plan.M_total));
  fflush(stdout);

  fftw_free(swap_sndft_adjoint);
  fftw_free(swap_sndft_trafo);

  /** finalise the one dimensional plan */
  nsfft_finalize(&my_plan);
}

int main()
{ 
  int J;

  system("clear");
  printf("1) computing a two dimensional nsdft, nsfft\n\n");
  simple_test_nsfft(2,5,8);

  printf("\npress a key"); getc(stdin);

  system("clear");
  printf("2) computing a three dimensional nsdft, nsfft\n\n");
  simple_test_nsfft(3,5,8);

  printf("\npress a key"); getc(stdin);

  system("clear");
  printf("3) accuracy, two dimensional nsdft, nsfft\n\n");
  for(J=1; J<10; J++)
    accuracy_nsfft(2, J, 1000);

  printf("\npress a key"); getc(stdin);

  system("clear");
  printf("4) accuracy, three dimensional nsdft, nsfft\n\n");
  for(J=1; J<9; J++)
    accuracy_nsfft(3, J, 1000);

  return 1;
}
