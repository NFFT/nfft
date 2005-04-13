#include "utils.h"
#include "nfft.h"

void accuracy(int d)
{
  int j,k,t,m;                         
  nfft_plan my_plan;                   
  fftw_complex *slow;

  int N[d],n[d];
  int M;
  M=10000;

  slow=(fftw_complex*)fftw_malloc(M*sizeof(fftw_complex));

  for(t=0; t<d; t++)
    {
      N[t]=(1<<(12/d));
      n[t]=2*N[t];
    }

  /** init a plan */
  for(m=0; m<15; m++)
    {
      nfft_init_specific(&my_plan, d, N, M, n, m,
			 PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F, 
			 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
      
      /** init pseudo random nodes */
      for(j=0; j<my_plan.M; j++)
	for(t=0; t<d; t++)
	  {
	    my_plan.x[d*j+t]=((double)rand())/RAND_MAX-0.5;
	  }
      
      /** precompute psi, the entries of the matrix B */
      if(my_plan.nfft_flags & PRE_PSI)
	nfft_precompute_psi(&my_plan);
      
      /** init pseudo random Fourier coefficients and show them */
      for(k=0;k<my_plan.N_L;k++)
	{
	  my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
	  my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
	}
      
      /** direct trafo and show the result */
      ndft_trafo(&my_plan);
      
      SWAPC(my_plan.f,slow);
      
      /** approx. trafo and show the result */
      nfft_trafo(&my_plan);
      
      printf("%e, %e\n",
	     E_2_error_c(slow, my_plan.f, M),
	     E_infty_error_c(slow, my_plan.f, M, my_plan.f_hat, my_plan.N_L));
      
      /** finalise the one dimensional plan */
      nfft_finalize(&my_plan);
    }
}

int main()
{ 
  int d;
  for(d=1; d<4; d++)
    accuracy(d);
 
  return 1;
}
