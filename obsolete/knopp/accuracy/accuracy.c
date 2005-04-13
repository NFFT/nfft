#include "utils.h"
#include "nfft.h"

void accuracy(int d)
{
  int j,k,t,m;                         
  nnfft_plan my_plan;                   
  fftw_complex *slow;

  int N[d],n[d];
  int M1,M2;
  M1=10000;M2=1;

  slow=(fftw_complex*)fftw_malloc(M1*sizeof(fftw_complex));

  for(t=0; t<d; t++)
    {
      N[t]=(1<<(12/d));
      n[t]=2*N[t];
      M2*=N[t];
    }

  /** init a plan */
  for(m=0; m<15; m++)
    {
      nnfft_init_specific(&my_plan, d, M1, M2, N, n, m, 
                          PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F );

      
      /** init pseudo random nodes */
      for(j=0; j<my_plan.M1; j++)
	for(t=0; t<d; t++)
	  {
	    my_plan.x[d*j+t]=((double)rand())/RAND_MAX-0.5;
	  }
          
     for(j=0; j<my_plan.M2; j++)
        for(t=0; t<d; t++)
          {
            my_plan.v[d*j+t]=((double)rand())/RAND_MAX-0.5;
          }
      
      /** precompute psi, the entries of the matrix B */
      if(my_plan.nnfft_flags & PRE_PSI)
	nnfft_precompute_psi(&my_plan);
        
       /** precompute psi, the entries of the matrix D */ 
      if(my_plan.nnfft_flags & PRE_PHI_HUT)
        nnfft_precompute_phi_hut(&my_plan);
      
      /** init pseudo random Fourier coefficients */
      for(k=0;k<my_plan.M2;k++)
	{
	  my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
	  my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
	}
      
      /** direct trafo and show the result */
      nndft_trafo(&my_plan);
      
      SWAPC(my_plan.f,slow);
      
      /** approx. trafo and show the result */
      nnfft_trafo(&my_plan);
        
      printf("%e, %e\n",
	     E_2_error_c(slow, my_plan.f, M1),
	     E_infty_error_c(slow, my_plan.f, M1, my_plan.f_hat, my_plan.M2));
      
      /** finalise the one dimensional plan */
      nnfft_finalize(&my_plan);
    }
}

int main()
{ 
  int d;
  for(d=1; d<4; d++)
    accuracy(d);
 
  return 1;
}
