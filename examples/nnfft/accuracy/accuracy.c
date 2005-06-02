#include "utils.h"
#include "nfft.h"

void accuracy(int d)
{
  int j,k,t,m;                         
  nnfft_plan my_plan;                   
  complex *slow;

  int N[d],n[d];
  int M_total,N_total;
  M_total=10000;N_total=1;

  slow=(complex*)fftw_malloc(M_total*sizeof(complex));

  for(t=0; t<d; t++)
    {
      N[t]=(1<<(12/d));
      n[t]=2*N[t];
      N_total*=N[t];
    }

  /** init a plan */
  for(m=0; m<15; m++)
    {
      nnfft_init_specific(&my_plan, d, N_total, M_total, N, n, m, 
                          PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F );

      
      /** init pseudo random nodes */
      for(j=0; j<my_plan.M_total; j++)
	for(t=0; t<d; t++)
	  {
	    my_plan.x[d*j+t]=((double)rand())/RAND_MAX-0.5;
	  }
          
     for(j=0; j<my_plan.N_total; j++)
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
      for(k=0;k<my_plan.N_total;k++)
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
	     E_2_error_c(slow, my_plan.f, M_total),
	     E_infty_error_c(slow, my_plan.f, M_total, my_plan.f_hat, my_plan.N_total));
      
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
