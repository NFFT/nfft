#include "utils.h"
#include "nfft.h"

double norm2(int k0, int k1)
{
  return sqrt(k0*k0+k1*k1);
}

void glacier(int N,int M)
{
  int j,k,k0,k1,l;                      /* nodes, freqencies iterations       */
  nfft_plan my_plan;                    /* plan for the one dimensional nfft  */
  infft_plan my_iplan;                  /* plan for the one dimensional infft */
  FILE* fp;                             /* input file                         */
  int my_N[2],my_n[2];
   
  /* initialise my_plan */
  my_N[0]=N; my_n[0]=next_power_of_2(N);
  my_N[1]=N; my_n[1]=next_power_of_2(N);
  nfft_init_specific(&my_plan, 2, my_N, M, my_n, 6, 
		     PRE_PHI_HUT| PRE_PSI| PRE_FULL_PSI| MALLOC_X| MALLOC_F_HAT|
		     MALLOC_F| FFT_OUT_OF_PLACE,
		     FFTW_MEASURE| FFTW_DESTROY_INPUT);
  //nfft_init_2d(&my_plan,N,N,M);

  /* initialise my_iplan, specific */
  infft_init_specific(&my_iplan,&my_plan, CGNE_R | PRECOMPUTE_DAMP);

  /* init nodes */
  fp=fopen("input_data.dat","r");
  for(j=0;j<my_plan.M;j++)
    {
      fscanf(fp,"%le %le %le",&my_plan.x[2*j+0],&my_plan.x[2*j+1],
	     &my_iplan.given_f[j][0]);
      my_iplan.given_f[j][1]=0.0;
    }
  fclose(fp);

  /* precompute psi */
  if(my_plan.nfft_flags & PRE_PSI)
    {
      nfft_precompute_psi(&my_plan);

      if(my_plan.nfft_flags & PRE_FULL_PSI)
	nfft_full_psi(&my_plan,pow(10,-10));
    }

  /* initialise damping factors */
  if(my_iplan.infft_flags & PRECOMPUTE_DAMP)
    for(k0=0;k0<my_plan.N[0];k0++)
      for(k1=0;k1<my_plan.N[1];k1++)
	{
	  my_iplan.w_hat[k0*my_plan.N[1]+k1]=
	    modified_multiquadric(1.3,0.8,norm2(k0-my_plan.N[0]/2,
						k1-my_plan.N[1]/2));
	}

  /* init some guess */
  for(k=0;k<my_plan.N_L;k++)
    {
      my_iplan.f_hat_iter[k][0]=0.0;
      my_iplan.f_hat_iter[k][1]=0.0;
    }

  /* inverse trafo */  
  infft_before_loop(&my_iplan);
  for(l=0;l<30;l++)
    { 
      fprintf(stderr,"%e,\n",sqrt(my_iplan.dot_r_iter));
      infft_loop_one_step(&my_iplan);
    }

  for(k=0;k<my_plan.N_L;k++)
    printf("%le %le\n",my_iplan.f_hat_iter[k][0],my_iplan.f_hat_iter[k][1]);

  infft_finalize(&my_iplan);  
  nfft_finalize(&my_plan);  
}



int main(int argc, char **argv)
{
  glacier(atoi(argv[1]),atoi(argv[2]));

  return 1;
}
