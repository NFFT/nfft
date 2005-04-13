#include "utils.h"
#include "nfft.h"

void timing_1d(int setup)
{
  int j,k;
  nfft_plan my_plan;
  int my_N,my_n;
  double t;
  
  for(my_N=1024; my_N<=(1<<20); my_N*=2)
    {
      my_n=next_power_of_2(2*my_N);

      switch(setup)
	{
	case 1:
	  {
	    nfft_init_specific(&my_plan, 1, &my_N, my_N, &my_n, 6,
			       MALLOC_X| MALLOC_F_HAT| MALLOC_F, 
			       FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
	    break;
	  }
	case 2:
	  {
	    nfft_init_1d(&my_plan,my_N,my_N);
	    break;
	  }
	case 3:
	  {
	    nfft_init_specific(&my_plan, 1, &my_N, my_N, &my_n, 6,
			       PRE_PHI_HUT| PRE_PSI| PRE_FULL_PSI| MALLOC_X|
			       MALLOC_F_HAT| MALLOC_F| FFT_OUT_OF_PLACE, 
			       FFTW_MEASURE | FFTW_DESTROY_INPUT);
	  }
	} 

      /**/

      for(j=0;j<my_plan.M;j++)
	{
	  my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
	}

      if(my_plan.nfft_flags & PRE_PSI)
	{
	  nfft_precompute_psi(&my_plan);
	  
	  if(my_plan.nfft_flags & PRE_FULL_PSI)
	    nfft_full_psi(&my_plan,pow(10,-20));
	}

      for(k=0;k<my_plan.N_L;k++)
	{
	  my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
	  my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
	}

      printf("%d\t",my_N);

      t=second();
      nfft_trafo(&my_plan);
      t=second()-t;
      printf("%e\t",t);

      if((my_N<=16384)&&(setup==1))
	{
	  t=second();
	  ndft_trafo(&my_plan);
	  t=second()-t;
	  printf("%e\n",t);
	}
      else
	printf("inf\n");

      fflush(stdout);

      nfft_finalize(&my_plan);  
    }
} 


void timing_2d()
{
  int j,k;
  nfft_plan my_plan;
  int my_N;
  double t;

  for(my_N=16; my_N<=1024; my_N*=2)
    {
      nfft_init_2d(&my_plan,my_N,my_N,my_N*my_N);

      for(j=0;j<my_plan.M;j++)
	{
	  my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
	  my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
	}

      if(my_plan.nfft_flags & PRE_PSI)
	nfft_precompute_psi(&my_plan);

      for(k=0;k<my_plan.N_L;k++)
	{
	  my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
	  my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
	}

      printf("%d\t",my_N);

      t=second();
      nfft_trafo(&my_plan);
      t=second()-t;
      printf("%e\t",t);

      if(my_N<=128)
	{
	  t=second();
	  ndft_trafo(&my_plan);
	  t=second()-t;
	  printf("%e\n",t);
	}
      else
	printf("inf\n");
      
      fflush(stdout);

      nfft_finalize(&my_plan);  
    }
} 

void timing_2d_fast()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N;

  int N[2],n[2];

#ifndef MEASURE_TIME
  fprintf(stderr,"\nWarning: MEASURE_TIME not enabled\n");
#endif


  for(my_N=128; my_N<=1024; my_N+=128)
    {
      N[0]=my_N; n[0]=next_power_of_2(2*my_N);
      N[1]=my_N; n[1]=next_power_of_2(2*my_N);

      /*nfft_init_2d(&my_plan,my_N,my_N,my_N*my_N);*/
      nfft_init_specific(&my_plan, 2, N, my_N*my_N, n, 6,
			 PRE_PHI_HUT| PRE_PSI| PRE_FULL_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F| FFT_OUT_OF_PLACE, 
			 FFTW_MEASURE | FFTW_DESTROY_INPUT);

      for(j=0;j<my_plan.M;j++)
	{
	  my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
	  my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
	}

      if(my_plan.nfft_flags & PRE_PSI)
	{
	  nfft_precompute_psi(&my_plan);
	  
	  if(my_plan.nfft_flags & PRE_FULL_PSI)
	    nfft_full_psi(&my_plan,pow(10,-20));
	}

      for(k=0;k<my_plan.N_L;k++)
	{
	  my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
	  my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
	}

      printf("%d\t",my_N);

      nfft_trafo(&my_plan);

      printf("\n");
      
      fflush(stdout);

      nfft_finalize(&my_plan);  
    }
} 

void timing_3d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N;
  double t;

  for(my_N=16; my_N<=64; my_N+=4)
    {
      nfft_init_3d(&my_plan,my_N,my_N,my_N,my_N*my_N*my_N);

      for(j=0;j<my_plan.M;j++)
	{
	  my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
	  my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
	  my_plan.x[2*j+2]=((double)rand())/RAND_MAX-0.5;
	}

      if(my_plan.nfft_flags & PRE_PSI)
	nfft_precompute_psi(&my_plan);

      for(k=0;k<my_plan.N_L;k++)
	{
	  my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
	  my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
	}

      printf("%d\t",my_N);

      t=second();
      nfft_trafo(&my_plan);
      t=second()-t;
      printf("%e\t",t);

      if(my_N<=32)
	{
	  t=second();
	  ndft_trafo(&my_plan);
	  t=second()-t;
	  printf("%e\n",t);
	}
      else
	printf("inf\n");

      fflush(stdout);

      nfft_finalize(&my_plan);  
    }
} 


int main()
{ 
  /*timing_1d(1);
    timing_1d(2);
    timing_1d(3);*/

  /*timing_2d();*/
  
  /*timing_3d();*/

  timing_2d_fast();

  return 1;
}
