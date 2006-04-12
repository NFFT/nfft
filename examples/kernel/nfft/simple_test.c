#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
#include "nfft3.h"

void simple_test_nfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies  */
  nfft_plan my_plan;                    /**< plan for the nfft               */

  int N[1],n[1];
  N[0]=14; n[0]=32;

  /** init an one dimensional plan */
 nfft_init_guru(&my_plan, 1, N, 19, n, 6,
                 PRE_PHI_HUT| PRE_LIN_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
    }

 /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,my_plan.N_total,"given Fourier coefficients, vector f_hat"); 

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** approx. adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,my_plan.N_total,"adjoint ndft, vector f_hat");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,my_plan.N_total,"adjoint nfft, vector f_hat");

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}

void simple_test_nfft_2d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */

  int N[2],n[2];
  N[0]=10; n[0]=32;
  N[1]=10; n[1]=64;

  /** init an one dimensional plan */
  //nfft_init_2d(&my_plan,12,12,19);
  nfft_init_guru(&my_plan, 2, N, 19, n, 4,
		 PRE_PHI_HUT| PRE_LIN_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT, FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
      my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
    }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,12,
	"given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** direct adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint ndft, vector f_hat (first 12 entries)");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint nfft, vector f_hat (first 12 entries)"); 

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}

void simple_test_nfft_2d_huge()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */

  int N[2],n[2];
  N[0]=2050; n[0]=4096;
  N[1]=2050; n[1]=4096;

  /** init an one dimensional plan */

  nfft_init_guru(&my_plan, 2, N, 10, n, 
		 6, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT| MALLOC_X| MALLOC_F |
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[2*j]=((double)rand())/RAND_MAX-0.5;
      my_plan.x[2*j+1]=((double)rand())/RAND_MAX-0.5;
    }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/RAND_MAX + I* ((double)rand())/RAND_MAX;

  vpr_complex(my_plan.f_hat,12,
	"given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  ndft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f"); 

  /** approx. trafo and show the result */
  nfft_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** direct adjoint and show the result */
  ndft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint ndft, vector f_hat (first 12 entries)");

  /** approx. adjoint and show the result */
  nfft_adjoint(&my_plan);
  vpr_complex(my_plan.f_hat,12,"adjoint nfft, vector f_hat (first 12 entries)"); 

  /** finalise the one dimensional plan */
  nfft_finalize(&my_plan);
}



void measure_time_nfft_1d()
{
  int j,k;                              /**< index for nodes and freqencies   */
  nfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N;
  double t;

  for(my_N=16; my_N<=1024*1024; my_N*=2)
    {
      nfft_init_1d(&my_plan,my_N,my_N);

      for(j=0;j<my_plan.M_total;j++)
	my_plan.x[j]=((double)rand())/RAND_MAX-0.5;

      if(my_plan.nfft_flags & PRE_ONE_PSI)
        nfft_precompute_one_psi(&my_plan);
  
      for(k=0;k<my_plan.N_total;k++)
	my_plan.f_hat[k]=((double)rand())/RAND_MAX + I*((double)rand())/RAND_MAX;

      if(my_N<=8192)
	{
	  t=second();
	  ndft_trafo(&my_plan);
	  t=second()-t;
	  printf("t_ndft=%e,\t",t);
	}
      else
	printf("t_ndft=nan\t");

      t=second();
      nfft_trafo(&my_plan);
      t=second()-t;
      printf("t_nfft=%e\t",t);
      
      printf("(N=M_total=%d)\n",my_N);

      nfft_finalize(&my_plan);  
    }
} 

void simple_test_infft_1d()
{
  int j,k,l;                            /**< index for nodes, freqencies,iter*/
  nfft_plan my_plan;                    /**< plan for the nfft               */
  infft_plan my_iplan;                  /**< plan for the inverse nfft       */
 
  /** initialise an one dimensional plan */
  nfft_init_1d(&my_plan, 8, 4);

  /** initialise my_iplan */
  infft_init(&my_iplan,&my_plan);

  /** init pseudo random nodes */
  vrand_shifted_unit_double(my_plan.x,my_plan.M_total);
  
  /** precompute psi, the entries of the matrix B */
  if(my_plan.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&my_plan);

  /** init pseudo random samples and show them */
  vrand_unit_complex(my_iplan.y,my_plan.M_total);

  vpr_complex(my_iplan.y,my_plan.M_total,"given data, vector y");

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_plan.N_total;k++)
      my_iplan.f_hat_iter[k]=0; 

  vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
	"approximate solution, vector f_hat_iter");

  /** solve the system */
  infft_before_loop(&my_iplan);
  for(l=0;l<3;l++)
    {
      printf("iteration l=%d\n",l);
      infft_loop_one_step(&my_iplan);
      vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
	    "approximate solution, vector f_hat_iter");
      
      SWAP_complex(my_iplan.f_hat_iter,my_plan.f_hat);
      nfft_trafo(&my_plan);
      vpr_complex(my_plan.f,my_plan.M_total,"fitting the data, vector f");
      SWAP_complex(my_iplan.f_hat_iter,my_plan.f_hat);
    }
  
  infft_finalize(&my_iplan);  
  nfft_finalize(&my_plan);  
}

void measure_time_nfft(int d, int N, unsigned test_ndft)
{
  int j,k,r;
  int M,NN[d],nn[d];

  nfft_plan p;
  fftw_plan p_fft;

  double t,t_fft,t_ndft,t_nfft;

  printf("$%d$&\t$%d$&\t",d,(int)(log(N)/log(2)+0.5));

  for(j=0,M=1;j<d;j++)
    {
      M=N*M;
      NN[j]=N;
      nn[j]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 4, 
		 PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
  /*  nfft_init(&p, d, NN, M);*/
  p_fft=fftw_plan_dft(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  /** init pseudo random nodes */
  for(j=0;j<p.d*p.M_total;j++)
      p.x[j]=drand48()-0.5;

  nfft_precompute_one_psi(&p);

  /** init pseudo random Fourier coefficients */
  for(k=0;k<p.N_total;k++)
    p.f_hat[k] = drand48() + I* drand48();

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=second();
      fftw_execute(p_fft);
      t=second()-t;
      t_fft+=t;
    }
  t_fft/=r;
  printf("$%.2e$&\t",t_fft);

  /** NDFT */
  if(test_ndft)
    {
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=second();
          ndft_trafo(&p);
          t=second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("$%.2e$&\t",t_ndft);
    }
  else
    printf("*\t\t");

  
  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=second();
      nfft_trafo(&p);
      t=second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("$%.2e$\\\\\n",t_nfft);

  fflush(stdout);

  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 

int main()
{
  int l,m,d,logN;

  for(d=1;d<=3;d++)
    for(logN=3;logN<=22/d;logN++)
      if(logN*d<=12)
	measure_time_nfft(d,(1U<< logN),1);
      else
	measure_time_nfft(d,(1U<< logN),0);
  
  exit(-1);

  system("clear");
  printf("1) computing an one dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_1d();

  getc(stdin);

  system("clear"); 
  printf("2) computing a two dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_2d();
/*
  getc(stdin);

  system("clear"); 
  printf("3) computing an one dimensional infft\n\n");
  simple_test_infft_1d();

  getc(stdin);
    
  system("clear"); 
  printf("3) computing times for one dimensional nfft\n\n");
  measure_time_nfft_1d();
*/
  return 1;
}
