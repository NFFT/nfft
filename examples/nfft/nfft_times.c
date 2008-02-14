#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

void measure_time_nfft(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d],j;
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  double auxC=pow(2,29);

  printf("\\verb+%d+&\t",(int)(log(N)/log(2)*d+0.5));
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 4, 
		 PRE_PHI_HUT| 
		 PRE_PSI| MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  p_fft=fftw_plan_dft(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  nfft_precompute_one_psi(&p);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=nfft_second();
      fftw_execute(p_fft);
      t=nfft_second()-t;
      t_fft+=t;
    }
  t_fft/=r;

  printf("\\verb+%.1e+ & \\verb+(%.1f)+ &\t",t_fft,t_fft/(p.N_total*(log(N)/log(2)*d))*auxC);

  /** NDFT */
  if(test_ndft)
    {
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_trafo(&p);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("\\verb+%.1e+ & \\verb+(%d)+&\t",t_ndft,(int)round(t_ndft/(p.N_total*p.N_total)*auxC));
    }
  else
    printf("\\verb+*+\t&\t&\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      switch(d)
	{
	  case 1: { nfft_trafo_1d(&p); break; }
	  case 2: { nfft_trafo_2d(&p); break; }
	  case 3: { nfft_trafo_3d(&p); break; }
          default: nfft_trafo(&p);
	}
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  
  printf("\\verb+%.1e+ & \\verb+(%d)+\\\\\n",t_nfft,(int)round(t_nfft/(p.N_total*(log(N)/log(2)*d))*auxC));

  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 

static int comp(const void *x,const void *y)
{
  return ((* (double*) x)<(* (double*) y)?-1:1);
}

void measure_time_nfft_XXX2(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  printf("%d\t",(int)(log(N)/log(2)*d+0.5)); fflush(stdout);
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 2, 
		 PRE_PHI_HUT|
		 PRE_PSI|
		 MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  p_fft=fftw_plan_dft(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  double _Complex *swapndft=(double _Complex*)fftw_malloc(p.M_total*sizeof(double _Complex));

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  qsort(p.x,p.M_total,sizeof(double),comp);
  //nfft_vpr_double(p.x,p.M_total,"nodes x"); 

  nfft_precompute_one_psi(&p);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=nfft_second();
      fftw_execute(p_fft);
      t=nfft_second()-t;
      t_fft+=t;
    }
  t_fft/=r;

  printf("%.1e\t",t_fft);

  /** NDFT */
  if(test_ndft)
    {
      NFFT_SWAP_complex(p.f,swapndft);
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_trafo(&p);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("%.1e\t",t_ndft);
      NFFT_SWAP_complex(p.f,swapndft);
    }
  else
    printf("\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_trafo(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f, p.M_total));

  /** NFFT_1d */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_trafo_1d(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f, p.M_total));

  printf("\n");

  fftw_free(swapndft);
  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 

void measure_time_nfft_XXX3(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  printf("%d\t",(int)(log(N)/log(2)*d+0.5)); fflush(stdout);
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 2, 
		 PRE_PHI_HUT|
		 PRE_PSI|
		 MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  p_fft=fftw_plan_dft(d, NN, p.f, p.f_hat, FFTW_BACKWARD, FFTW_MEASURE);

  double _Complex *swapndft=(double _Complex*)fftw_malloc(p.N_total*sizeof(double _Complex));

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  qsort(p.x,p.M_total,sizeof(double),comp);
  //nfft_vpr_double(p.x,p.M_total,"nodes x"); 

  nfft_precompute_one_psi(&p);

  /** init pseudo random samples */
  nfft_vrand_unit_complex(p.f, p.N_total);

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=nfft_second();
      fftw_execute(p_fft);
      t=nfft_second()-t;
      t_fft+=t;
    }
  t_fft/=r;

  printf("%.1e\t",t_fft);

  /** NDFT */
  if(test_ndft)
    {
      NFFT_SWAP_complex(p.f_hat,swapndft);
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_adjoint(&p);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("%.1e\t",t_ndft);
      NFFT_SWAP_complex(p.f_hat,swapndft);
    }
  else
    printf("\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_adjoint(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f_hat, p.N_total));

  /** NFFT_1d */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_adjoint_1d(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f_hat, p.N_total));

  printf("\n");

  fftw_free(swapndft);
  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 




void measure_time_nfft_XXX4(int d, int N, unsigned test_ndft)
{
  int j,r, M, NN[d], nn[d];
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  printf("%d\t",(int)(log(N)/log(2)*d+0.5)); fflush(stdout);
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 4, 
		 PRE_PHI_HUT|
		 PRE_PSI|
		 MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  p_fft=fftw_plan_dft(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  double _Complex *swapndft=(double _Complex*)fftw_malloc(p.M_total*sizeof(double _Complex));

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  //for(j=0;j<2*M;j+=2)
  //   p.x[j]=0.5*p.x[j]+0.25*(p.x[j]>=0?1:-1);

  //qsort(p.x,p.M_total,sizeof(double),comp);
  //nfft_vpr_double(p.x,p.M_total,"nodes x"); 

  nfft_precompute_one_psi(&p);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=nfft_second();
      fftw_execute(p_fft);
      t=nfft_second()-t;
      t_fft+=t;
    }
  t_fft/=r;

  printf("%.1e\t",t_fft);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** NDFT */
  if(test_ndft)
    {
      NFFT_SWAP_complex(p.f,swapndft);
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_trafo(&p);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("%.1e\t",t_ndft);

      //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

      NFFT_SWAP_complex(p.f,swapndft);
    }
  else
    printf("\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_trafo(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  /** NFFT_2d */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_trafo_2d(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  printf("\n");

  fftw_free(swapndft);
  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 


void measure_time_nfft_XXX5(int d, int N, unsigned test_ndft)
{
  int j,r, M, NN[d], nn[d];
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  printf("%d\t",(int)(log(N)/log(2)*d+0.5)); fflush(stdout);
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 4, 
		 PRE_PHI_HUT|
		 PRE_PSI|
		 MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  p_fft=fftw_plan_dft(d, NN, p.f, p.f_hat, FFTW_FORWARD, FFTW_MEASURE);

  double _Complex *swapndft=(double _Complex*)fftw_malloc(p.N_total*sizeof(double _Complex));

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  //sort_nodes(p.x,p.d,p.M_total,

  nfft_precompute_one_psi(&p);

  /** init pseudo random samples */
  nfft_vrand_unit_complex(p.f, p.M_total);

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=nfft_second();
      fftw_execute(p_fft);
      t=nfft_second()-t;
      t_fft+=t;
    }
  t_fft/=r;

  printf("%.1e\t",t_fft);

  /** init pseudo random samples */
  nfft_vrand_unit_complex(p.f, p.M_total);

  /** NDFT */
  if(test_ndft)
    {
      NFFT_SWAP_complex(p.f_hat,swapndft);
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_adjoint(&p);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("%.1e\t",t_ndft);

      //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

      NFFT_SWAP_complex(p.f_hat,swapndft);
    }
  else
    printf("\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_adjoint(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  /** NFFT_2d */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_adjoint_2d(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  printf("\n");

  fftw_free(swapndft);
  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 


void measure_time_nfft_XXX6(int d, int N, unsigned test_ndft)
{
  int j,r, M, NN[d], nn[d];
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  printf("%d\t",(int)(log(N)/log(2)*d+0.5)); fflush(stdout);
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 2, 
		 PRE_PHI_HUT|
		 FG_PSI|
		 MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  p_fft=fftw_plan_dft(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  double _Complex *swapndft=(double _Complex*)fftw_malloc(p.M_total*sizeof(double _Complex));

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  //qsort(p.x,p.M_total,sizeof(double),comp);
  //nfft_vpr_double(p.x,p.M_total,"nodes x"); 

  nfft_precompute_one_psi(&p);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=nfft_second();
      fftw_execute(p_fft);
      t=nfft_second()-t;
      t_fft+=t;
    }
  t_fft/=r;

  printf("%.1e\t",t_fft);

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(p.f_hat, p.N_total);

  /** NDFT */
  if(test_ndft)
    {
      NFFT_SWAP_complex(p.f,swapndft);
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_trafo(&p);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("%.1e\t",t_ndft);

      //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

      NFFT_SWAP_complex(p.f,swapndft);
    }
  else
    printf("\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_trafo(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  /** NFFT_3d */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_trafo_3d(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  printf("\n");

  fftw_free(swapndft);
  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 


void measure_time_nfft_XXX7(int d, int N, unsigned test_ndft)
{
  int j,r, M, NN[d], nn[d];
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  printf("%d\t",(int)(log(N)/log(2)*d+0.5)); fflush(stdout);
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 2, 
		 PRE_PHI_HUT|
		 FG_PSI|
		 MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  p_fft=fftw_plan_dft(d, NN, p.f, p.f_hat, FFTW_FORWARD, FFTW_MEASURE);

  double _Complex *swapndft=(double _Complex*)fftw_malloc(p.N_total*sizeof(double _Complex));

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x, p.d*p.M_total);

  //sort_nodes(p.x,p.d,p.M_total,

  nfft_precompute_one_psi(&p);

  /** init pseudo random samples */
  nfft_vrand_unit_complex(p.f, p.M_total);

  /** FFT */
  t_fft=0;
  r=0;
  while(t_fft<0.1)
    {
      r++;
      t=nfft_second();
      fftw_execute(p_fft);
      t=nfft_second()-t;
      t_fft+=t;
    }
  t_fft/=r;

  printf("%.1e\t",t_fft);

  /** init pseudo random samples */
  nfft_vrand_unit_complex(p.f, p.M_total);

  /** NDFT */
  if(test_ndft)
    {
      NFFT_SWAP_complex(p.f_hat,swapndft);
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_adjoint(&p);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;
      printf("%.1e\t",t_ndft);

      //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

      NFFT_SWAP_complex(p.f_hat,swapndft);
    }
  else
    printf("\t");

  /** NFFT */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_adjoint(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  /** NFFT_3d */
  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_adjoint_3d(&p);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;
  printf("%.1e\t",t_nfft);
  if(test_ndft)
    printf("(%.1e)\t",nfft_error_l_2_complex(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  printf("\n");

  fftw_free(swapndft);
  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 

int main2()
{
  int l,d,logIN;

  for(l=3;l<=6;l++)
    {
      d=3;
      logIN=d*l;
      if(logIN<=15)
	{
	  measure_time_nfft_XXX6(d,(1U<< (logIN/d)),1);
	  measure_time_nfft_XXX7(d,(1U<< (logIN/d)),1);
	}
      else
	{
	  measure_time_nfft_XXX6(d,(1U<< (logIN/d)),0);
	  measure_time_nfft_XXX7(d,(1U<< (logIN/d)),0);
	}
    }
  
  exit(-1);


  for(l=7;l<=12;l++)
    {
      d=2;
      logIN=d*l;
      if(logIN<=15)
	{
	  measure_time_nfft_XXX4(d,(1U<< (logIN/d)),1);
	  measure_time_nfft_XXX5(d,(1U<< (logIN/d)),1);
	}
      else
	{
	  measure_time_nfft_XXX4(d,(1U<< (logIN/d)),0);
	  measure_time_nfft_XXX5(d,(1U<< (logIN/d)),0);
	}
    }
  
  exit(-1);

  for(l=3;l<=12;l++)
    {
      logIN=l;
      if(logIN<=15)
	{
	  measure_time_nfft_XXX2(1,(1U<< (logIN)),1);
	  measure_time_nfft_XXX3(1,(1U<< (logIN)),1);
	}
      else
	{
	  measure_time_nfft_XXX2(1,(1U<< (logIN)),0);
	  measure_time_nfft_XXX3(1,(1U<< (logIN)),0);
	} 
    }

  exit(-1);
}

int main()
{
  int l,d,logIN;

  printf("\\hline $l_N$ & FFT & $/(2^{l_N} l_N)$ & NDFT & $/4^{l_N}$ & NFFT $/(2^{l_N} l_N)$ \\\\\n");
  printf("\\hline \\hline \\multicolumn{7}{|c|}{$d=1$} \\\\ \\hline\n");
  for(l=3;l<=22;l++)
    {
      d=1;
      logIN=l;
      if(logIN<=15)
	measure_time_nfft(d,(1U<< (logIN/d)),1);
      else
	measure_time_nfft(d,(1U<< (logIN/d)),0);
      
      fflush(stdout);
    }
  
  printf("\\hline $l_N$ & FFT & $/(2^{l_N} l_N)$ & NDFT & $/4^{l_N}$ & NFFT $/(2^{l_N} l_N)$ \\\\ \n");
  printf("\\hline \\hline \\multicolumn{7}{|c|}{$d=2$} \\\\ \\hline\n");
  for(l=3;l<=11;l++)
    {
      d=2;
      logIN=d*l;
      if(logIN<=15)
	measure_time_nfft(d,(1U<< (logIN/d)),1);
      else
	measure_time_nfft(d,(1U<< (logIN/d)),0);

      fflush(stdout);
    }

  printf("\\hline \\hline \\multicolumn{7}{|c|}{$d=3$} \\\\ \\hline\n");
  for(l=3;l<=7;l++)
    {
      d=3;
      logIN=d*l;
      if(logIN<=15)
	measure_time_nfft(d,(1U<< (logIN/d)),1);
      else
	measure_time_nfft(d,(1U<< (logIN/d)),0);

      fflush(stdout);
    }

  return 1;
}
