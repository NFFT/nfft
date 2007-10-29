#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

void measure_time_nfft(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  double t, t_fft, t_ndft, t_nfft;

  nfft_plan p;
  fftw_plan p_fft;

  printf("$%d$&\t",(int)(log(N)/log(2)*d+0.5));
     
  for(r=0,M=1;r<d;r++)
    {
      M=N*M;
      NN[r]=N;
      nn[r]=2*N;
    }

  nfft_init_guru(&p, d, NN, M, nn, 4, 
		 PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

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

  printf("$%.1e$&\t",t_fft);

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
      printf("$%.1e$&\t",t_ndft);
    }
  else
    printf("*&\t\t");

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
  
  if(d==1)
    printf("$%.1e$&\t",t_nfft);
  else
    printf("$%.1e$ \\\\\n",t_nfft);

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

  nfft_init_guru(&p, d, NN, M, nn, 5, 
		 PRE_PHI_HUT|
		 PRE_FG_PSI|
		 MALLOC_F_HAT| MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

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

  printf("\n",t_nfft);

  fftw_free(swapndft);
  fftw_destroy_plan(p_fft);
  nfft_finalize(&p);
} 

int main()
{
  int l,d,logIN;

  for(l=3;l<=12;l++)
    {
      logIN=l;
      if(logIN<=15)
	measure_time_nfft_XXX2(1,(1U<< (logIN)),1);
      else
	measure_time_nfft_XXX2(1,(1U<< (logIN)),0);
    }

  exit(-1);

  printf("\\multicolumn{4}{c|}{$d=1$}&\t\\multicolumn{4}{c|}{$d=2$}\\\\\n");
  for(l=3;l<=22;l++)
    {
      d=1;
      logIN=l;
      if(logIN<=15)
	measure_time_nfft(d,(1U<< (logIN/d)),1);
      else
	measure_time_nfft(d,(1U<< (logIN/d)),0);

      if(l<12)
	{
	  d=2;
	  logIN=d*l;
	  if(logIN<=15)
	    measure_time_nfft(d,(1U<< (logIN/d)),1);
	  else
	    measure_time_nfft(d,(1U<< (logIN/d)),0);
	}

      if(l==12)
	printf("\\multicolumn{4}{c|}{$d=3$}\\\\\n");
       
      if((l>12)&&(l<=17))
	{
	  d=3;
	  logIN=d*(l-10);
	  if(logIN<=15)
	    measure_time_nfft(d,(1U<< (logIN/d)),1);
	  else
	    measure_time_nfft(d,(1U<< (logIN/d)),0);
	}

      if(l>17)
	printf(" & & & \\\\\n");
	
	fflush(stdout);
      }

  return 1;
}
