#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <time.h>
#include <complex.h>

#include "nfft3.h"
#include "util.h"

/** @defgroup fast gauss transform
 * 
 * @{ 
 */
#define DGT_PRE_CEXP     (1U<< 0)
#define FGT_NDFT         (1U<< 1)
#define FGT_APPROX_B     (1U<< 2)

typedef struct fgt_plan_ 
{
  int N;                                /**< number of source knots          */
  int M;                                /**< number of target knots          */

  complex *alpha;                       /**< source coefficients             */
  complex *f;                           /**< target evaluations              */

  unsigned flags;                       /**< flags precomp. and approx.type  */

  complex sigma;                        /**< parameter of the Gaussian       */

  double *x;                            /**< source knots in [-1/4,1/4]      */
  double *y;                            /**< target knots in [-1/4,1/4]      */

  complex *pre_cexp;                    /**< precomputed values for dgt      */

  int n;                                /**< expansion degree                */
  double p;                             /**< period, at least 1              */

  complex *b;                           /**< expansion coefficients          */

  nfft_plan *nplan1;                    /**< source nfft plan                */
  nfft_plan *nplan2;                    /**< target nfft plan                */

} fgt_plan;

void dgt_trafo(fgt_plan *ths)
{
  int j,k,l;
  
  for(j=0; j<ths->M; j++)
    ths->f[j] = 0;
  
  if(ths->flags & DGT_PRE_CEXP)
    for(j=0,l=0; j<ths->M; j++)
      for(k=0; k<ths->N; k++,l++)
        ths->f[j] += ths->alpha[k]*ths->pre_cexp[l];
  else
    for(j=0; j<ths->M; j++)
      for(k=0; k<ths->N; k++)
        ths->f[j] += ths->alpha[k]*cexp(-ths->sigma*(ths->y[j]-ths->x[k])*
					(ths->y[j]-ths->x[k]));
}

void fgt_trafo(fgt_plan *ths)
{
  int l;

  if(ths->flags & FGT_NDFT)
    {
      ndft_adjoint(ths->nplan1);

      for(l=0; l<ths->n; l++)
        ths->nplan1->f_hat[l] *= ths->b[l];
  
      ndft_trafo(ths->nplan2);
    }
  else
    {
      nfft_adjoint(ths->nplan1);

      for(l=0; l<ths->n; l++)
        ths->nplan1->f_hat[l] *= ths->b[l];
  
      nfft_trafo(ths->nplan2);
    }
}

void fgt_init_guru(fgt_plan *ths, int N, int M, complex sigma, int n, double p, int m, unsigned flags)
{
  int j,k,l,n_fftw;
  fftw_plan fplan;

  ths->M = M;
  ths->N = N;
  ths->sigma = sigma;
  ths->flags = flags;

  ths->x = (double*)fftw_malloc(ths->N*sizeof(double));
  ths->y = (double*)fftw_malloc(ths->M*sizeof(double));
  ths->alpha = (complex*)fftw_malloc(ths->N*sizeof(complex));
  ths->f = (complex*)fftw_malloc(ths->M*sizeof(complex));



/* TODO: ESTIMATE n, p for simple init !!!!!!!!!!!!!!!!!!!!!!*/
  ths->n = n;
  ths->p = p;

  ths->b = (complex*)fftw_malloc(ths->n*sizeof(complex));

  ths->nplan1 = (nfft_plan*) fftw_malloc(sizeof(nfft_plan));
  ths->nplan2 = (nfft_plan*) fftw_malloc(sizeof(nfft_plan));

  n_fftw=next_power_of_2(2*ths->n);

/* TODO: PRE_PHI_HUT only once?, FFTW_INIT only once?, PRE_LIN_PSI instead of PRE_PSI and only once? */
  nfft_init_guru(ths->nplan1, 1, &(ths->n), ths->N, &n_fftw, m, PRE_PHI_HUT|
                 PRE_PSI| MALLOC_X| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
  nfft_init_guru(ths->nplan2, 1, &(ths->n), ths->M, &n_fftw, m, PRE_PHI_HUT|
                 PRE_PSI| MALLOC_X| FFTW_INIT, FFTW_MEASURE);

  ths->nplan1->f = ths->alpha;
  ths->nplan2->f_hat = ths->nplan1->f_hat;
  ths->nplan2->f = ths->f;

  if(ths->flags & FGT_APPROX_B)
    {
      fplan = fftw_plan_dft_1d(ths->n, ths->b, ths->b, FFTW_FORWARD,
                               FFTW_MEASURE);

      for(j=0; j<ths->n; j++)
	ths->b[j] = cexp(-ths->p*ths->p*ths->sigma*(j-ths->n/2)*(j-ths->n/2)/
                          ((double)ths->n*ths->n)) / ths->n;
      
      fftshift_complex(ths->b, 1, &ths->n);  
      fftw_execute(fplan);
      fftshift_complex(ths->b, 1, &ths->n);
      
      fftw_destroy_plan(fplan);
    }
  else
    {
      for(j=0; j<ths->n; j++)
	ths->b[j] = 1.0/ths->p * csqrt(PI/ths->sigma)*
	  cexp(-PI*PI*(j-ths->n/2)*(j-ths->n/2)/
	       (ths->p*ths->p*ths->sigma));
    }
}

void fgt_init(fgt_plan *ths, int N, int M, complex sigma, double eps)
{
  double p;
  int n;

  p=0.5+sqrt(-log(eps)/creal(sigma));
  if(p<1)
    p=1;

  n=2*((int)ceil(p*cabs(sigma)/PI * sqrt(-log(eps)/creal(sigma))));

  if(N*M<=((int)(1U<<20)))
    fgt_init_guru(ths, N, M, sigma, n, p, 7, DGT_PRE_CEXP);
  else
   fgt_init_guru(ths, N, M, sigma, n, p, 7, 0);
}

void fgt_init_knot_dependent(fgt_plan *ths)
{
  int j,k,l;

  if(ths->flags & DGT_PRE_CEXP)
   {
     ths->pre_cexp=(complex*)fftw_malloc(ths->M*ths->N*sizeof(complex));

     for(j=0,l=0; j<ths->M; j++)
       for(k=0; k<ths->N; k++,l++)
         ths->pre_cexp[l]=cexp(-ths->sigma*(ths->y[j]-ths->x[k])*
                                (ths->y[j]-ths->x[k]));
   }

  for(j=0; j<ths->nplan1->M_total; j++)
    ths->nplan1->x[j] = ths->x[j]/ths->p;
  for(j=0; j<ths->nplan2->M_total; j++)
    ths->nplan2->x[j] = ths->y[j]/ths->p;
  
  if(ths->nplan1->nfft_flags & PRE_PSI)
    nfft_precompute_psi(ths->nplan1);
  if(ths->nplan2->nfft_flags & PRE_PSI)
    nfft_precompute_psi(ths->nplan2);
}

void fgt_finalize(fgt_plan *ths)
{
  nfft_finalize(ths->nplan2);
  nfft_finalize(ths->nplan1);

  fftw_free(ths->nplan2);
  fftw_free(ths->nplan1);

  fftw_free(ths->b);

  fftw_free(ths->f);
  fftw_free(ths->y);

  fftw_free(ths->alpha);
  fftw_free(ths->x);
}
/** @} 
 */

/** @defgroup test gauss transform
 * contains the tests
 * @{ 
 */

void fgt_test_init_rand(fgt_plan *ths)
{
  int j,k;

  for(k=0; k<ths->N; k++)
    ths->x[k] = (double)rand()/(2.0*RAND_MAX)-1.0/4.0;

  for(j=0; j<ths->M; j++)
    ths->y[j] = (double)rand()/(2.0*RAND_MAX)-1.0/4.0;

  for(k=0; k<ths->N; k++)
    ths->alpha[k] =   (double)rand()/(RAND_MAX)-1.0/2.0
	          + I*(double)rand()/(RAND_MAX)-I/2.0;
}

double fgt_test_measure_time(fgt_plan *ths, unsigned dgt)
{
  int r;
  double t_out,t; 
  double tau=0.01;

  t_out=0;
  r=0; 
  while(t_out<tau)
    {
      r++;
      if(dgt)
        {
          t=second();
          dgt_trafo(ths);
        }
      else
        {
          t=second();
          fgt_trafo(ths);
        }
      t_out+=second()-t;
    }
  t_out/=r;

  return t_out;
}

void fgt_test_simple(int N, int M, complex sigma, double eps)
{
  fgt_plan my_plan;
  complex *swap_dgt;
     
  fgt_init(&my_plan, N, M, sigma, eps);
  swap_dgt = (complex*)fftw_malloc(my_plan.M*sizeof(complex));

  fgt_test_init_rand(&my_plan);
  fgt_init_knot_dependent(&my_plan);   

  SWAP_complex(swap_dgt,my_plan.f);
  dgt_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M,"discrete gauss transform");
  SWAP_complex(swap_dgt,my_plan.f);

  fgt_trafo(&my_plan);
  vpr_complex(my_plan.f,my_plan.M,"fast gauss transform");

  printf("\n relative error: %1.3e\n", error_l_infty_1_complex(swap_dgt,
         my_plan.f, my_plan.M, my_plan.alpha, my_plan.N));

  fftw_free(swap_dgt);
  fgt_finalize(&my_plan);
}

/** fgt_andersson_test similar to the test in
 *  F. Andersson and G. Beylkin.
 *  The fast Gauss transform with complex parameters.
 *  J. Comput. Physics 203 (2005) 274-286
 */
void fgt_test_andersson()
{
  fgt_plan my_plan;
  complex *swap_dgt;
  int N;
  int r,k,j;

  complex sigma=4*(138+I*100);
  int n=128;
  int N_dgt_pre_exp=(int)(1U<<11);
  int N_dgt=(int)(1U<<19);

  printf("n=%d, sigma=%1.3e+i%1.3e\n",n,creal(sigma),cimag(sigma));

  for(N=((int)(1U<<6)); N<((int)(1U<<22)); N=N<<1)
    {
      printf("$%d$\t & ",N);

      if(N<N_dgt_pre_exp)
        fgt_init_guru(&my_plan, N, N, sigma, n, 1, 7, DGT_PRE_CEXP);
      else
        fgt_init_guru(&my_plan, N, N, sigma, n, 1, 7, 0);

      swap_dgt = (complex*)fftw_malloc(my_plan.M*sizeof(complex));

      fgt_test_init_rand(&my_plan);
      
      fgt_init_knot_dependent(&my_plan);      

      if(N<N_dgt)
	{
          SWAP_complex(swap_dgt,my_plan.f);
          if(N<N_dgt_pre_exp)
            my_plan.flags^=DGT_PRE_CEXP;
 
	  printf("$%1.1e$\t & ",fgt_test_measure_time(&my_plan, 1));
          if(N<N_dgt_pre_exp)
            my_plan.flags^=DGT_PRE_CEXP;
          SWAP_complex(swap_dgt,my_plan.f);  
	}
      else
	printf("\t\t & ");	

      if(N<N_dgt_pre_exp)
	printf("$%1.1e$\t & ",fgt_test_measure_time(&my_plan, 1));
      else
	printf("\t\t & ");	

      my_plan.flags^=FGT_NDFT;
      printf("$%1.1e$\t & ",fgt_test_measure_time(&my_plan, 0));
      my_plan.flags^=FGT_NDFT;

      printf("$%1.1e$\t & ",fgt_test_measure_time(&my_plan, 0));

      printf("$%1.1e$\t \\\\ \n", error_l_infty_1_complex(swap_dgt, my_plan.f,
             my_plan.M, my_plan.alpha, my_plan.N));
      fflush(stdout);

      fftw_free(swap_dgt);
      fgt_finalize(&my_plan);
      fftw_cleanup();
    }
}

void fgt_test_error()
{
  fgt_plan my_plan;
  complex *swap_dgt;
  int n,mi;

  complex sigma=4*(138+I*100);
  int N=1000;
  int M=1000;
  int m[2]={7,3};
  
  printf("N=%d;\tM=%d;\nsigma=%1.3e+i*%1.3e;\n",N,M,creal(sigma),cimag(sigma));
  printf("error=[\n");

  swap_dgt = (complex*)fftw_malloc(M*sizeof(complex));

  for(n=8; n<=128; n+=4)
    {
      printf("%d\t",n);
      for(mi=0;mi<2;mi++)
        {
          fgt_init_guru(&my_plan, N, M, sigma, n, 1, m[mi], 0);
          fgt_test_init_rand(&my_plan);    
          fgt_init_knot_dependent(&my_plan);    

          SWAP_complex(swap_dgt,my_plan.f);
          dgt_trafo(&my_plan);
          SWAP_complex(swap_dgt,my_plan.f);

          fgt_trafo(&my_plan);

          printf("%1.3e\t", error_l_infty_1_complex(swap_dgt, my_plan.f,
                 my_plan.M, my_plan.alpha, my_plan.N));
          fflush(stdout);

          fgt_finalize(&my_plan);
          fftw_cleanup();
        }
      printf("\n");
    }
  printf("];\n");

  fftw_free(swap_dgt);
}

void fgt_test_error_p()
{
  fgt_plan my_plan;
  complex *swap_dgt;
  int n,pi;

  complex sigma=20+40I;
  int N=1000;
  int M=1000;
  double p[3]={1,1.5,2};
  
  printf("N=%d;\tM=%d;\nsigma=%1.3e+i*%1.3e;\n",N,M,creal(sigma),cimag(sigma));
  printf("error=[\n");

  swap_dgt = (complex*)fftw_malloc(M*sizeof(complex));

  for(n=8; n<=128; n+=4)
    {
      printf("%d\t",n);
      for(pi=0;pi<3;pi++)
        {
          fgt_init_guru(&my_plan, N, M, sigma, n, p[pi], 7, 0);
          fgt_test_init_rand(&my_plan);    
          fgt_init_knot_dependent(&my_plan);    

          SWAP_complex(swap_dgt,my_plan.f);
          dgt_trafo(&my_plan);
          SWAP_complex(swap_dgt,my_plan.f);

          fgt_trafo(&my_plan);

          printf("%1.3e\t", error_l_infty_1_complex(swap_dgt, my_plan.f,
                 my_plan.M, my_plan.alpha, my_plan.N));
          fflush(stdout);

          fgt_finalize(&my_plan);
          fftw_cleanup();
        }
      printf("\n");
    }
  printf("];\n");  
}
/** @} 
 */

int main()
{
  /**  pipe to output_andersson.tex */
  fgt_test_simple();

  /**  pipe to output_andersson.tex */
//fgt_test_andersson();

  /** pipe to output_error.m */
//fgt_test_error();

  /** pipe to output_error_p.m */
//fgt_test_error_p();

  return 1;
}
