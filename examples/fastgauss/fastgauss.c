#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <time.h>
#include <complex.h>

#include "nfft3.h"
#include "util.h"

/** @defgroup fast gauss transform
 * contains the gauss transforms with and without precomputation
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

void fgt_init(fgt_plan *ths, int N, int M, complex sigma, int n, double p, unsigned flags)
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

  if(ths->flags & DGT_PRE_CEXP)
   {
     ths->pre_cexp=(complex*)fftw_malloc(my_plan.M*my_plan.N*sizeof(complex));

     for(j=0,l=0; j<ths->M; j++)
       for(k=0; k<ths->N; k++,l++)
         ths->pre_cexp[l]=cexp(-ths->sigma*(ths->y[j]-ths->x[k])*
                                (ths->y[j]-ths->x[k]));
   }

/* TODO: ESTIMATE n, p for simple init !!!!!!!!!!!!!!!!!!!!!!*/
  ths->n = n;
  ths->p = p;

  ths->b = (complex*)fftw_malloc(ths->n*sizeof(complex));

  ths->nplan1 = (nfft_plan*) fftw_malloc(sizeof(nfft_plan));
  ths->nplan2 = (nfft_plan*) fftw_malloc(sizeof(nfft_plan));

  n_fftw=next_power_of_2(2*ths->n);

/* TODO: PRE_PHI_HUT only once?, FFTW_INIT only once?, PRE_LIN_PSI instead of PRE_PSI and only once? */
  nfft_init_guru(ths->nplan1, 1, &(ths->n), ths->N, &n_fftw, 7, PRE_PHI_HUT|
                 PRE_PSI| MALLOC_X| MALLOC_F_HAT| FFTW_INIT, FFTW_MEASURE);
  nfft_init_guru(ths->nplan2, 1, &(ths->n), ths->M, &n_fftw, 7, PRE_PHI_HUT|
                 PRE_PSI| MALLOC_X| FFTW_INIT, FFTW_MEASURE);

  for(j=0; j<ths->nplan1->M_total; j++)
    ths->nplan1->x[j] = ths->x[j]/ths->p;
  for(j=0; j<ths->nplan2->M_total; j++)
    ths->nplan2->x[j] = ths->y[j]/ths->p;
  
  if(ths->nplan1->nfft_flags & PRE_PSI)
    nfft_precompute_psi(ths->nplan1);
  if(ths->nplan2->nfft_flags & PRE_PSI)
    nfft_precompute_psi(ths->nplan2);

  ths->nplan1->f = ths->alpha;
  ths->nplan2->f_hat = ths->nplan1->f_hat;
  ths->nplan2->f = ths->f;

  if(ths->flags & FGT_APPROX_B)
    {
      ths->fplan = fftw_plan_dft_1d(ths->n, ths->b, ths->b, FFTW_FORWARD,
                                    FFTW_MEASURE);

      for(j=0; j<ths->n; j++)
	ths->b[j] = cexp(-ths->p*ths->p*ths->sigma*(j-ths->n/2)*(j-ths->n/2)/
                          ((double)ths->n*ths->n)) / ths->n;
      
      fftshift_complex(ths->b, ths->n);  
      fftw_execute(ths->fplan);
      fftshift_complex(ths->b, ths->n);
      
      fftw_destroy_plan(ths->fplan);
    }
  else
    {
      for(j=0; j<ths->n; j++)
	ths->b[j] = 1.0/ths->p * csqrt(PI/ths->sigma)*
	  cexp(-PI*PI*(j-ths->n/2)*(j-ths->n/2)/
	       (ths->p*ths->p*ths->sigma));
    }

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

void TODO()
{
  for(k=0; k<ths->N; k++)
    ths->x[k] = (double)rand()/(2.0*RAND_MAX)-1.0/4.0;

  for(j=0; j<ths->M; j++)
    ths->y[j] = (double)rand()/(2.0*RAND_MAX)-1.0/4.0;

  for(k=0; k<ths->N; k++)
    {
      ths->alpha[k] = (double)rand()/(RAND_MAX)-1.0/2.0
	             + I*(double)rand()/(RAND_MAX)-I/2.0;
    }
}

/** fgt_andersson_test similar to the test in
 *  F. Andersson and G. Beylkin.
 *  The fast Gauss transform with complex parameters.
 *  J. Comput. Physics, to appear.
 */
void fgt_andersson_test()
{
  fgt_plan my_plan;
  double t;
  complex *direct;
  int N;
  int r;

  complex sigma=4*(138+I*100);
  int n=128;

  double tau=0.01;

  printf("n=%d, sigma=%1.3e+i%1.3e\n",n,creal(sigma),cimag(sigma));

  for(N=((int)(1U<<6)); N<((int)(1U<<22)); N=N<<1)
    {
      printf("$%d$\t & ",N);
      fgt_init(&my_plan, N, N, sigma, n, 1, DGT_PRE_CEXP);

      if(N<((int)(1U<<19)))
	{
          direct = (complex*)fftw_malloc(my_plan.M*sizeof(complex));
          SWAP_complex(direct,my_plan.f);
          my_plan.flags^=DGT_PRE_CEXP;
          t_dgt=0;
          r=0; 
          while(t_dgt<0.1)
            {
              r++;
              t=second();
              dgt_trafo(&my_plan);
              t_dgt+=second()-t;
            }
          t_dgt/=r;

	  printf("$%1.1e$\t & ",t_dgt);
          my_plan.flags^=DGT_PRE_CEXP;
          SWAP_complex(direct,my_plan.f);  
	}
      else
	printf("\t\t & ");	

      if(N<((int)(1U<<14)))
	{
          t_dgt_pre=0;
          r=0; 
          while(t_dgt_pre<0.1)
            {
              r++;
              t=second();
              dgt_trafo(&my_plan);
              t_dgt_pre+=second()-t;
            }
          t_dgt_pre/=r;

	  printf("$%1.1e$\t & ",t_dgt_pre);
	}
      else
	printf("\t\t & ");	

       my_plan.flags^=FGT_NDFT;
       t_fgt_ndft=0;
       r=0; 
       while(t_fgt_ndft<0.1)
         {
           r++;
           t=second();
           dgt_trafo(&my_plan);
           t_fgt_ndft+=second()-t;
         }
       t_fgt_ndft/=r;

       printf("$%1.1e$\t & ",t_fgt_ndft);
       my_plan.flags^=FGT_NDFT;
      
       t_fgt_nfft=0;
       r=0; 
       while(t_fgt_nfft<0.1)
         {
           r++;
           t=second();
           dgt_trafo(&my_plan);
           t_fgt_nfft+=second()-t;
         }
       t_fgt_nfft/=r;

       printf("$%1.1e$\t & ",t_fgt_nfft);


      printf("$%1.1e$\t \\\\ \n", error_l_infty_1_complex(direct, my_plan.f,
             my_plan.M, my_plan.alpha, my_plan.N));
      fflush(stdout);

      fgt_finalize(&my_plan);
      fftw_cleanup();
    }
}


int main()
{
  /**  pipe to output_andersson.tex */
  fgt_andersson_test();

  /** pipe to output_error.m */
  //fgt_error();

  /** pipe to output_error_p.m */
  //fgt_error_p();

  return 1;
}
