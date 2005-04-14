#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <time.h>
#include <complex.h>

#include "nfft.h"
#include "utils.h"

#define SAMPLED_FT       (1U<< 0)
#define DFT_APPROX       (1U<< 1)

typedef struct gauss_plan_ 
{
  int N,M;

  double p;

  fftw_complex delta;

  double *x;
  double *y;

  fftw_complex *alpha;
  fftw_complex *f;

  /** precomputed values of the Gaussian for the direct transform */
  fftw_complex *pre_cexp;

  /** data for the fast transform */
  int n;

  fftw_complex *b;

  fftw_plan my_fftw_plan;
  
  nfft_plan *my_nfft_plan1;
  nfft_plan *my_nfft_plan2;

} gauss_plan;

/** @defgroup utilities
 * contains utility routines
 * @{ 
 */
void fft_shift(fftw_complex *a, int n)
{
  int j;
  fftw_complex temp;

  for(j=0; j<n/2; j++)
    {
      temp=a[j];
      a[j]=a[j+n/2];
      a[j+n/2]=temp;
    }
}

double l_infty_error_relative_to_alpha(fftw_complex *x_0, fftw_complex *x,
				       int M, fftw_complex *alpha, int N)
{
  int j;
  double error,l1_alpha;

  error=0;
  for(j=0; j<M; j++) 
    if(cabs(x_0[j]-x[j])>error)
      error=cabs(x_0[j]-x[j]);

  l1_alpha=0;
  for(j=0; j<N; j++)
    l1_alpha += cabs(alpha[j]);
  
  return error/l1_alpha;
}
/** @} 
 */ 

/** @defgroup gauss transform
 * contains the gauss transforms with and without precomputation
 * @{ 
 */
void precompute_cexp(gauss_plan *this)
{
  int j,k,l;
  
  for(j=0,l=0; j<this->M; j++)
    for(k=0; k<this->N; k++,l++)
      this->pre_cexp[l]=cexp(-this->delta*(this->y[j]-this->x[k])*
			     (this->y[j]-this->x[k]));
}

void gauss_trafo(gauss_plan *this)
{
  int j,k;
  
  for(j=0; j<this->M; j++)
    this->f[j] = 0;
  
  for(j=0; j<this->M; j++)
    for(k=0; k<this->N; k++)
      this->f[j] += this->alpha[k]*cexp(-this->delta*(this->y[j]-this->x[k])*
					(this->y[j]-this->x[k]));
}

void gauss_trafo_pre(gauss_plan *this)
{
  int j,k,l;
  
  for(j=0; j<this->M; j++)
    this->f[j] = 0;
  
  for(j=0,l=0; j<this->M; j++)
    for(k=0; k<this->N; k++,l++)
      this->f[j] += this->alpha[k]*this->pre_cexp[l];
}
/** @} 
 */

/** @defgroup fast gauss transform
 * contains the fast gauss transforms ndft and nfft based
 * @{ 
 */
void fast_gauss_trafo_ndft(gauss_plan *this)
{
  int l;

  ndft_adjoint(this->my_nfft_plan1);

  for(l=0; l<this->n; l++)
    this->my_nfft_plan1->f_hat[l] *= this->b[l];
  
  ndft_trafo(this->my_nfft_plan2);
}

void fast_gauss_trafo_nfft(gauss_plan *this)
{
  int l;

  nfft_adjoint(this->my_nfft_plan1);

  for(l=0; l<this->n; l++)
    this->my_nfft_plan1->f_hat[l] *= this->b[l];
  
  nfft_trafo(this->my_nfft_plan2);
}
/** @} 
 */

/** @defgroup init gauss transform
 * contains the initialization and finalization
 * @{ 
 */
void init_gauss_trafo(gauss_plan *this, int N, int M, fftw_complex delta,
		      int n, double p)
{
  int j,k,sigma_n;

  this->M = M;
  this->N = N;
  this->p = p;
  this->delta = delta;

  this->x = (double*)fftw_malloc(this->N*sizeof(double));
  this->alpha = (fftw_complex*)fftw_malloc(this->N*sizeof(fftw_complex));

  this->y = (double*)fftw_malloc(this->M*sizeof(double));
  this->f = (fftw_complex*)fftw_malloc(this->M*sizeof(fftw_complex));

  for(k=0; k<this->N; k++)
    this->x[k] = (double)rand()/(2.0*RAND_MAX)-1.0/4.0;

  for(j=0; j<this->M; j++)
    this->y[j] = (double)rand()/(2.0*RAND_MAX)-1.0/4.0;

  for(k=0; k<this->N; k++)
    {
      this->alpha[k] = (double)rand()/(RAND_MAX)-1.0/2.0
	             + I*(double)rand()/(RAND_MAX)-I/2.0;
    }

  /* ... and for the fast trafo */
  this->n = n;

  this->b = (fftw_complex*)fftw_malloc(this->n*sizeof(fftw_complex));

  this->my_nfft_plan1 = (nfft_plan*) fftw_malloc(sizeof(nfft_plan));
  this->my_nfft_plan2 = (nfft_plan*) fftw_malloc(sizeof(nfft_plan));

  sigma_n=next_power_of_2(2*this->n);
  if(0)
    {
      nfft_init_specific(this->my_nfft_plan1, 1, &(this->n), this->N, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT,
			 FFTW_MEASURE);
      nfft_init_specific(this->my_nfft_plan2, 1, &(this->n), this->M, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI, FFTW_MEASURE);
      this->my_nfft_plan1->x = this->x;
      this->my_nfft_plan2->x = this->y;
    }
  else
    {
      nfft_init_specific(this->my_nfft_plan1, 1, &(this->n), this->N, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT,
			 FFTW_MEASURE);
      nfft_init_specific(this->my_nfft_plan2, 1, &(this->n), this->M, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_X,
			 FFTW_MEASURE);
      for(j=0; j<this->my_nfft_plan1->M; j++)
	this->my_nfft_plan1->x[j] = this->x[j]/this->p;
      for(j=0; j<this->my_nfft_plan2->M; j++)
	this->my_nfft_plan2->x[j] = this->y[j]/this->p;
    }
  
  if(this->my_nfft_plan1->nfft_flags & PRE_PSI)
    nfft_precompute_psi(this->my_nfft_plan1);

  this->my_nfft_plan1->f = this->alpha; 

  if(this->my_nfft_plan2->nfft_flags & PRE_PSI)
    nfft_precompute_psi(this->my_nfft_plan2);

  this->my_nfft_plan2->f_hat = this->my_nfft_plan1->f_hat;

  this->my_nfft_plan2->f = this->f;
}

void finalize_gauss_trafo(gauss_plan *this)
{
  nfft_finalize(this->my_nfft_plan2);
  nfft_finalize(this->my_nfft_plan1);

  fftw_free(this->my_nfft_plan2);
  fftw_free(this->my_nfft_plan1);

  fftw_free(this->b);

  fftw_free(this->f);
  fftw_free(this->y);

  fftw_free(this->alpha);
  fftw_free(this->x);
}

void init_gauss_trafo_approx_type(gauss_plan *this, unsigned approx_type)
{
  int j;

  if(approx_type==SAMPLED_FT)
    {
      for(j=0; j<this->n; j++)
	this->b[j] = 1.0/this->p * csqrt(PI/this->delta)*
	  cexp(-PI*PI*(j-this->n/2)*(j-this->n/2)/
	       (this->p*this->p*this->delta));
    }

  if(approx_type==DFT_APPROX)
    {
      this->my_fftw_plan = fftw_plan_dft_1d(this->n, this->b, this->b,
					    FFTW_FORWARD, FFTW_MEASURE);
      for(j=0; j<this->n; j++)
	this->b[j] = cexp(-this->p*this->p*this->delta*(j-this->n/2)*
			  (j-this->n/2)/((double)this->n*this->n)) / this->n;
      
      fft_shift(this->b, this->n);  
      fftw_execute(this->my_fftw_plan);
      fft_shift(this->b, this->n);
      
      fftw_destroy_plan(this->my_fftw_plan);
    }
}

void init_gauss_trafo_nfft_change_m(gauss_plan *this, int m)
{
  int sigma_n,j;

  nfft_finalize(this->my_nfft_plan1);
  nfft_finalize(this->my_nfft_plan2);

  sigma_n=next_power_of_2(2*this->n);
  if(0)
    {
      nfft_init_specific(this->my_nfft_plan1, 1, &(this->n), this->N, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT,
			 FFTW_MEASURE);
      nfft_init_specific(this->my_nfft_plan2, 1, &(this->n), this->M, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI, FFTW_MEASURE);
      this->my_nfft_plan1->x = this->x;
      this->my_nfft_plan2->x = this->y;
    }
  else
    {
      nfft_init_specific(this->my_nfft_plan1, 1, &(this->n), this->N, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT,
			 FFTW_MEASURE);
      nfft_init_specific(this->my_nfft_plan2, 1, &(this->n), this->M, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_X,
			 FFTW_MEASURE);
      for(j=0; j<this->my_nfft_plan1->M; j++)
	this->my_nfft_plan1->x[j] = this->x[j]/this->p;
      for(j=0; j<this->my_nfft_plan2->M; j++)
	this->my_nfft_plan2->x[j] = this->y[j]/this->p;
    }
  
  if(this->my_nfft_plan1->nfft_flags & PRE_PSI)
    nfft_precompute_psi(this->my_nfft_plan1);

  this->my_nfft_plan1->f = this->alpha; 

  if(this->my_nfft_plan2->nfft_flags & PRE_PSI)
    nfft_precompute_psi(this->my_nfft_plan2);

  this->my_nfft_plan2->f_hat = this->my_nfft_plan1->f_hat;

  this->my_nfft_plan2->f = this->f;
}

void init_gauss_trafo_nfft_change_p(gauss_plan *this, double p)
{
  int sigma_n,j;

  nfft_finalize(this->my_nfft_plan1);
  nfft_finalize(this->my_nfft_plan2);

  this->p=p;

  sigma_n=next_power_of_2(2*this->n);
  if(p==1)
    {
      nfft_init_specific(this->my_nfft_plan1, 1, &(this->n), this->N, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT,
			 FFTW_MEASURE);
      nfft_init_specific(this->my_nfft_plan2, 1, &(this->n), this->M, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI, FFTW_MEASURE);
      this->my_nfft_plan1->x = this->x;
      this->my_nfft_plan2->x = this->y;
    }
  else
    {
      nfft_init_specific(this->my_nfft_plan1, 1, &(this->n), this->N, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT,
			 FFTW_MEASURE);
      nfft_init_specific(this->my_nfft_plan2, 1, &(this->n), this->M, &sigma_n,
			 7, PRE_PHI_HUT| PRE_PSI| MALLOC_X,
			 FFTW_MEASURE);
      for(j=0; j<this->my_nfft_plan1->M; j++)
	this->my_nfft_plan1->x[j] = this->x[j]/this->p;
      for(j=0; j<this->my_nfft_plan2->M; j++)
	this->my_nfft_plan2->x[j] = this->y[j]/this->p;
    }
  
  if(this->my_nfft_plan1->nfft_flags & PRE_PSI)
    nfft_precompute_psi(this->my_nfft_plan1);

  this->my_nfft_plan1->f = this->alpha; 

  if(this->my_nfft_plan2->nfft_flags & PRE_PSI)
    nfft_precompute_psi(this->my_nfft_plan2);

  this->my_nfft_plan2->f_hat = this->my_nfft_plan1->f_hat;

  this->my_nfft_plan2->f = this->f;
}
/** @} 
 */

/** @defgroup test gauss transform
 * contains the tests
 * @{ 
 */

/** gauss_andersson_test similar to the test in
 *  F. Andersson and G. Beylkin.
 *  The fast Gauss transform with complex parameters.
 *  J. Comput. Physics, to appear.
 */
void gauss_andersson_test()
{
  gauss_plan my_plan;
  double t;
  fftw_complex *direct;
  int N;

  fftw_complex delta=4*(138+I*100);
  int n=128;
  int r,r_max;
  r_max=1000;

  printf("n=%d, delta=%1.3e+i%1.3e\n",n,creal(delta),cimag(delta));

  for(N=((int)(1U<<6)); N<((int)(1U<<22)); N=N<<1)
    {
      printf("$%d$\t & ",N);
      init_gauss_trafo(&my_plan, N, N, delta, n, 1);
      init_gauss_trafo_approx_type(&my_plan, SAMPLED_FT);
      direct = (fftw_complex*)fftw_malloc(my_plan.M*sizeof(fftw_complex));
      
      SWAPC(direct,my_plan.f);

      if(N<((int)(1U<<19)))
	{
	  t=second();
	  gauss_trafo(&my_plan);
	  t=second()-t;
	  if(t>=0.1)
	    printf("$%1.1e$\t & ",t);
	  else
	    {
	      t=second();
	      for(r=0;r<r_max;r++)
		gauss_trafo(&my_plan);
	      t=second()-t;
	      printf("$%1.1e$\t & ",t/r_max);
	    }
	}
      else
	printf("\t\t & ");	

      if(N<((int)(1U<<14)))
	{
	  my_plan.pre_cexp=(fftw_complex*)fftw_malloc(my_plan.M*my_plan.N*
						      sizeof(fftw_complex));
	  precompute_cexp(&my_plan);
	  t=second();
	  gauss_trafo_pre(&my_plan);
	  t=second()-t;
	  if(t>=0.1)
	    printf("$%1.1e$\t & ",t);
	  else
	    {
	      t=second();
	      for(r=0;r<r_max;r++)
		gauss_trafo_pre(&my_plan);
	      t=second()-t;
	      printf("$%1.1e$\t & ",t/r_max);
	    }
	}
      else
	printf("\t\t & ");	

      t=second();
      fast_gauss_trafo_ndft(&my_plan);
      t=second()-t;
      if(t>=0.1)
	printf("$%1.1e$\t & ",t);
      else
	{
	  t=second();
	  for(r=0;r<r_max;r++)
	    fast_gauss_trafo_ndft(&my_plan);
	  t=second()-t;
	  printf("$%1.1e$\t & ",t/r_max);
	}
      
      SWAPC(direct,my_plan.f);
      
      t=second();
      fast_gauss_trafo_nfft(&my_plan);
      t=second()-t;
      if(t>=0.1)
	printf("$%1.1e$\t & ",t);
      else
	{
	  t=second();
	  for(r=0;r<r_max;r++)
	    fast_gauss_trafo_nfft(&my_plan);
	  t=second()-t;
	  printf("$%1.1e$\t & ",t/r_max);
	}
      
      printf("$%1.1e$\t \\\\ \n",
	     l_infty_error_relative_to_alpha(direct,my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));
      fflush(stdout);

      finalize_gauss_trafo(&my_plan);
      fftw_cleanup();
    }
}

void gauss_error()
{
  gauss_plan my_plan;
  fftw_complex *direct;
  int n;

  fftw_complex delta=4*(138+I*100);
  int N=1000;
  int M=1000;
  
  printf("N=%d;\tM=%d;\ndelta=%1.3e+i*%1.3e;\n",N,M,creal(delta),cimag(delta));
  printf("error=[\n");

  for(n=8; n<=128; n+=4)
    {
      printf("%d\t",n);
      
      init_gauss_trafo(&my_plan, N, M, delta, n, 1);

      direct = (fftw_complex*)fftw_malloc(my_plan.M*sizeof(fftw_complex));
      
      SWAPC(direct,my_plan.f);
      gauss_trafo(&my_plan);
      SWAPC(direct,my_plan.f);

      /** fast gauss trafo with sampled continuous Fourier transform for b */
      init_gauss_trafo_approx_type(&my_plan, SAMPLED_FT);

      fast_gauss_trafo_ndft(&my_plan);
      printf("%1.3e\t",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      init_gauss_trafo_nfft_change_m(&my_plan, 7);
      fast_gauss_trafo_nfft(&my_plan);
      printf("%1.3e\t",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      init_gauss_trafo_nfft_change_m(&my_plan, 3);
      fast_gauss_trafo_nfft(&my_plan);
      printf("%1.3e\t",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));
     
      /** fast gauss trafo with truncated kernel and DFT-coefficients for b */
      init_gauss_trafo_approx_type(&my_plan, DFT_APPROX);
      fast_gauss_trafo_ndft(&my_plan);
      printf("%1.3e\t",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      init_gauss_trafo_nfft_change_m(&my_plan, 7);
      fast_gauss_trafo_nfft(&my_plan);
      printf("%1.3e\t",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      init_gauss_trafo_nfft_change_m(&my_plan, 3);
      fast_gauss_trafo_nfft(&my_plan);
      printf("%1.3e\n",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      finalize_gauss_trafo(&my_plan);
      fftw_cleanup();
    }
  printf("];\n");
}

void gauss_error_p()
{
  gauss_plan my_plan;
  fftw_complex *direct;
  int n;

  fftw_complex delta=20+40I;
  int N=1000;
  int M=1000;
  
  printf("N=%d;\tM=%d;\ndelta=%1.3e+i*%1.3e;\n",N,M,creal(delta),cimag(delta));
  printf("error=[\n");

  for(n=8; n<=128; n+=4)
    {
      printf("%d\t",n);
      
      /** p=1 */
      init_gauss_trafo(&my_plan, N, M, delta, n, 1);

      direct = (fftw_complex*)fftw_malloc(my_plan.M*sizeof(fftw_complex));
      
      SWAPC(direct,my_plan.f);
      gauss_trafo(&my_plan);
      SWAPC(direct,my_plan.f);

      /** fast gauss trafo with sampled continuous Fourier transform for b */
      init_gauss_trafo_approx_type(&my_plan, SAMPLED_FT);
      fast_gauss_trafo_nfft(&my_plan);
      printf("%1.3e\t",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      /** p=1.5 */
      init_gauss_trafo_nfft_change_p(&my_plan, 1.5);
      init_gauss_trafo_approx_type(&my_plan, SAMPLED_FT);
      fast_gauss_trafo_nfft(&my_plan);
      printf("%1.3e\t",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      /** p=2 */
      init_gauss_trafo_nfft_change_p(&my_plan, 2);
      init_gauss_trafo_approx_type(&my_plan, SAMPLED_FT);
      fast_gauss_trafo_nfft(&my_plan);
      printf("%1.3e\n",
	     l_infty_error_relative_to_alpha(direct, my_plan.f, my_plan.M,
					     my_plan.alpha, my_plan.N));

      finalize_gauss_trafo(&my_plan);
      fftw_cleanup();
    }
  printf("];\n");
}
/** @} 
 */

int main()
{
  /**  pipe to output_andersson.tex */
  gauss_andersson_test();

  /** pipe to output_error.m */
  //gauss_error();

  /** pipe to output_error_p.m */
  //gauss_error_p();

  return 1;
}
