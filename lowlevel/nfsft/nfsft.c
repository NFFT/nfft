/** 
 * \file nfsft.c
 * \brief Implementation file for the NFSFT module
 * \author Jens Keiner
 */

/* Include library header. */
#include "nfft3.h"

/* Include utilities header. */
#include "util.h"

//#include "util.h"
//#include "u.h"
//#include "direct.h"
//#include "legendre.h"
//#include "flft.h"
//#include "c2f.h"
//#include <stdlib.h>

#define NFSFT_DEFAULT_NFFT_CUTOFF 6

/** Global structure for wisdom. */
static struct nfsft_wisdom wisdom = {0,0};

void nfsft_init(nfsft_plan *plan, int N, int M)
{
  /* Call nfsft_init_guro with no flags. */
  nfsft_init_advanced(plan, N, M, 0U);
} 
  
void nfsft_init_advanced(nfsft_plan* plan, int N, int M, unsigned int nfsft_flags)
{
  /* Call nfsft_init_guro with no flags and default NFFT cut-off. */
  nfsft_init_guru(plan, N, M, 0U, NFSFT_DEFAULT_NFFT_CUTOFF);
}

void nfsft_init_guru(nfsft_plan* plan, int N, int M, unsigned int flags, 
                     int nfft_cutoff)
{
  /** Array for NFFT sizes */
  int nfft_size[2] = {0,0};
  /** Array for FFTW sizes */
  int fftw_size[2] = {0,0};

  /** Allocate memory for the plan. */
  plan = (nfsft_plan*) malloc(sizeof(struct nfsft_plan));
  
  /* Save the flags in the plan. */
  plan->flags = flags;
  
  /* Save the bandwidth N and the number of samples M in the plan. */
  plan->N = N;
  plan->M_total = M;
  
  /* Calculate the next greater power of two with respect to the bandwidth N and
   * the corresponding exponent. */
  plan->NPT = next_power_of_2(plan->N);
  plan->t = log((double)plan->NPT)/log(2.0);
  
  /* Save length of array of Fourier coefficients. Owing to the data layout the 
   * length is (2N+1)(2NPT) */
  plan->N_total = (2*plan->N+1)*2*plan->NPT;
  
  /* Allocate memory for spherical Fourier coefficients, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_F_HAT)
  {
    plan->f_hat = (complex*) malloc(plan->N_total*sizeof(complex));
  }
  
  /* Allocate memory for samples, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_F)
  {
    plan->f = (complex*) malloc(plan->M_total*sizeof(complex));
  }
  
  /* Allocate memory for nodes, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_X)
  {
    plan->x = (double*) malloc(2*plan->M_total*sizeof(double));
  }
  
  /* Initialize NFFT plan. */
  /** \todo Check NFFT and FFTW sizes. */
  nfft_size[0] = 2*(plan->NPT+1);
  nfft_size[1] = 2*(plan->NPT+1);
  fftw_size[0] = 4*plan->NPT;
  fftw_size[1] = 4*plan->NPT;
  /** \todo Check NFFT flags. */
  nfft_init_guru(&plan->plan_nfft, 2, nfft_size, plan->M, fftw_size, 
                     nfft_cutoff, PRE_PHI_HUT | PRE_PSI | FFT_OUT_OF_PLACE | 
                     FFTW_INIT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  /* Assign angle array. */
  plan->plan_nfft.x = plan->x;
  /* Assign result array. */
  plan->plan_nfft.f = plan->f;
  /* Precompute PSI, if necessary. */
  if (plan->plan_nfft.nfft_flags & PRE_PSI) 
  {  
    nfft_precompute_psi(&plan->plan_nfft); 
  } 
}

void nfsft_precompute(int M, double threshold, nfsft_precompute_flags_old flags)
{
  int i;
  int ti;
  
  /*  Check if already initialized. */
  if (wisdom.initialized == 1 )
  {
    return;
  }
  
  /* Set the flags. */
  wisdom.flags = flags;

  /* Set t. */
  wisdom.t = ngpt(M);
  /* Set N. */
  wisdom.N = 1<<wisdom.t;

  /* Precompute three-term recurrence coefficients. */
  wisdom.alpha = (double*) malloc((BW_MAX+1)*(BW_MAX+1)*sizeof(double));    
  wisdom.beta = (double*) malloc((BW_MAX+1)*(BW_MAX+1)*sizeof(double));    
  wisdom.gamma = (double*) malloc((BW_MAX+1)*(BW_MAX+1)*sizeof(double));    
  alpha_al_all_old(wisdom.alpha,BW_MAX);
  beta_al_all_old(wisdom.beta,BW_MAX);
  gamma_al_all_old(wisdom.gamma,BW_MAX);
  
  /* Set the threshold. */
  wisdom.threshold  = threshold;

  /* Precompute. */
  wisdom.U = precomputeU_old(wisdom.t, wisdom.threshold, wisdom.alpha, wisdom.beta, 
                         wisdom.gamma, 
                         (wisdom.flags & NFSFT_BW_WINDOW_OLD) == NFSFT_BW_WINDOW_OLD);
  
  /* Delete coefficients direct algorithms are deactivated. */
  if ((wisdom.flags & NFSFT_FAST_ONLY_OLD) == NFSFT_FAST_ONLY_OLD)
  {
    free(wisdom.alpha);
    free(wisdom.beta);
    free(wisdom.gamma);
  }
  
  /* Check, if bandwidth big enough. */
  if (wisdom.N >= 4)
  {  
    wisdom.work       = (complex*) calloc(3*(wisdom.N+1),sizeof(complex));   
    wisdom.old        = (complex*) malloc(2*(wisdom.N)*sizeof(complex));
    wisdom.ergeb      = (complex*) calloc(3*(wisdom.N+1),sizeof(complex));
    wisdom.vec1       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec2       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec3       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec4       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.a2         = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.b2         = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.z          = (complex*) malloc(wisdom.N*sizeof(complex));
    wisdom.plans_dct3 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(wisdom.t-1));
    wisdom.plans_dct2 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(wisdom.t-1));
    wisdom.kinds      = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    wisdom.kinds[0]   = FFTW_REDFT01;
    wisdom.kinds[1]   = FFTW_REDFT01;
    wisdom.kindsr     = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    wisdom.kindsr[0]  = FFTW_REDFT10;
    wisdom.kindsr[1]  = FFTW_REDFT10;
    wisdom.lengths    = (int*) malloc((wisdom.t-1)*sizeof(int));
    for (i = 0, ti = 4; i < wisdom.t-1; i++, ti<<=1)
    {
      wisdom.lengths[i] = ti;
      wisdom.plans_dct3[i] = 
        fftw_plan_many_r2r(1, &wisdom.lengths[i], 2, (double*)wisdom.vec1, NULL, 
                           2, 1, (double*)wisdom.vec1, NULL, 2, 1, wisdom.kinds, 
                           0);
      wisdom.plans_dct2[i] = 
        fftw_plan_many_r2r(1, &wisdom.lengths[i], 2, (double*)wisdom.vec1, NULL, 
                           2, 1, (double*)wisdom.vec1, NULL, 2, 1,wisdom.kindsr, 
                           0);
    }  
  }
  
  wisdom.initialized = 1;
}

void nfsft_forget()
{
  int i;
  
  if (wisdom.initialized == 1)
  {
    if (wisdom.N >= 4)
    {  
      forgetU_old(wisdom.U,wisdom.N,wisdom.t,
              (wisdom.flags & NFSFT_BW_WINDOW_OLD) == NFSFT_BW_WINDOW_OLD);
      
      fftw_free(wisdom.work);
      fftw_free(wisdom.old);
      fftw_free(wisdom.ergeb);
      fftw_free(wisdom.vec1);
      fftw_free(wisdom.vec2);
      fftw_free(wisdom.vec3);
      fftw_free(wisdom.vec4);
      fftw_free(wisdom.a2);
      fftw_free(wisdom.b2);
      free(wisdom.z);
      
      for(i = 0; i < wisdom.t-1; i++)
      {
        fftw_destroy_plan(wisdom.plans_dct3[i]);
        fftw_destroy_plan(wisdom.plans_dct2[i]);
      }  
      
      free(wisdom.plans_dct3);
      free(wisdom.plans_dct2);
      free(wisdom.kinds);
      free(wisdom.kindsr);
      free(wisdom.lengths);
    }
    if ((wisdom.flags & NFSFT_FAST_ONLY_OLD) != NFSFT_FAST_ONLY_OLD)
    {
      free(wisdom.alpha);
      free(wisdom.beta);
      free(wisdom.gamma);
    }
    wisdom.initialized = 0;
  }
}
                     
void nfsft_finalize(nfsft_plan_old plan)
{
  nfft_finalize(&plan->plan_nfft);  
  free(plan);
}

void ndsft_trafo(nfsft_plan_old plan)
{
  /* Compute direct transform. */
  if (wisdom.flags & NFSFT_FAST_ONLY_OLD)
  {  
  }
  else
  {
    /** Check for normalization. */
    if (plan->flags & NFSFT_NORMALIZED_OLD)
    {
      normalize_f_hat_old(plan->f_hat, plan->M);
    }  
    ndsft_old(plan->D, plan->angles, plan->f, plan->M, plan->f_hat, &wisdom);
  }
}

void ndsft_adjoint(nfsft_plan_old plan)
{ 
  if (wisdom.flags & NFSFT_FAST_ONLY_OLD)
  {  
  }
  else
  {
    adjoint_ndsft_old(plan->D, plan->angles, plan->f, plan->M, plan->f_hat, 
                  &wisdom);
    /** Check for normalization. */
    if (plan->flags & NFSFT_NORMALIZED_OLD)
    {
      normalize_f_hat_old(plan->f_hat, plan->M);
    }  
  }
}

void nfsft_trafo(nfsft_plan_old plan)
{
  /** Counter for loops */
  int i,n;
	
  if (wisdom.initialized == 0 || plan->M > wisdom.N)
  {
    return;
  }  
  
  if (plan->M < 3)
  {
    ndsft_trafo_old(plan);
  }
  else if ((wisdom.flags & NFSFT_BW_WINDOW_OLD) == 0U || plan->M > 1<<(wisdom.t-1))
  {
    plan->plan_nfft.f = plan->f;
    
    /** Check for normalization. */
    if (plan->flags & NFSFT_NORMALIZED_OLD)
    {
      normalize_f_hat_old(plan->f_hat, plan->M);
    }  

    /* Compute FLFT. */
    for (n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {
//		  fprintf(stderr,"flft: n = %d\n",n);
//			fflush(stderr);
      flft_old(plan->M, plan->t, abs(n), plan->f_hat[i], &wisdom,&nstab,&ntotal);
    }
   
    /* Convert Chebyshev coefficients to Fourier coefficients. */
//    fprintf(stderr,"cheb2exp_old\n",n);
//	  fflush(stderr);
    cheb2exp_old(plan->plan_nfft.f_hat, plan->f_hat, plan->M, plan->N); 

		//norm = norm_nfft_1_old(plan->plan_nfft.f_hat, 2*(plan->N+1));
		//fprintf(stderr,"nfft_norm_1 = %.4E\n",norm);
		fflush(stderr);
		
    /* Execute NFFT. */
		if (plan->flags & NFSFT_USE_NDFT_OLD)
		{
      ndft_trafo(&plan->plan_nfft);
		}
		else
		{
      nfft_trafo(&plan->plan_nfft);
		}
		
  } 
}

void nfsft_adjoint(nfsft_plan_old plan)
{
  /** Counter for loops */
  int i,n;
    
  if (wisdom.initialized == 0 || plan->M > wisdom.N)
  {
    return;
  }
  
  if (plan->M < 3)
  {
    ndsft_adjoint_old(plan);
  }
  else if ((wisdom.flags & NFSFT_BW_WINDOW_OLD) == 0U || plan->M > 1<<(wisdom.t-1))
  {
    plan->plan_nfft.f = plan->f;

    /* Execute adjoint NFFT. */
		if (plan->flags & NFSFT_USE_NDFT_OLD)
		{
      ndft_adjoint(&plan->plan_nfft);
		}
		else
		{
      nfft_adjoint(&plan->plan_nfft);
		}
    
    /* Convert Chebyshev coefficients to Fourier coefficients. */
    cheb2exp_adjoint_old(plan->plan_nfft.f_hat, plan->f_hat, plan->M, plan->N); 
        
    /* Calculate FLFT. */
    for(n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {
      flft_adjoint_old(plan->M, plan->t, abs(n), plan->f_hat[i], &wisdom);
    }          

    /** Check for normalization. */
    if (plan->flags & NFSFT_NORMALIZED_OLD)
    {
      normalize_f_hat_old(plan->f_hat, plan->M);
    }  
  }    
}
