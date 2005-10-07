/** 
* \file nfsft.c
* \brief Implementation file for the nfsft library
* \author Jens Keiner
*/

#include "api.h"
#include "util.h"
#include "u.h"
#include "direct.h"
#include "legendre.h"
#include "flft.h"
#include "c2f.h"
#include <stdlib.h>

/** Global structure for wisdom. */
static struct nfsft_wisdom_old wisdom = {0,0};
static int nstab;
static int ntotal;

void nfsft_precompute_old(int M, double threshold, nfsft_precompute_flags_old flags)
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

void nfsft_forget_old()
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

nfsft_plan_old nfsft_init_old(int m, int d, fftw_complex **f_hat, double *angles, 
                      fftw_complex *f, nfsft_flags_old flags)
{
  return nfsft_init_guru_old(m, d, f_hat, angles, f, flags, 6);
} 
  
nfsft_plan_old nfsft_init_guru_old(int m, int d, fftw_complex **f_hat, double *angles, 
                           fftw_complex *f, nfsft_flags_old flags, int cutoff)
{
  int nfft_size[2] = {0,0};
  int fftw_size[2] = {0,0};

  nfsft_plan_old plan = (nfsft_plan_old) malloc(sizeof(struct nfsft_plan_s_old));
  
  plan->flags = flags;
  plan->D = d;
  plan->M = m;
  plan->t = ngpt(plan->M);
  plan->N = 1<<plan->t;
  plan->angles = angles;
  plan->f_hat = f_hat;
  plan->f = f;
  
  /* Initialize NFFT. */
  nfft_size[0] = 2*(plan->N+1);
  nfft_size[1] = 2*(plan->N+1);
  fftw_size[0] = 4*plan->N;
  fftw_size[1] = 4*plan->N;
	/* MS */
/*  nfft_init_guru(&plan->plan_nfft, 2, nfft_size, plan->D, fftw_size, 
                     cutoff, PRE_PHI_HUT | PRE_PSI | MALLOC_F_HAT | FFT_OUT_OF_PLACE, 
                     FFTW_ESTIMATE| FFTW_DESTROY_INPUT);*/
  nfft_init_guru(&plan->plan_nfft, 2, nfft_size, plan->D, fftw_size, 
                     cutoff, PRE_PHI_HUT | PRE_PSI | MALLOC_F_HAT | FFT_OUT_OF_PLACE | FFTW_INIT, 
                     FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
	/* MS */
  /* Assign angle array. */
  plan->plan_nfft.x = plan->angles;
  /* Assign result array. */
  plan->plan_nfft.f = plan->f;
  if (plan->plan_nfft.nfft_flags & PRE_PSI) 
  {  
    nfft_precompute_psi(&plan->plan_nfft); 
  } 
  
  return(plan);
}

void nfsft_finalize_old(nfsft_plan_old plan)
{
  nfft_finalize(&plan->plan_nfft);  
  free(plan);
}

void ndsft_trafo_old(nfsft_plan_old plan)
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

void ndsft_adjoint_old(nfsft_plan_old plan)
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

void nfsft_trafo_old(nfsft_plan_old plan)
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

void nfsft_adjoint_old(nfsft_plan_old plan)
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


int nfsft_trafo_old_stab(nfsft_plan_old plan)
{
  /** Counter for loops */
  int i,n;
  
  if (wisdom.initialized == 0 || plan->M > wisdom.N)
  {
    return 0;
  }
  
  nstab = 0; 
  ntotal = 0;

  if (plan->M < 3)
  {
    return 0;
  }
  else
  {
    /* Compute FLFT. */
    for (n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {
      flft_stab_old(plan->M, plan->t, abs(n), plan->f_hat[i], &wisdom,&nstab,&ntotal);
    }
  } 
  return nstab;
}

nfsft_plan_old nfsft_init_stab(int m, int d, fftw_complex **f_hat, double *angles, 
                      fftw_complex *f, nfsft_flags_old flags)
{
  nfsft_plan_old plan = (nfsft_plan_old) malloc(sizeof(struct nfsft_plan_s_old));
  
  plan->M = m;
  plan->t = ngpt(plan->M);
  plan->N = 1<<plan->t;
  plan->f_hat = f_hat;
  return(plan);
}

void nfsft_precompute_stab(int M, double threshold)
{
  int i;
  int ti;
  
  /*  Check if already initialized. */
  if (wisdom.initialized == 1 )
  {
    return;
  }
  
  wisdom.t = ngpt(M);
  wisdom.N = 1<<wisdom.t;
  
  /* Precompute three-term recurrence coefficients. */
  wisdom.alpha = (double*) malloc((BW_MAX+1)*(BW_MAX+1)*sizeof(double));    
  wisdom.beta = (double*) malloc((BW_MAX+1)*(BW_MAX+1)*sizeof(double));    
  wisdom.gamma = (double*) malloc((BW_MAX+1)*(BW_MAX+1)*sizeof(double));    
  //wisdom.gamma_m1 = (double*) malloc((BW_MAX+1)*sizeof(double));    
  alpha_al_all_old(wisdom.alpha,BW_MAX);
  beta_al_all_old(wisdom.beta,BW_MAX);
  gamma_al_all_old(wisdom.gamma,BW_MAX);
  //gamma_al_m1_all(wisdom.gamma_m1,BW_MAX);
  
  if (wisdom.N >= 4)
  {  
    wisdom.threshold  = threshold;
    wisdom.work       = (complex*) calloc(3*(wisdom.N+1),sizeof(complex));   
    wisdom.old        = (complex*) malloc(2*(wisdom.N)*sizeof(complex));
    wisdom.ergeb      = (complex*) calloc(3*(wisdom.N+1),sizeof(complex));
    wisdom.vec1       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec2       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec3       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.vec4       = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.a2         = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.b2         = (complex*) fftw_malloc(sizeof(complex)*(wisdom.N+1));
    wisdom.kinds      = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    wisdom.kinds[0]   = FFTW_REDFT01;
    wisdom.kinds[1]   = FFTW_REDFT01;
    wisdom.kindsr     = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    wisdom.kindsr[0]  = FFTW_REDFT10;
    wisdom.kindsr[1]  = FFTW_REDFT10;
    wisdom.plans_dct3 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(wisdom.t-1));
    wisdom.plans_dct2 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(wisdom.t-1));
    wisdom.lengths    = (int*) malloc((wisdom.t-1)*sizeof(int));
    for (i = 0, ti = 4; i < wisdom.t-1; i++, ti<<=1)
    {
      wisdom.lengths[i] = ti;
      wisdom.plans_dct3[i] = fftw_plan_many_r2r(1, &wisdom.lengths[i], 2, (double*)wisdom.vec1, NULL, 2, 1,
                                                (double*)wisdom.vec1, NULL, 2, 1, wisdom.kinds, 0);
      wisdom.plans_dct2[i] = fftw_plan_many_r2r(1, &wisdom.lengths[i], 2, (double*)wisdom.vec1, NULL, 2, 1,
                                                (double*)wisdom.vec1, NULL, 2, 1, wisdom.kindsr, 0);
    }  
    wisdom.z = (complex*) malloc(wisdom.N*sizeof(complex));
    wisdom.U = precomputeU_stab_old(wisdom.t, wisdom.threshold, wisdom.alpha, wisdom.beta, wisdom.gamma);
  }
  
  wisdom.initialized = 1;
}

void nfsft_forget_stab_old()
{
  int i;
  
  if (wisdom.initialized == 1)
  {
    if (wisdom.N >= 4)
    {  
      forgetU_stab_old(wisdom.U,wisdom.N,wisdom.t);
      
      free(wisdom.z);
      fftw_free(wisdom.work);
      fftw_free(wisdom.old);
      fftw_free(wisdom.ergeb);
      fftw_free(wisdom.vec1);
      fftw_free(wisdom.vec2);
      fftw_free(wisdom.vec3);
      fftw_free(wisdom.vec4);
      fftw_free(wisdom.a2);
      fftw_free(wisdom.b2);
      
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
    free(wisdom.alpha);
    free(wisdom.beta);
    free(wisdom.gamma);
    //free(wisdom.gamma_m1);
    wisdom.initialized = 0;
  }
}

void nfsft_finalize_stab_old(nfsft_plan_old plan)
{
  free(plan);
}

void nfsft_get_stat(int *nstabp, int *ntotalp)
{
  *nstabp = nstab;
  *ntotalp = ntotal;
}

void nfsft_flft_trafo(complex *in,int n, int M)
{
  flft_old(M, ngpt(M), abs(n), in, &wisdom, &nstab, &ntotal);
}

void nfsft_flft_trafo_adjoint(complex *in,int n, int M)
{ 
  flft_adjoint_old(M, ngpt(M), abs(n), in, &wisdom);
}

