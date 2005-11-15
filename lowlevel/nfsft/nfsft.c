/** 
 * \file nfsft.c
 * \brief Implementation file for the NFSFT module
 * \author Jens Keiner
 */

/* Include library header. */
#include "nfft3.h"

/* Include private API header. */
#include "api.h"

/* Include utilities header. */
#include "util.h"

/* Include standard C headers. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

//#include "util.h"
//#include "u.h"
//#include "direct.h"
#include "legendre.h"
#include "dpt.h"
//#include "flft.h"
//#include "c2f.h"
//#include <stdlib.h>

#define NFSFT_DEFAULT_NFFT_CUTOFF 6

/** Global structure for wisdom. */
static struct nfsft_wisdom wisdom = {false,0U};

void nfsft_init(nfsft_plan *plan, int N, int M)
{
  /* Call nfsft_init_advanced with flags to allocate memory. */
  nfsft_init_advanced(plan, N, M, NFSFT_MALLOC_X | NFSFT_MALLOC_F | 
                      NFSFT_MALLOC_F_HAT);
} 
  
void nfsft_init_advanced(nfsft_plan* plan, int N, int M, 
                         unsigned int flags)
{
  /* Call nfsft_init_guro with the flags and default NFFT cut-off. */
  nfsft_init_guru(plan, N, M, flags, NFSFT_DEFAULT_NFFT_CUTOFF);
}

void nfsft_init_guru(nfsft_plan* plan, int N, int M, unsigned int flags, 
                     int nfft_cutoff)
{
  /** Array for NFFT sizes */
  int *nfft_size;
  /** Array for FFTW sizes */
  int *fftw_size;

  /* Save the flags in the plan. */
  plan->flags = flags;
  
  /* Save the bandwidth N and the number of samples M in the plan. */
  plan->N = N;
  plan->M_total = M;
  
  /* Calculate the next greater power of two with respect to the bandwidth N and
   * the corresponding exponent. */
  plan->NPT = next_power_of_2(plan->N);
  plan->t = (int)(log((double)plan->NPT)/log(2.0));
  
  /* Save length of array of Fourier coefficients. Owing to the data layout the 
   * length is (2N+1)(2NPT) */
  plan->N_total = (2*plan->N+1)*(2*plan->NPT+1);
  
  /* Allocate memory for spherical Fourier coefficients, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_F_HAT)
  {
    plan->f_hat = (complex*) calloc(plan->N_total,sizeof(complex));
  }
  
  /* Allocate memory for samples, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_F)
  {
    plan->f = (complex*) calloc(plan->M_total,sizeof(complex));
  }
  
  /* Allocate memory for nodes, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_X)
  {
    plan->x = (double*) calloc(2*plan->M_total,sizeof(double));
  }
  
  /* Check if fast algorithm is activated. */
  if (plan->flags & NFSFT_NO_FAST_ALGORITHM)
  {
  }
  else
  {
    /* Initialize NFFT plan. */
    nfft_size = (int*) malloc(2*sizeof(int));
    fftw_size = (int*) malloc(2*sizeof(int));
    
    /** TODO Check NFFT and FFTW sizes. */
    nfft_size[0] = 2*(plan->NPT+1);
    nfft_size[1] = 2*(plan->NPT+1);
    fftw_size[0] = 4*plan->NPT;
    fftw_size[1] = 4*plan->NPT;
    /** \todo Check NFFT flags. */
    nfft_init_guru(&plan->plan_nfft, 2, nfft_size, plan->M_total, fftw_size, 
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
    
    free(nfft_size);
    free(fftw_size);
  }
}

void nfsft_precompute(int N, double kappa, 
                      unsigned int flags)
{ 
  /*  Check if already initialized. */
  if (wisdom.initialized == true)
  {
    return;
  }
  
  /* Save the precomputation flags. */
  wisdom.flags = flags;

  /* Compute and save /f$N_{\text{max}} = 2^{\ceil{\log_2 N}}/f$ as next greater 
   * power of two with respect to /f$N/f$. */
  wisdom.N_MAX = next_power_of_2(N);
  /* Save t. */
  wisdom.T_MAX = (int) log((double)wisdom.N_MAX)/log(2.0);

  /* Check, if precomputation for direct algorithms needs to be performed. */    
  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM)
  {
  }
  else
  {
    /*fprintf(stdout,"Precomputation for ndsft_trafo done.\n");*/
    /* Precompute three-term recurrence coefficients. */
    wisdom.alpha = (double*) malloc((wisdom.N_MAX+1)*(wisdom.N_MAX+2)*
      sizeof(double));    
    wisdom.beta = (double*) malloc((wisdom.N_MAX+1)*(wisdom.N_MAX+2)*
      sizeof(double));    
    wisdom.gamma = (double*) malloc((wisdom.N_MAX+1)*(wisdom.N_MAX+2)*
      sizeof(double));
    /** \todo Change to functions which compute only for fixed order n. */       
    alpha_al_all(wisdom.alpha,wisdom.N_MAX+1);
    beta_al_all(wisdom.beta,wisdom.N_MAX+1);
    gamma_al_all(wisdom.gamma,wisdom.N_MAX+1);
  }
  
  /* Check, if precomputation for direct algorithms needs to be performed. */    
  if (wisdom.flags & NFSFT_NO_FAST_ALGORITHM)
  {
  }
  else
  { 
    /* Set the threshold. */
    //wisdom.threshold  = threshold;
  
    /* Precompute. */
    //wisdom.U = precomputeU_old(wisdom.t, wisdom.threshold, wisdom.alpha, wisdom.beta, 
    //                      wisdom.gamma, 
    //                       (wisdom.flags & NFSFT_BW_WINDOW_OLD) == NFSFT_BW_WINDOW_OLD);
    
    /* Delete coefficients direct algorithms are deactivated. */
    /*if ((wisdom.flags & NFSFT_FAST_ONLY_OLD) == NFSFT_FAST_ONLY_OLD)
    {
      free(wisdom.alpha);
      free(wisdom.beta);
      free(wisdom.gamma);
    }*/
    
    /* Check, if bandwidth big enough. */
    /*if (wisdom.N >= 4)
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
    }
    
    wisdom.initialized = 1;*/
  }

  /* Wisdom has been initialised. */
  wisdom.initialized = true;
}
 
void nfsft_forget()
{
  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM == 0)
  {
    /* Free arrays holding three-term recurrence coefficients. */
    free(wisdom.alpha);   
    free(wisdom.beta);    
    free(wisdom.gamma);
  }

  /* Wisdom is now uninitialised. */
  wisdom.initialized = false;
}
                     
void nfsft_finalize(nfsft_plan* plan)
{
}

void ndsft_trafo(nfsft_plan* plan)
{
  /** Node index */
  int m;
  /** Legendre index k */
  int k;
  /** Legendre index n */
  int n;
  /** n_abs = |n| */
  int n_abs;
  
  /**
   * Pointer to current three-term recurrence coefficient \f$\alpha_k^n\f$
   * for associated Legendre functions \f$P_k^n\f$. 
   */
  double *alpha;
  /**
   * Pointer to current three-term recurrence coefficient \f$\beta_k^n\f$
   * for associated Legendre functions \f$P_k^n\f$. 
   */
  double *gamma;
  
  /** Index used in Clenshaw algorithm. */
  int index;
  /** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *a;
  /** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *temp;
  
  /** The final result for a single node. */
  complex f_m;
  
  /** Used to store the angles theta_m */
  double *theta;
  /** Used to store the angles phi_m */
  double *phi;
  
  /* Allocate memory for auxilliary arrays. */
  temp = (complex*) malloc (sizeof(complex) * (plan->N+1));
  theta = (double*) malloc (sizeof(double) * plan->M_total);    
  phi = (double*) malloc (sizeof(double) * plan->M_total);    
  
  if (plan->flags & NFSFT_NORMALIZED)
  {
    for (k = 0; k <= plan->N; k++)
    {
      for (n = -k; n <= k; n++)
      {
        plan->f_hat[(2*plan->NPT+1)*(n+plan->N)+plan->NPT+k] *= 
          sqrt((2*k+1)/(4.0*PI));
      }
    }
  }
   
   /* Distinguish by bandwidth M. */
  if (plan->N == 0)
  {
    /* Constant function */
    for (m = 0; m < plan->M_total; m++)
    {
      plan->f[m] = plan->f_hat[plan->NPT+0];
    }
  }
  else
  { 
    /* Scale angles phi_m from [-0.5,0.5] to [-pi,pi] and angles \theta_m from 
     * [0,0.5] to and [0,pi], respectively. */ 
    for (m = 0; m < plan->M_total; m++)
    {
      theta[m] = 2.0*PI*plan->x[2*m];
      phi[m] = 2.0*PI*plan->x[2*m+1];
    }

    /* Apply cosine to angles theta_m. */
    for (m = 0; m < plan->M_total; m++)
    {
       theta[m] = cos(theta[m]);
    }
    
    /* Evaluate 
     * \[
     *   \sum_{k=0}^N \sum_{n=-k}^k a_k^n P_k^{|n|}(\cos\theta_m) e^{i n \phi_m} 
     *     = \sum_{n=-N}^N \sum_{k=|n|}^N a_k^n P_k^{|n|}(\cos\theta_m) 
     *     e^{i n \phi_m}. 
     * \] */
    for (m = 0; m < plan->M_total; m++)
    {
      /* Initialize result for current node. */
      f_m = 0.0;

      /* For n = -N,...,N, evaluate 
       * \[
       *   b_n := \sum_{k=|n|}^N a_k^n P_k^{|n|}(\cos\vartheta_m)
       * using Clenshaw's algorithm. */
      for (n = -plan->N; n <= plan->N; n++)
      {
        /* Get Fourier coefficients vector. */
        a = &(plan->f_hat[(n+plan->N)*(2*plan->NPT+1)+plan->NPT]);
        
        /* Take absolute value of n. */
        n_abs = abs(n);
        
        /* Get three-term recurrence coefficients vectors. */
        //alpha = wisdom.alpha + plan->N - n_abs;
        //gamma = wisdom.gamma + plan->N - n_abs;
        alpha = &(wisdom.alpha[ROWK(n_abs)]);
        gamma = &(wisdom.gamma[ROWK(n_abs)]);
         
        /* Make copy of array a. */ 
        memcpy(temp,a,(plan->N+1)*sizeof(complex));
        
        /* Clenshaw's algorithm */        
        for (k = plan->N; k > n_abs + 1; k--)
        {
          index = k - n_abs;
          temp[k-1] += temp[k] * alpha[index] * theta[m]; 
          temp[k-2] += temp[k] * gamma[index];
        }
        //fprintf(stdout,"\n");
        
        /* Compute final step if neccesary. */
        if (n_abs < plan->N)
        {  
          temp[0+n_abs] += temp[1+n_abs] * wisdom.alpha[ROWK(n_abs)+1] * theta[m];
        }
        
        /* Write final result b_n of multiplication by normalization constant to 
         * array b  = (b_{-M},...,b_M). */
        f_m += temp[0+n_abs] * wisdom.gamma[ROW(n_abs)] *
          pow(1- theta[m] * theta[m], 0.5*n_abs) * cexp(I*n*phi[m]);
      }
            
      /* Write result to vector f. */
      plan->f[m] = f_m;                
    }
  }  
  
  /* Free auxilliary arrays. */
  free(phi);
  free(theta);
  free(temp);
}

void ndsft_adjoint(nfsft_plan* plan)
{ 
  /** Node index */
  int m;
   /** Legendre index k */
  int k;
  /** Legendre index n */
  int n;
  /** nleg = |n| */
  int n_abs;
  
  /**
   * Pointer to current three-term recurrence coefficient \f$\alpha_k^n\f$
   * for associated Legendre functions \f$P_k^n\f$. 
   */
  double *alpha;
  /**
   * Pointer to current three-term recurrence coefficient \f$\beta_k^n\f$
   * for associated Legendre functions \f$P_k^n\f$. 
   */
  double *gamma;
  
  /** Index used in Clenshaw algorithm. */
  int index;
  /** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *temp;
   
  /** Used to store the angles theta_d */
  double *theta;
  /** Used to store the angles phi_d */
  double *phi;
  
  /* Allocate memory for auxilliary arrays. */
  temp = (complex*) malloc (sizeof(complex) * (plan->N+1));
  theta = (double*) malloc (sizeof(double) * plan->M_total);    
  phi = (double*) malloc (sizeof(double) * plan->M_total);    
  
  /* Scale angles phi_j from [-0.5,0.5) to [-pi,pi) and angles \theta_j from 
   * [0,0.5] to and [0,pi], respectively. */ 
  for (m = 0; m < plan->M_total; m++)
  {
    theta[m] = 2.0*PI*plan->x[2*m];
    phi[m] = 2.0*PI*plan->x[2*m+1];
  }

  for (n = -plan->N; n <= plan->N; n++)
  {
    for (k = abs(n); k <= plan->N; k++)
    {
      plan->f_hat[(2*plan->NPT+1)*(n+plan->N)+plan->NPT+k] = 0.0;
    }  
  }
  
  /* Distinguish by bandwidth N. */
  if (plan->N == 0)
  {
    /* Constant function */
    for (m = 0; m < plan->M_total; m++)
    {
      plan->f_hat[1] += plan->f[m];
    }
  }
  else
  { 
    /* Apply cosine to angles theta_j. */
    for (m = 0; m < plan->M_total; m++)
    {
       theta[m] = cos(theta[m]);
    }
    
    for (m = 0; m < plan->M_total; m++)
    {
      for (n = -plan->N; n <= plan->N; n++)
      {
        /* Take absolute value of n. */
        n_abs = abs(n);

        /* Get three-term recurrence coefficients vectors. */
        alpha = &(wisdom.alpha[ROWK(n_abs)]);
        gamma = &(wisdom.gamma[ROWK(n_abs)]);
        
        temp[n_abs] = plan->f[m] * cexp(-I*n*phi[m]) * wisdom.gamma[ROW(n_abs)] *
          pow(1- theta[m] * theta[m], 0.5*n_abs);
                
        /* Compute final step if neccesary. */
        if (n_abs < plan->N)
        {  
          temp[n_abs+1] = temp[n_abs] * wisdom.alpha[ROWK(n_abs)+1] * theta[m];
        }

        /* Clenshaw's algorithm */        
        for (k = n_abs+2; k <= plan->N; k++)
        {
          index = k - n_abs;
          temp[k] = alpha[index] * theta[m] * temp[k-1] + gamma[index] * temp[k-2];
        }
        
        /* Copy result */
        for (k = n_abs; k <= plan->N; k++)
        {
          plan->f_hat[(n+plan->N)*(2*plan->NPT+1)+plan->NPT+k] += temp[k];
        }  
      }      
    }
  }  
  
  if (plan->flags & NFSFT_NORMALIZED)
  {
    for (k = 0; k <= plan->N; k++)
    {
      for (n = -k; n <= k; n++)
      {
        plan->f_hat[(2*plan->NPT+1)*(n+plan->N)+plan->NPT+k] *= 
          sqrt((2*k+1)/(4.0*PI));
      }
    }
  }
    
  /* Free auxilliary arrays. */
  free(phi);
  free(theta);
  free(temp);
}

void nfsft_trafo(nfsft_plan* plan)
{
  /** Counter for loops */
  /*int i,n;
  
  if (wisdom.initialized == 0 || plan->M > wisdom.N)
  {
    return;
  }  
  
  if (plan->M <= 3)
  {
    ndsft_trafo(plan);
  }
  else if ((wisdom.flags & NFSFT_BW_WINDOW) == 0U || plan->M > 1<<(wisdom.t-1))
  {
    plan->plan_nfft.f = plan->f;
    */
    /** Check for normalization. */
/*    if (plan->flags & NFSFT_NORMALIZED_OLD)
    {
      normalize_f_hat_old(plan->f_hat, plan->M);
    }*/  

    /* Compute FLFT. */
/*    for (n = -plan->M, i = 0; n <= plan->M; n++, i++) 
    {*/
//      fprintf(stderr,"flft: n = %d\n",n);
//      fflush(stderr);
/*      flft_old(plan->M, plan->t, abs(n), plan->f_hat[i], &wisdom,&nstab,&ntotal);
    }*/
   
    /* Convert Chebyshev coefficients to Fourier coefficients. */
//    fprintf(stderr,"cheb2exp_old\n",n);
//    fflush(stderr);
    //cheb2exp_old(plan->plan_nfft.f_hat, plan->f_hat, plan->M, plan->N); 

    //norm = norm_nfft_1_old(plan->plan_nfft.f_hat, 2*(plan->N+1));
    //fprintf(stderr,"nfft_norm_1 = %.4E\n",norm);
    //fflush(stderr);
    
    /* Execute NFFT. */
    /*if (plan->flags & NFSFT_USE_NDFT_OLD)
    {
      ndft_trafo(&plan->plan_nfft);
    }
    else
    {
      nfft_trafo(&plan->plan_nfft);
    }*/
    
  //} 
}

void nfsft_adjoint(nfsft_plan* plan)
{
}
