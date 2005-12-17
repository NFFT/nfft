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
#include <string.h>

#define NFSFT_DEFAULT_NFFT_CUTOFF 12
#define NFSFT_DEFAULT_THRESHOLD 1000
#define NFSFT_BREAKTHROUGH 5

/** Global structure for wisdom. */
static struct nfsft_wisdom wisdom = {false,0U};

inline void c2e(nfsft_plan *plan)
{
  int k;             /**< Degree \f$k\f$                           */
  int n;             /**< Order \f$k\f$                            */
  complex last;      /**< Stores temporary values                  */
  complex act;       /**< Stores temporary values                  */
  complex *xp;       /**< Auxilliary pointer                       */
  complex *xm;       /**< Auxilliary pointer                       */
  int low, up;
  int lowe, upe;
  
  /* Determine lower and upper bounds for loops processing even and odd terms. */
  lowe = -plan->N + (plan->N%2);
  upe = -lowe;

  //for (n = -plan->N; n <= plan->N; n++)
  for (n = lowe; n <= upe; n += 2)
  {
    xm = &(plan->f_hat[NFSFT_INDEX(-1,n,plan)]);
    xp = &(plan->f_hat[NFSFT_INDEX(+1,n,plan)]);       
    for(k = 1; k <= plan->N; k++)
    {
      *xp *= 0.5;
      *xm-- = *xp++;
    }
  }
  
  /* Determine lower and upper bounds for loops processing even and odd terms. */
  low = -plan->N + (1-plan->N%2);
  up = -low;
  
  /* Process odd terms. */
  /* Incorporate sine term. */
  for (n = low; n <= up; n += 2)
  {
    plan->f_hat[NFSFT_INDEX(0,n,plan)] *= 2.0;
    xp = &(plan->f_hat[NFSFT_INDEX(-plan->N,n,plan)]);
    xm = &(plan->f_hat[NFSFT_INDEX(plan->N,n,plan)]);
    last = *xm;
    *xm = 0.5 * I * (0.5*xm[-1]);
    *xp++ = -(*xm--);
    for (k = plan->N-1; k > 0; k--)
    {
      act = *xm;
      *xm = 0.5 * I * (0.5*(xm[-1] - last));
      *xp++ = -(*xm--);
      last = act;
    }
    *xm = 0.0;
  }
}

inline void c2e_transposed(nfsft_plan *plan)
{
  int k;             /**< Degree \f$k\f$                           */
  int n;             /**< Order \f$k\f$                            */
  complex last;      /**< Stores temporary values                  */
  complex act;       /**< Stores temporary values                  */
  complex *xp;       /**< Auxilliary pointer                       */
  complex *xm;       /**< Auxilliary pointer                       */
  int low, up;
  int lowe, upe;
  
  /* Determine lower and upper bounds for loops processing even and odd terms. */
  lowe = -plan->N + (plan->N%2);
  upe = -lowe;

  //for (n = -plan->N; n <= plan->N; n++)
  for (n = lowe; n <= upe; n += 2)
  {
    xm = &(plan->f_hat[NFSFT_INDEX(-1,n,plan)]);
    xp = &(plan->f_hat[NFSFT_INDEX(+1,n,plan)]);       
    for(k = 1; k <= plan->N; k++)
    {
      *xp += *xm--;
      *xp++ *= 0.5;;
    }
  }
  
  /* Determine lower and upper bounds for loops processing even and odd terms. */
  low = -plan->N + (1-plan->N%2);
  up = -low;
  
  /* Process odd terms. */
  /* Incorporate sine term. */
  for (n = low; n <= up; n += 2)
  {
    xm = &(plan->f_hat[NFSFT_INDEX(-1,n,plan)]);
    xp = &(plan->f_hat[NFSFT_INDEX(+1,n,plan)]);       
    for(k = 1; k <= plan->N; k++)
    {
      *xp++ -= *xm--;
    }

    plan->f_hat[NFSFT_INDEX(0,n,plan)] = -0.25*I*plan->f_hat[NFSFT_INDEX(1,n,plan)];
    last = plan->f_hat[NFSFT_INDEX(1,n,plan)];
    plan->f_hat[NFSFT_INDEX(1,n,plan)] = -0.25*I*plan->f_hat[NFSFT_INDEX(2,n,plan)];
    xp = &(plan->f_hat[NFSFT_INDEX(2,n,plan)]);
    for (k = 2; k < plan->N; k++)
    {
      act = *xp;
      *xp++ = -0.25 * I * (xp[1] - last);
      last = act;
    }
    *xp = 0.25 * I * last;

    plan->f_hat[NFSFT_INDEX(0,n,plan)] *= 2.0;
  }
}

void nfsft_init(nfsft_plan *plan, int N, int M)
{
  /* Call nfsft_init_advanced with flags to allocate memory. */
  nfsft_init_advanced(plan, N, M, NFSFT_MALLOC_X | NFSFT_MALLOC_F | NFSFT_MALLOC_F_HAT);
} 
  
void nfsft_init_advanced(nfsft_plan* plan, int N, int M, 
                         unsigned int flags)
{
  /* Call nfsft_init_guro with the flags and default NFFT cut-off. */
  nfsft_init_guru(plan, N, M, flags, NFSFT_DEFAULT_NFFT_CUTOFF);
}

void nfsft_init_guru(nfsft_plan *plan, int N, int M, unsigned int flags, 
                     int nfft_cutoff)
{
  #ifdef NFSFT_OPTIMIZED
    /** NFFT size */
    int nfft_size;
    /** NFFT size */
    int nfct_size;
    /** FFTW size */
    int fftw_size;
    /** FFTW size */
    int fctw_size;
  #else
    /** NFFT size */
    int *nfft_size;    
    /** FFTW size */
    int *fftw_size;
  #endif  

  /* Save the flags in the plan. */
  plan->flags = flags;
  
  /* Save the bandwidth N and the number of samples M in the plan. */
  plan->N = N;
  plan->M_total = M;
  
  /* Calculate the next greater power of two with respect to the bandwidth N and
   * the corresponding exponent. */
  plan->NPT = next_power_of_2(plan->N);
  plan->t = (int)(log((double)plan->NPT)/log(2.0));
  
  #ifdef NFSFT_OPTIMIZED
    plan->maxMN = max(plan->M_total, plan->N+1);
  #endif
  
  /* Save length of array of Fourier coefficients. Owing to the data layout the 
   * length is (2N+1)(2NPT) */
  #ifdef NFSFT_OPTIMIZED
    plan->N_total = (2*plan->N+1)*(plan->maxMN);
  #else
    plan->N_total = (2*plan->N+2)*(2*plan->N+2);
  #endif
  
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
    plan->x = (double*) calloc(plan->M_total,2*sizeof(double));
  }
  
  /* Check if fast algorithm is activated. */
  if (plan->flags & NFSFT_NO_FAST_ALGORITHM)
  {
  }
  else
  {
    #ifdef NFSFT_OPTIMIZED
      /** TODO Check NFFT and FFTW sizes. */
      nfft_size = 2*plan->N+2;
      fftw_size = 4*plan->N;
      nfct_size = plan->N+2;
      fctw_size = 2*plan->N;
  
      /** \todo Check NFFT flags. */
      nfft_init_guru(&plan->plan_nfft, 1, &nfft_size, plan->M_total, &fftw_size, 
                         nfft_cutoff, PRE_PHI_HUT | PRE_PSI | FFT_OUT_OF_PLACE | 
                         FFTW_INIT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      nfct_init_guru(&plan->plan_nfct, 2, &nfct_size, plan->M_total, &fctw_size, 
                         nfft_cutoff, PRE_PHI_HUT | PRE_PSI | FFT_OUT_OF_PLACE | 
                         FFTW_INIT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      nfst_init_guru(&plan->plan_nfst, 2, &nfct_size, plan->M_total, &fctw_size, 
                         nfft_cutoff, PRE_PHI_HUT | PRE_PSI | FFT_OUT_OF_PLACE | 
                         FFTW_INIT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      /* Assign angle array. */
      plan->plan_nfft.x = plan->phi;
      /* Assign angle array. */
      plan->plan_nfct.x = plan->theta;
      /* Assign angle array. */
      plan->plan_nfst.x = plan->theta;
      /* Assign result array. */
      plan->plan_nfft.f = plan->f;
      /* Assign result array. Will be updated later anyway. */
      plan->plan_nfct.f = plan->f;
      /* Assign result array. Will be updated later anyway.  */
      plan->plan_nfst.f = plan->f;
      /* Precompute PSI, if necessary. */
    #else
      nfft_size = (int*) malloc(2*sizeof(int));
      fftw_size = (int*) malloc(2*sizeof(int));
    
      /** TODO Check NFFT and FFTW sizes. */
      nfft_size[0] = 2*plan->N+2;
      nfft_size[1] = 2*plan->N+2;
      fftw_size[0] = 4*plan->N;
      fftw_size[1] = 4*plan->N;
  
      /** \todo Check NFFT flags. */
      nfft_init_guru(&plan->plan_nfft, 2, nfft_size, plan->M_total, fftw_size, 
                         nfft_cutoff, FFT_OUT_OF_PLACE | 
                         FFTW_INIT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      /* Assign angle array. */
      plan->plan_nfft.x = plan->x;
      /* Assign result array. */
      plan->plan_nfft.f = plan->f;
      /* Assign Fourier coefficients array. */
      plan->plan_nfft.f_hat = plan->f_hat;
      /* Precompute. */     
      //if (plan->plan_nfft.nfft_flags & PRE_PSI) 
      //{  
        // TODO Add precomputation if neccessary and possible.
        //nfft_precompute_psi(&plan->plan_nfft); 
      //} 
      free(nfft_size);
      free(fftw_size);
    #endif
  }
}

void nfsft_precompute(int N, double kappa, 
                      unsigned int flags)
{ 
  int n;
  int k;
  
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
  wisdom.T_MAX = (int) (log((double)wisdom.N_MAX)/log(2.0));
  
  //fprintf(stdout,"\nnfsft_precompute: N = %d, N_MAX = %d, T_MAX = %d\n",N,wisdom.
  //  N_MAX,wisdom.T_MAX);

  /* Check, if precomputation for direct algorithms needs to be performed. */    
  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM)
  {
    wisdom.alpha = NULL;
    wisdom.beta = NULL;
    wisdom.gamma = NULL;
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
    //fflush(stdout);    
    alpha_al_all(wisdom.alpha,wisdom.N_MAX);
    beta_al_all(wisdom.beta,wisdom.N_MAX);
    gamma_al_all(wisdom.gamma,wisdom.N_MAX);
  }
  
  /* Check, if precomputation for direct algorithms needs to be performed. */ 
  /*fprintf(stdout,"wisdom.flags & NFSFT_NO_FAST_ALGORITHM = %d\n",
    wisdom.flags & NFSFT_NO_FAST_ALGORITHM);   
  fprintf(stdout,"wisdom.N_MAX = %d\n",wisdom.N_MAX);
  fprintf(stdout,"NFSFT_BREAKTHROUGH = %d\n",NFSFT_BREAKTHROUGH);*/
  if (wisdom.flags & NFSFT_NO_FAST_ALGORITHM)
  {
  }
  else if (wisdom.N_MAX >= NFSFT_BREAKTHROUGH)
  { 
    /* Precompute data for DPT. */
    
    /* Check, if recursion coefficients have already been calculated. */
    if (wisdom.alpha != NULL)
    {
      /* Use the recursion coefficients to precompute DPT data using persistent 
       * arrays. */
      wisdom.set = dpt_init(wisdom.N_MAX+1,wisdom.T_MAX,DPT_PERSISTENT_DATA);
      for (n = 0; n <= wisdom.N_MAX; n++)
      {
        //fprintf(stdout,"DPT precomputation: n = %d\n",n);
        /*for (k = 0; k < wisdom.N_MAX+2; k++)
        {
          /*fprintf(stdout,"alpha[%d] = %le, beta[%d] = %le, gamma[%d] = %le\n",
            k,
            wisdom.alpha[ROW(n)+k],
            k,
            wisdom.beta[ROW(n)+k],
            k,
            wisdom.gamma[ROW(n)+k]);
        }*/
        dpt_precompute(wisdom.set,n,&wisdom.alpha[ROW(n)],&wisdom.beta[ROW(n)],
          &wisdom.gamma[ROW(n)],n,kappa);
      } 
    }
    else
    {
      wisdom.alpha = (double*) malloc((wisdom.N_MAX+2)*sizeof(double));
      wisdom.beta = (double*) malloc((wisdom.N_MAX+2)*sizeof(double));
      wisdom.gamma = (double*) malloc((wisdom.N_MAX+2)*sizeof(double));
      wisdom.set = dpt_init(wisdom.N_MAX+1,wisdom.T_MAX,0U);
      for (n = 0; n <= wisdom.N_MAX; n++)
      {
        alpha_al_row(wisdom.alpha,wisdom.N_MAX,n);
        beta_al_row(wisdom.beta,wisdom.N_MAX,n);
        gamma_al_row(wisdom.gamma,wisdom.N_MAX,n);
        dpt_precompute(wisdom.set,n,wisdom.alpha,wisdom.beta,wisdom.gamma,n,kappa);
      } 
      free(wisdom.alpha);
      free(wisdom.beta);
      free(wisdom.gamma);
      wisdom.alpha = NULL;
      wisdom.beta = NULL;
      wisdom.gamma = NULL;
    }    
  }

  /* Wisdom has been initialised. */
  wisdom.initialized = true;
}
 
void nfsft_forget()
{
  if (wisdom.flags & NFSFT_NO_DIRECT_ALGORITHM)
  {
  }
  else
  {
    /* Free arrays holding three-term recurrence coefficients. */
    free(wisdom.alpha);   
    free(wisdom.beta);    
    free(wisdom.gamma);
    wisdom.alpha = NULL;
    wisdom.beta = NULL;
    wisdom.gamma = NULL;    
  }
  
  if (wisdom.flags & NFSFT_NO_FAST_ALGORITHM)
  {
  }
  else if (wisdom.N_MAX >= NFSFT_BREAKTHROUGH)
  {
    dpt_finalize(wisdom.set); 
  }    

  /* Wisdom is now uninitialised. */
  wisdom.initialized = false;
}
                     
void nfsft_finalize(nfsft_plan *plan)
{
  nfft_finalize(&plan->plan_nfft);
  
  /* Allocate memory for spherical Fourier coefficients, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_F_HAT)
  {
    free(plan->f_hat);
  }
  
  /* Allocate memory for samples, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_F)
  {
    free(plan->f);
  }
  
  /* Allocate memory for nodes, if neccesary. */
  if (plan->flags & NFSFT_MALLOC_X)
  {
    free(plan->x);
  }  
}

void ndsft_trafo(nfsft_plan *plan)
{
  /** Node index */
  int m;
  /** Legendre index k */
  int k;
  /** Legendre index n */
  int n;
  /** n_abs = |n| */
  int n_abs;
  int j;
  
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
  //complex *temp;
  complex it1;
  complex it2;
  complex temp;
  
  /** The final result for a single node. */
  complex f_m;
  
  double stheta;
  double sphi;
  
  if (plan->flags & NFSFT_NORMALIZED)
  {
    for (k = 0; k <= plan->N; k++)
    {
      for (n = -k; n <= k; n++)
      {
        plan->f_hat[NFSFT_INDEX(k,n,plan)] *= 
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
      plan->f[m] = plan->f_hat[NFSFT_INDEX(0,0,plan)];
    }
  }
  else
  {   
    /* Evaluate 
     * \[
     *   \sum_{k=0}^N \sum_{n=-k}^k a_k^n P_k^{|n|}(\cos\theta_m) e^{i n \phi_m} 
     *     = \sum_{n=-N}^N \sum_{k=|n|}^N a_k^n P_k^{|n|}(\cos\theta_m) 
     *     e^{i n \phi_m}. 
     * \] */
    for (m = 0; m < plan->M_total; m++)
    {
      stheta = cos(2.0*PI*plan->x[2*m+1]);
      sphi = 2.0*PI*plan->x[2*m];
      
      /* Initialize result for current node. */
      f_m = 0.0;

      /* For n = -N,...,N, evaluate 
       * \[
       *   b_n := \sum_{k=|n|}^N a_k^n P_k^{|n|}(\cos\vartheta_m)
       * using Clenshaw's algorithm. */
      for (n = -plan->N; n <= plan->N; n++)
      {
        /* Get Fourier coefficients vector. */
        a = &(plan->f_hat[NFSFT_INDEX(0,n,plan)]);
        
        /* Take absolute value of n. */
        n_abs = abs(n);
        
        /* Get three-term recurrence coefficients vectors. */
        alpha = &(wisdom.alpha[ROW(n_abs)]);
        gamma = &(wisdom.gamma[ROW(n_abs)]);
         
        /* Clenshaw's algorithm */        
        it2 = a[plan->N];
        it1 = a[plan->N-1];
        for (k = plan->N; k > n_abs + 1; k--)
        {
          temp = a[k-2] + it2 * gamma[k];
          it2 = it1 + it2 * alpha[k] * stheta;
          it1 = temp; 
        }
        
        /* Compute final step if neccesary. */
        if (n_abs < plan->N)
        {  
          //temp[0+n_abs] += temp[1+n_abs] * wisdom.alpha[ROWK(n_abs)+1] * stheta;
          it2 = it1 + it2 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
        }

        /* Write final result b_n of multiplication by normalization constant to 
         * array b  = (b_{-M},...,b_M). */
        f_m += it2 * wisdom.gamma[ROW(n_abs)] *
          pow(1- stheta * stheta, 0.5*n_abs) * cexp(I*n*sphi);
      }
            
      /* Write result to vector f. */
      plan->f[m] = f_m;                
    }
  }  
}

void ndsft_adjoint(nfsft_plan *plan)
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
  complex *a;
  /** Pointer to auxilliary array for Clenshaw algorithm. */
  complex it1;
  complex it2;
  complex temp;
  
  double stheta;
  double sphi;
  
  /*for (n = -plan->N; n <= plan->N; n++)
  {
    for (k = abs(n); k <= plan->N; k++)
    {
      plan->f_hat[NFSFT_INDEX(k,n,plan)] = 0.0;
    }  
  }*/
  memset(plan->f_hat,0U,plan->N_total*sizeof(complex));
  
  /* Distinguish by bandwidth N. */
  if (plan->N == 0)
  {
    /* Constant function */
    for (m = 0; m < plan->M_total; m++)
    {
      plan->f_hat[NFSFT_INDEX(0,0,plan)] += plan->f[m];
    }
  }
  else
  {   
    for (m = 0; m < plan->M_total; m++)
    {
      stheta = cos(2.0*PI*plan->x[2*m+1]);
      sphi = 2.0*PI*plan->x[2*m];
      
      for (n = -plan->N; n <= plan->N; n++)
      {
        /* Take absolute value of n. */
        n_abs = abs(n);

        /* Get three-term recurrence coefficients vectors. */
        alpha = &(wisdom.alpha[ROW(n_abs)]);
        gamma = &(wisdom.gamma[ROW(n_abs)]);
        
        it1 = plan->f[m] * wisdom.gamma[ROW(n_abs)] *
          pow(1 - stheta * stheta, 0.5*n_abs) * cexp(-I*n*sphi);
        plan->f_hat[NFSFT_INDEX(n_abs,n,plan)] += it1; 
        /*fprintf(stdout,"\nndsft_adjoint: cos(theta) = %le\n",stheta);
        fprintf(stdout,"\nndsft_adjoint: phi = %le\n",sphi);
        fprintf(stdout,"\nndsft_adjoint: gamma = %le\n",wisdom.gamma[ROW(n_abs)]);
        fprintf(stdout,"\nndsft_adjoint: f[%d] = %le + I*%le\n",m,plan->f[m]);
        fprintf(stdout,"\nndsft_adjoint: (1-%le^2)^0.5*%d = %le\n",stheta,n_abs,
          pow(1 - stheta * stheta, 0.5*n_abs));
        fprintf(stdout,"\nndsft_adjoint: e^-I*%d*%le = %le + I*%le\n",n,sphi,
          creal(cexp(-I*n*sphi)),
          cimag(cexp(-I*n*sphi)));*/
                
        /* Compute final step if neccesary. */
        if (n_abs < plan->N)
        {  
          it2 = it1 * wisdom.alpha[ROWK(n_abs)+1] * stheta;
          plan->f_hat[NFSFT_INDEX(n_abs+1,n,plan)] += it2; 
        } 

        /* Transposed Clenshaw algorithm */        
        for (k = n_abs+2; k <= plan->N; k++)
        {
          temp = it2;
          it2 = alpha[k] * stheta * it2 + gamma[k] * it1;
          it1 = temp;
          plan->f_hat[NFSFT_INDEX(k,n,plan)] += it2; 
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
        plan->f_hat[NFSFT_INDEX(k,n,plan)] *= 
          sqrt((2*k+1)/(4.0*PI));
      }
    }
  }
}

void nfsft_trafo(nfsft_plan *plan)
{
  /** Counter for loops */
  int k,n;
  
  if (wisdom.initialized == 0 || plan->N > wisdom.N_MAX)
  {
    return;
  }  
  
  if (plan->N < NFSFT_BREAKTHROUGH)
  {
    ndsft_trafo(plan);
  }
  else if (((wisdom.flags & NFSFT_BANDWIDTH_WINDOW) == 0U || 
    plan->N > (wisdom.N_MAX>>1)) && plan->N <= wisdom.N_MAX)
  {
    plan->plan_nfft.x = plan->x;
    plan->plan_nfft.f = plan->f;
    plan->plan_nfft.f_hat = plan->f_hat;
    
    /** Check for normalization. */
    if (plan->flags & NFSFT_NORMALIZED)
    {
      for (k = 0; k <= plan->N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          plan->f_hat[NFSFT_INDEX(k,n,plan)] *= 
            sqrt((2*k+1)/(4.0*PI));
        }
      }
    }
    
    if (plan->flags & NFSFT_USE_DPT)
    {
      /* Compute DPT. */
      for (n = -plan->N; n <= plan->N; n++) 
      {
        dpt_trafo(wisdom.set,abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
    }
    else
    {
      /* Compute FPT. */
      for (n = -plan->N; n <= plan->N; n++) 
      {
        fpt_trafo(wisdom.set,abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
    }

    /* Convert Chebyshev coefficients to Fourier coefficients. */
    c2e(plan); 
    
    /* Execute NFFT. */
    if (plan->flags & NFSFT_USE_NDFT)
    {
      ndft_trafo(&plan->plan_nfft);
    }
    else
    {
      nfft_trafo(&plan->plan_nfft);
    }    
  } 
}

void nfsft_adjoint(nfsft_plan *plan)
{
    /** Counter for loops */
  int k,n;
  
  if (wisdom.initialized == 0 || plan->N > wisdom.N_MAX)
  {
    return;
  }  
  
  if (plan->N < NFSFT_BREAKTHROUGH)
  {
    ndsft_adjoint(plan);
  }
  else if (((wisdom.flags & NFSFT_BANDWIDTH_WINDOW) == 0U || 
    plan->N > (wisdom.N_MAX>>1)) && plan->N <= wisdom.N_MAX)
  {
    plan->plan_nfft.x = plan->x;
    plan->plan_nfft.f = plan->f;
    plan->plan_nfft.f_hat = plan->f_hat;
    
    /* Execute NFFT. */
    if (plan->flags & NFSFT_USE_NDFT)
    {
      ndft_adjoint(&plan->plan_nfft);
    }
    else
    {
      nfft_adjoint(&plan->plan_nfft);
    }    

    /* Convert Chebyshev coefficients to Fourier coefficients. */
    c2e_transposed(plan); 
    
    if (plan->flags & NFSFT_USE_DPT)
    {
      /* Compute DPT. */
      for (n = -plan->N; n <= plan->N; n++) 
      {
        dpt_transposed(wisdom.set,abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
    }
    else
    {
      /* Compute FPT. */
      for (n = -plan->N; n <= plan->N; n++) 
      {
        fpt_transposed(wisdom.set,abs(n),
          &plan->f_hat[NFSFT_INDEX(abs(n),n,plan)],
          &plan->f_hat[NFSFT_INDEX(0,n,plan)],
          plan->N,0U);
      }
    }    

    /** Check for normalization. */
    if (plan->flags & NFSFT_NORMALIZED)
    {
      for (k = 0; k <= plan->N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          plan->f_hat[NFSFT_INDEX(k,n,plan)] *= 
            sqrt((2*k+1)/(4.0*PI));
        }
      }
    }
    
    if (plan->flags & NFSFT_ZERO_F_HAT)
    {
      for (n = -plan->N; n <= plan->N+1; n++)
      {
        memset(&plan->f_hat[NFSFT_INDEX(-plan->N-1,n,plan)],0U,(plan->N+1+abs(n))*
          sizeof(complex));
        /*for (k = -plan->N-1; k < abs(n); k++)
        {
          plan->f_hat[NFSFT_INDEX(k,n,plan)] = 0.0;
        }*/
      }
    }  
  } 
}
