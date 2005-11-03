#in                                clude "dpt.h"

#include <fftw3.h>

dpt_set dpt_init(const int M, const int t, const int flags)
{
  /** Polynomial length */
  int plength;
  /** Cascade level */
  int tau;
  /** FFTW plan */
  fftw_plan plan;  
  
  /* Allocate memory for new DPT set. */
  dpt_set_s *set = malloc(sizeof(dpt_set_s));
  
  /* Save parameters in structure. */
  set->flags = flags;
  set->M = M;
  set->t = t;
  set->N = 1<<t;
  
  /* Allocate memory for L transforms. */
  set->dpt = malloc(L*sizeof(dpt_data));
 
  /* Create arrays with Chebyshev nodes. */  
  
  /* Initialize array with Chebyshev coefficients for the polynomial x. This 
   * would be trivially an array containing a 1 as second entry with all other 
   * coefficients set to zero. In order to compensate for the multiplicative 
   * factor 2 introduced by the DCT-III, we set this coefficient to 0.5 here. */
  set->xc = (double *) calloc(1<<set->t,sizeof(double));
  set->xc[1] = 0.5;
  
  /* Allocate memory for array of pointers to node arrays. */
  set->xvecs = (double**) malloc((set->t-1)*sizeof(double*));
  /* For each polynomial length starting with 4, compute the Chebyshev nodes 
   * using a DCT-III. */
  plength = 4;
  for (tau = 1; tau < t; tau++)
  {
    /* Allocate memory for current array. */
    set->xvecs[tau-1] = (double*) malloc(plength*sizeof(double));
    /* Create plan for DCT-III. */
    plan = fftw_plan_r2r_1d(plength, set->xc, set->xvecs[tau-1], FFTW_REDFT01, 
                            FFTW_PRESERVE_INPUT);
    /* Execute it. */
    fftw_execute(plan);
    /* Destroy the plan. */
    fftw_destroy_plan(plan);
    /* Increase length to next power of two. */
    plength = plength << 1;
  }
 
  /* Return the newly created DPT set. */ 
  return set;
}

void dpt_precompute(dpt_set set, const int m, double const* alpha, 
                    double const* beta, double const* gamma, int k_start)
{
  
  int tau;          /**< Cascade level                                       */
  int l;            /**< Level index                                         */
  int plength;      /**< Length of polynomials for the next level in the 
                         cascade                                             */
  int degree;       /**< Degree of polynomials for the current level in the 
                         cascade                                             */
  int firstl;       /**< First index l for current cascade level             */ 
  int lastl;        /**< Last index l for current cascade level and current  */  
  int tau_stab;     /**< Cascade level for stabilization                     */
  int plength_stab; /**< Length of polynomials for the next level in the 
                         cascade for stabilization                           */
  int degree_stab;  /**< Degree of polynomials for the current level in the 
                         cascade for stabilization                           */
  double *a11;      /**< Array containing function values of the 
                         (1,1)-component of U_k^n.                           */
  double *a12;      /**< Array containing function values of the 
                         (1,2)-component of U_k^n.                           */
  double *a21;      /**< Array containing function values of the 
                         (2,1)-component of U_k^n.                           */
  double *a22;      /**< Array containing function values of the 
                         (2,2)-component of U_k^n.                           */
  double *calpha;
  double *cbeta;
  double *cgamma;
  int needstab = 0; /**< Used to indicate that stabilization is neccessary.  */
  
  dpt_data *data;
  
  /* Allocate memory for DPT transform data. */
  set->dpt[m] = malloc(sizeof(dpt_data));

  /* Get pointer to DPT data. */
  data = &(set->dpt[m]);
    
  int k_start_tilde = max(min(k_start,set->N-2));
  int N_tilde = set->N-1;
    
  // printf("k_start = %d, k_start_tilde = %d, N = %d, N_tilde = %d\n",
  // k_start,k_start_tilde,set->N,N_tilde);

  /* Allocate memory for the cascade with t = log_2(N) many levels. */
  data->steps = (dpt_step**) fftw_malloc(sizeof(struct dpt_step *) * t);
    
  /* For tau = 1,...t compute the matrices U_{n,tau,l}. */
  plength = 4;
  for (tau = 1; tau < t; tau++)
  {     
    /* Compute auxilliary values. */
    degree = plength>>1;     
    /* Compute first l. */
    firstl = FIRST_L;
    /* Compute last l. */
    lastl = LAST_L;
    //printf("tau = %d: plength = %d, degree = %d, firstl = %d, lastl = %d\n",
    //  tau, plength, degree, firstl, lastl);

    /* Compute number of matrices for this level. */
    //nsteps = lastl - firstl + 1;
      
    /* Allocate memory for current level. This level will contain 2^{t-tau-1} 
     * many matrices. */
    data->steps[tau] = (struct dpt_step*) fftw_malloc(sizeof(struct dpt_step) 
                       * (lastl+1)); 
    
    /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
    for (l = firstl; l <= lastl; l++)
    {        
      //fprintf(stderr,"n = %d [%d,%d], tau = %d [%d,%d], l = %d [%d,%d]\n",n,0,M,tau,1,t-1,l,firstl,lastl);

      /* Allocate memory for the components of U_{n,tau,l}. */
      a11 = (double*) fftw_malloc(sizeof(double)*plength);
      a12 = (double*) fftw_malloc(sizeof(double)*plength);
      a21 = (double*) fftw_malloc(sizeof(double)*plength);
      a22 = (double*) fftw_malloc(sizeof(double)*plength);               
      
      /* Evaluate the associated polynomials at the 2^{tau+1} Chebyshev 
       * nodes. */
      
      /* Get the pointers to the three-term recurrence coeffcients. */
      calpha = &(alpha[plength*l+1+1]);
      cbeta = &(beta[plength*l+1+1]);
      cgamma = &(gamma[plength*l+1+1]);
      
      if (set->flags & DPT_NO_STABILIZATION)
      {
        /* Evaluate P_{2^{tau}-2}^n(\cdot,2^{tau+1}l+2). */
        eval_clenshaw(xvecs[tau-1], a11, plength, degree-2, calpha, cbeta, 
          cgamma);
        eval_clenshaw(xvecs[tau-1], a12, plength, degree-1, calpha, cbeta, 
          cgamma);
        calpha--;
        cbeta--;
        cgamma--;
        eval_clenshaw(xvecs[tau-1], a21, plength, degree-1, calpha, cbeta, 
          cgamma);
        eval_clenshaw(xvecs[tau-1], a22, plength, degree, calpha, cbeta, 
          cgamma);
        needstab = 0;  
      }
      else
      {
        needstab = eval_clenshaw_thresh(xvecs[tau-1], a11, plength, degree-2, 
          calpha, cbeta, cgamma, threshold);
        if (needstab == 0)
        {
          /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+2). */
          needstab = eval_clenshaw_thresh(xvecs[tau-1], a12, plength, degree-1, 
            calpha, cbeta, cgamma, threshold);
          if (needstab == 0)
          { 
            calpha--;
            cbeta--;
            cgamma--;
            /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1). */
            needstab = eval_clenshaw_thresh(xvecs[tau-1], a21, plength, 
              degree-1, calpha, cbeta, cgamma, threshold);
            if (needstab == 0)
            { 
              /* Evaluate P_{2^{tau}}^n(\cdot,2^{tau+1}l+1). */
              needstab = eval_clenshaw_thresh(xvecs[tau-1], a22, plength, 
                degree, calpha, cbeta, cgamma, threshold);
            }
          }
        }
      }        
      
      /* Check if stabilization needed. */
      if (needstab == 0)
      {  
        data->steps[tau][l]->a11 = (double**) fftw_malloc(sizeof(double*)); 
        data->steps[tau][l]->a12 = (double**) fftw_malloc(sizeof(double*)); 
        data->steps[tau][l]->a21 = (double**) fftw_malloc(sizeof(double*)); 
        data->steps[tau][l]->a22 = (double**) fftw_malloc(sizeof(double*)); 
        /* No stabilization needed. */
        data->steps[tau][l]->a11[0] = a11;
        data->steps[tau][l]->a12[0] = a12;
        data->steps[tau][l]->a21[0] = a21;
        data->steps[tau][l]->a22[0] = a22;
        data->steps[tau][l]->stable = true;
      }          
      else 
      {    
        /* Stabilize. */
        degree_stab = degree*(2*l+1);          
        
        /* Old arrays are to small. */
        fftw_free(a11);
        fftw_free(a12);
        fftw_free(a21);
        fftw_free(a22);

        if (set->flags & DPT_BANDWIDTH_WINDOW)
        {
          data->steps[tau][l]->a11 = (double**) fftw_malloc(sizeof(double*)); 
          data->steps[tau][l]->a12 = (double**) fftw_malloc(sizeof(double*)); 
          data->steps[tau][l]->a21 = (double**) fftw_malloc(sizeof(double*)); 
          data->steps[tau][l]->a22 = (double**) fftw_malloc(sizeof(double*)); 

          plength_stab = 1<<t;

          /* Allocate memory for arrays. */
          a11 = (double*) fftw_malloc(sizeof(double)*plength_stab);
          a21 = (double*) fftw_malloc(sizeof(double)*plength_stab);
          a21 = (double*) fftw_malloc(sizeof(double)*plength_stab);
          a22 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            
          /* Get the pointers to the three-term recurrence coeffcients. */
          calpha = &(alpha[2]);
          cbeta = &(beta[2]);
          cgamma = &(gamma[2]);         
          /* Evaluate P_{2^{tau}(2l+1)-2}^n(\cdot,2). */
          eval_clenshaw(xvecs[t-2], a11, plength_stab, degree_stab-2, calpha, 
            cbeta, cgamma);
          /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
          eval_clenshaw(xvecs[t-2], a12, plength_stab, degree_stab-1, calpha, 
            cbeta, cgamma);
          calpha--;
          cbeta--;
          cgamma--;
          /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
          eval_clenshaw(xvecs[t-2], a21, plength_stab, degree_stab-1, calpha, 
            cbeta, cgamma);
          /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
          eval_clenshaw(xvecs[t-2], a22, plength_stab, degree_stab+0, calpha, 
            cbeta, cgamma);
          
          data->steps[tau][l]->a11[0] = a11;
          data->steps[tau][l]->a12[0] = a12;
          data->steps[tau][l]->a21[0] = a21;
          data->steps[tau][l]->a22[0] = a22;
          data->steps[tau][l]->stable = false;
        }  
        else
        {
          data->steps[tau][l]->a11 = (double**) fftw_malloc((t-tau)*
            sizeof(double*)); 
          data->steps[tau][l]->a12 = (double**) fftw_malloc((t-tau)*
            sizeof(double*)); 
          data->steps[tau][l]->a21 = (double**) fftw_malloc((t-tau)*
            sizeof(double*)); 
          data->steps[tau][l]->a22 = (double**) fftw_malloc((t-tau)*
            sizeof(double*)); 

          for (tau_stab = tau-1; tau_stab <= t-2; tau_stab++)
          {
            //tau_stab = t-2;
            plength_stab = 1<<(tau_stab+2);
            /* Allocate memory for arrays. */
            a11 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            a12 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            a21 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            a22 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            
            /* Get the pointers to the three-term recurrence coeffcients. */
            calpha = &(alpha[2]);
            cbeta = &(beta[2]);
            cgamma = &(gamma[2]);         
            /* Evaluate P_{2^{tau}(2l+1)-2}^n(\cdot,2). */
            eval_clenshaw(xvecs[tau_stab], a11, plength_stab, degree_stab-2, 
              calpha, cbeta, cgamma);
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
            eval_clnshaw(xvecs[tau_stab], a12, plength_stab, degree_stab-1, 
              calpha, cbeta, cgamma);
            calpha--;
            cbeta--;
            cgamma--;
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
            eval_clenshaw(xvecs[tau_stab], a21, plength_stab, degree_stab-1, 
              calpha, cbeta, cgamma);
            /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
            eval_clenshaw(xvecs[tau_stab], a22, plength_stab, degree_stab+0, 
              calpha, cbeta, cgamma);
            
            data->steps[tau][l]->a11[tau_stab-tau+1] = a11;
            data->steps[tau][l]->a12[tau_stab-tau+1] = a12;
            data->steps[tau][l]->a21[tau_stab-tau+1] = a21;
            data->steps[tau][l]->a22[tau_stab-tau+1] = a22;
            data->steps[tau][l]->stable = false;
          }
        }
      }
    }
    /** Increase polynomial degree to next power of two. */
    plength = plength << 1;
  }  
}

void dpt_finalize(dpt_set set)
{
  /* TODO Clean up DPT transform data structures. */
  
  /* Delete array of DPT transform data. */
  free(set->dpt);
  
  /* Delete arrays of Chebyshev nodes. */
  free(set->xc);
  for (tau = 1; tau < t; tau++)
  {
    free(set->xcvecs[tau-1]);
  }
  free(set->xvecs);  
  
  /* Free DPT set structure. */
  free(set);
}