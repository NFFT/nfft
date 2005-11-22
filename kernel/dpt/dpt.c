#include "dpt.h"

#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "util.h" 

#define K_START_TILDE(x,y) (MAX(MIN(x,y-2),0))
#define N_TILDE(y) (y-1)
#define K_END_TILDE(x,y) MIN(x,y-1)
#define FIRST_L(x,y) ((int)floor((x)/(double)y))
#define LAST_L(x,y) ((int)ceil(((x)+1)/(double)y)-1)

typedef struct dpt_step_
{
  bool stable;                            /**< Indicates if the values 
                                               contained represent a fast or 
                                               a slow stabilized step.       */
  double **a11,**a12,**a21,**a22;         /**< The matrix components         */
  double *gamma;                          /**<                               */
} dpt_step;

typedef struct dpt_data_
{
  dpt_step **steps;                       /**< The cascade summation steps   */
  int k_start;
  double *alphaN;
  double *betaN;
  double *gammaN;
  double alpha_0;
  double beta_0;
  double gamma_m1;
  /* Data for direct transform. */
  double *alpha;
  double *beta;
  double *gamma;
} dpt_data;

typedef struct dpt_set_s_
{
  int flags;                              /**< The flags                     */
  int M;                                  /**< The number of DPT transforms  */
  int N;                                  /**< The transform length. Must be 
                                               a power of two.               */
  int t;                                  /**< The exponent of N             */
  dpt_data *dpt;                          /**< The DPT transform data        */
  double **xcvecs;                        /**< Array of pointers to arrays 
                                               containing the Chebyshev 
                                               nodes                         */
  double *xc;                             /**< Array for Chebychev-nodes.    */ 
  complex *work;                          /**< */
  complex *result;                        /**< */
  complex *vec3;
  complex *vec4;
  complex *z;
  fftw_plan *plans_dct3;                  /**< Transform plans for the fftw 
                                               library                       */
  fftw_plan *plans_dct2;                  /**< Transform plans for the fftw 
                                               library                       */  
  fftw_r2r_kind *kinds;                   /**< Transform kinds for fftw 
                                               library                       */
  fftw_r2r_kind *kindsr;                  /**< Transform kinds for fftw 
                                               library                       */
  
  int *lengths; /**< Transform lengths for fftw library */  
  
  /* Data for slow transforms. */
  double *xc_slow;
} dpt_set_s;

void auvxpwy(double a, complex* u, complex* x, double* v, complex* y, 
  double* w, int n)
{
  int l;
  complex *u_ptr, *x_ptr, *y_ptr;
  double *v_ptr, *w_ptr;
  
  u_ptr = u;
  x_ptr = x;
  v_ptr = v;
  y_ptr = y;
  w_ptr = w;
  
  for (l = 0; l < n; l++)
  {
    *u++ = a * ((*v++) * (*x++) + (*w++) * (*y++));
  }
}

void dpt_do_step(complex  *a, complex *b, double *a11, double *a12, double *a21, 
                 double *a22, double gamma, int tau, dpt_set set)
{ 
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  double norm = 1.0/(length<<1);
  
  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;   
  
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)b,(double*)b);
  
  /* Check, if gamma is zero. */
  if (gamma == 0.0)
  {
    /* Perform multiplication only for second row. */
    auvxpwy(norm,b,b,a22,a,a21,length);
  }
  else 
  {
    /* Perform multiplication for both rows. */
    auvxpwy(norm,set->z,b,a22,a,a21,length);
    auvxpwy(norm*gamma,a,a,a11,b,a12,length);
    memcpy(b,set->z,length*sizeof(complex));    
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a);   
    /* Compensate for factors introduced by a raw DCT-II. */    
    a[0] *= 0.5;
  }  

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b);  
  /* Compensate for factors introduced by a raw DCT-II. */      
  b[0] *= 0.5;  
}

void eval_clenshaw(double *x, double *y, int size, int k, double *alpha, 
  double *beta, double *gamma)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector 
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  double a,b,x_val_act,a_old;
  double *x_act, *y_act;  
  double *alpha_act, *beta_act, *gamma_act;
  
  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;
    
    if (k == 0)
    {  
      *y_act = 1.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gamma[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));           
         b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);                  
    }
    x_act++;
    y_act++;
  }
}

int eval_clenshaw_thresh(double *x, double *y, int size, int k, double *alpha, 
  double *beta, double *gamma, double threshold)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector 
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  double a,b,x_val_act,a_old;
  double *x_act, *y_act;
  double *alpha_act, *beta_act, *gamma_act;
  
  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;
    
    if (k == 0)
    {  
     *y_act = 1.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gamma[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));           
         b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);                  
      if (fabs(*y_act) > threshold)
      {
        return 1;
      }
    }
    x_act++;
    y_act++;
  }
  return 0;
}

dpt_set dpt_init(const int M, const int t, const unsigned int flags)
{
  /** Polynomial length */
  int plength;
  /** Cascade level */
  int tau;
  /** Index m */
  int m;
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
  set->dpt = malloc((M+1)*sizeof(dpt_data));
  
  /* Initialize with NULL pointer. */
  for (m = 0; m <= set->M; m++)
  {
    set->dpt[m].steps = (dpt_step**)NULL;
  }
 
  /* Check if fast transform is activated. */
  if (set->flags & DPT_NO_FAST_TRANSFORM)
  {
  }
  else
  {
    /* Create arrays with Chebyshev nodes. */  
    
    /* Initialize array with Chebyshev coefficients for the polynomial x. This 
     * would be trivially an array containing a 1 as second entry with all other 
     * coefficients set to zero. In order to compensate for the multiplicative 
     * factor 2 introduced by the DCT-III, we set this coefficient to 0.5 here. */
    set->xc = (double *) calloc(1<<set->t,sizeof(double));
    set->xc[1] = 0.5;
    
    /* Allocate memory for array of pointers to node arrays. */
    set->xcvecs = (double**) malloc((set->t-1)*sizeof(double*));
    /* For each polynomial length starting with 4, compute the Chebyshev nodes 
     * using a DCT-III. */
    plength = 4;
    for (tau = 1; tau < t; tau++)
    {
      /* Allocate memory for current array. */
      set->xcvecs[tau-1] = (double*) malloc(plength*sizeof(double));
      /* Create plan for DCT-III. */
      plan = fftw_plan_r2r_1d(plength, set->xc, set->xcvecs[tau-1], FFTW_REDFT01, 
                              FFTW_PRESERVE_INPUT);
      /* Execute it. */
      fftw_execute(plan);
      /* Destroy the plan. */
      fftw_destroy_plan(plan);
      plan = NULL;
      /* Increase length to next power of two. */
      plength = plength << 1;
    }
    
    /** Allocate memory for auxilliary arrays. */
    set->work = (complex*) malloc(2*set->N*sizeof(complex));
    set->result = (complex*) malloc(2*set->N*sizeof(complex));
    set->vec3 = (complex*) malloc(set->N*sizeof(complex));
    set->vec4 = (complex*) malloc(set->N*sizeof(complex));
    set->z = (complex*) malloc(set->N*sizeof(complex));
    
    /** Initialize FFTW plans. */
    set->plans_dct3 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(set->t-1));
    set->plans_dct2 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(set->t-1)); 
    set->kinds      = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    set->kinds[0]   = FFTW_REDFT01;
    set->kinds[1]   = FFTW_REDFT01;
    set->kindsr     = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    set->kindsr[0]  = FFTW_REDFT10;
    set->kindsr[1]  = FFTW_REDFT10;
    set->lengths    = (int*) malloc((set->t-1)*sizeof(int));
    for (tau = 0, plength = 4; tau < set->t-1; tau++, plength<<=1)
    {
      set->lengths[tau] = plength;
      set->plans_dct3[tau] = 
        fftw_plan_many_r2r(1, &set->lengths[tau], 2, (double*)set->vec3, NULL, 
                           2, 1, (double*)set->vec4, NULL, 2, 1, set->kinds, 
                           0);
      set->plans_dct2[tau] = 
        fftw_plan_many_r2r(1, &set->lengths[tau], 2, (double*)set->vec3, NULL, 
                           2, 1, (double*)set->vec4, NULL, 2, 1,set->kindsr, 
                           0);
    }
  }  
  
  if (set->flags & DPT_NO_SLOW_TRANSFORM)
  {
  }
  else
  { 
    set->xc_slow = (double*) malloc((set->N+1)*sizeof(double));
  }         
 
  /* Return the newly created DPT set. */ 
  return set;
}

void dpt_precompute(dpt_set set, const int m, const double const* alpha, 
                    const double const* beta, const double const* gamma, int k_start,
                    const double threshold)
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
  int k_start_tilde;
  int N_tilde;
  
  dpt_data *data;
  
  /* Allocate memory for DPT transform data. */
  //set->dpt[m] = (dpt_data*) malloc(sizeof(dpt_data));

  /* Get pointer to DPT data. */
  data = &(set->dpt[m]);
  
  /* Check, if already precomputed. */
  if (data->steps != NULL)
  {
    return;
  }
  
  /* Save k_start. */
  data->k_start = k_start;
  
  /* Check if fast transform is activated. */
  if (set->flags & DPT_NO_FAST_TRANSFORM)
  {
  }
  else
  { 
    /* Save recursion coefficients. */
    data->alphaN = (double*) malloc((set->t-1)*sizeof(complex));
    data->betaN = (double*) malloc((set->t-1)*sizeof(complex));
    data->gammaN = (double*) malloc((set->t-1)*sizeof(complex));
    for (tau = 2; tau <= set->t; tau++)
    {
      data->alphaN[tau-2] = alpha[1<<tau]; 
      data->betaN[tau-2] = beta[1<<tau]; 
      data->gammaN[tau-2] = gamma[1<<tau];
    }   
    data->alpha_0 = alpha[1];
    data->beta_0 = beta[1];
    data->gamma_m1 = gamma[0];
      
    k_start_tilde = K_START_TILDE(data->k_start,set->N);
    N_tilde = N_TILDE(set->N);
      
    // printf("k_start = %d, k_start_tilde = %d, N = %d, N_tilde = %d\n",
    // k_start,k_start_tilde,set->N,N_tilde);
  
    /* Allocate memory for the cascade with t = log_2(N) many levels. */
    data->steps = (dpt_step**) malloc(sizeof(dpt_step*)*set->t);
      
    /* For tau = 1,...t compute the matrices U_{n,tau,l}. */
    plength = 4;
    for (tau = 1; tau < set->t; tau++)
    {     
      /* Compute auxilliary values. */
      degree = plength>>1;     
      /* Compute first l. */
      firstl = FIRST_L(k_start_tilde,plength);
      /* Compute last l. */
      lastl = LAST_L(N_tilde,plength);
      //printf("tau = %d: plength = %d, degree = %d, firstl = %d, lastl = %d\n",
      //  tau, plength, degree, firstl, lastl);
  
      /* Compute number of matrices for this level. */
      //nsteps = lastl - firstl + 1;
        
      /* Allocate memory for current level. This level will contain 2^{t-tau-1} 
       * many matrices. */
      data->steps[tau] = (dpt_step*) fftw_malloc(sizeof(dpt_step) 
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
          eval_clenshaw(set->xcvecs[tau-1], a11, plength, degree-2, calpha, cbeta, 
            cgamma);
          eval_clenshaw(set->xcvecs[tau-1], a12, plength, degree-1, calpha, cbeta, 
            cgamma);
          calpha--;
          cbeta--;
          cgamma--;
          eval_clenshaw(set->xcvecs[tau-1], a21, plength, degree-1, calpha, cbeta, 
            cgamma);
          eval_clenshaw(set->xcvecs[tau-1], a22, plength, degree, calpha, cbeta, 
            cgamma);
          needstab = 0;  
        }
        else
        {
          needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a11, plength, degree-2, 
            calpha, cbeta, cgamma, threshold);
          if (needstab == 0)
          {
            /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+2). */
            needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a12, plength, degree-1, 
              calpha, cbeta, cgamma, threshold);
            if (needstab == 0)
            { 
              calpha--;
              cbeta--;
              cgamma--;
              /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1). */
              needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a21, plength, 
                degree-1, calpha, cbeta, cgamma, threshold);
              if (needstab == 0)
              { 
                /* Evaluate P_{2^{tau}}^n(\cdot,2^{tau+1}l+1). */
                needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a22, plength, 
                  degree, calpha, cbeta, cgamma, threshold);
              }
            }
          }
        }        
        
        /* Check if stabilization needed. */
        if (needstab == 0)
        {  
          data->steps[tau][l].a11 = (double**) fftw_malloc(sizeof(double*)); 
          data->steps[tau][l].a12 = (double**) fftw_malloc(sizeof(double*)); 
          data->steps[tau][l].a21 = (double**) fftw_malloc(sizeof(double*)); 
          data->steps[tau][l].a22 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].gamma = (double*) fftw_malloc(sizeof(double)); 
          /* No stabilization needed. */
          data->steps[tau][l].a11[0] = a11;
          data->steps[tau][l].a12[0] = a12;
          data->steps[tau][l].a21[0] = a21;
          data->steps[tau][l].a22[0] = a22;
          data->steps[tau][l].gamma[0] = gamma[plength*l+1+1];
          data->steps[tau][l].stable = true;
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
            data->steps[tau][l].a11 = (double**) fftw_malloc(sizeof(double*)); 
            data->steps[tau][l].a12 = (double**) fftw_malloc(sizeof(double*)); 
            data->steps[tau][l].a21 = (double**) fftw_malloc(sizeof(double*)); 
            data->steps[tau][l].a22 = (double**) fftw_malloc(sizeof(double*)); 
            data->steps[tau][l].gamma = (double*) fftw_malloc(sizeof(double)); 
  
            plength_stab = 1<<set->t;
  
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
            eval_clenshaw(set->xcvecs[set->t-2], a11, plength_stab, degree_stab-2, 
              calpha, cbeta, cgamma);
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
            eval_clenshaw(set->xcvecs[set->t-2], a12, plength_stab, degree_stab-1, 
              calpha, cbeta, cgamma);
            calpha--;
            cbeta--;
            cgamma--;
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
            eval_clenshaw(set->xcvecs[set->t-2], a21, plength_stab, degree_stab-1, 
              calpha, cbeta, cgamma);
            /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
            eval_clenshaw(set->xcvecs[set->t-2], a22, plength_stab, degree_stab+0, 
              calpha, cbeta, cgamma);
            
            data->steps[tau][l].a11[0] = a11;
            data->steps[tau][l].a12[0] = a12;
            data->steps[tau][l].a21[0] = a21;
            data->steps[tau][l].a22[0] = a22;
            data->steps[tau][l].gamma[0] =  gamma[1+1];
            data->steps[tau][l].stable = false;
          }  
          else
          {
            data->steps[tau][l].a11 = (double**) fftw_malloc((set->t-tau)*
              sizeof(double*)); 
            data->steps[tau][l].a12 = (double**) fftw_malloc((set->t-tau)*
              sizeof(double*)); 
            data->steps[tau][l].a21 = (double**) fftw_malloc((set->t-tau)*
              sizeof(double*)); 
            data->steps[tau][l].a22 = (double**) fftw_malloc((set->t-tau)*
              sizeof(double*)); 
            data->steps[tau][l].gamma = (double*) fftw_malloc((set->t-tau)*
              sizeof(double)); 
  
            for (tau_stab = tau-1; tau_stab <= set->t-2; tau_stab++)
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
              eval_clenshaw(set->xcvecs[tau_stab], a11, plength_stab, degree_stab-2, 
                calpha, cbeta, cgamma);
              /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
              eval_clenshaw(set->xcvecs[tau_stab], a12, plength_stab, degree_stab-1, 
                calpha, cbeta, cgamma);
              calpha--;
              cbeta--;
              cgamma--;
              /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
              eval_clenshaw(set->xcvecs[tau_stab], a21, plength_stab, degree_stab-1, 
                calpha, cbeta, cgamma);
              /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
              eval_clenshaw(set->xcvecs[tau_stab], a22, plength_stab, degree_stab+0, 
                calpha, cbeta, cgamma);
              
              data->steps[tau][l].a11[tau_stab-tau+1] = a11;
              data->steps[tau][l].a12[tau_stab-tau+1] = a12;
              data->steps[tau][l].a21[tau_stab-tau+1] = a21;
              data->steps[tau][l].a22[tau_stab-tau+1] = a22;
              data->steps[tau][l].gamma[tau_stab-tau+1] = gamma[1+1];
            }
            data->steps[tau][l].stable = false;
          }
        }
      }
      /** Increase polynomial degree to next power of two. */
      plength = plength << 1;
    }  
  }  

  if (set->flags & DPT_NO_SLOW_TRANSFORM)
  {
  }
  else
  { 
    /* Check, if recurrence coefficients must be copied. */
    if (set->flags & DPT_PERSISTENT_DATA)
    {
      data->alpha = alpha;
      data->beta = beta;
      data->gamma = gamma;
    }
    else
    {
      data->alpha = (double*) malloc((set->N+1)*sizeof(double));
      data->beta = (double*) malloc((set->N+1)*sizeof(double));
      data->gamma = (double*) malloc((set->N+1)*sizeof(double));
      memcpy(data->alpha,alpha,(set->N+1)*sizeof(double));
      memcpy(data->beta,beta,(set->N+1)*sizeof(double));
      memcpy(data->gamma,gamma,(set->N+1)*sizeof(double));
    }
  }
}

void dpt_trafo(dpt_set set, const int m, const complex const* x, complex *y, 
  const int k_end, const unsigned int flags)
{
  int k;
  
  if (set->flags & DPT_NO_SLOW_TRANSFORM)
  {
    return;
  }
  
  /* Fill array with Chebyshev nodes. */
  for (k = 0; k <= k_end; k++)
  {
    set->xc_slow[k] = cos((PI*(k+0.5))/(k_end+1));
  }
  
}

void fpt_trafo(dpt_set set, const int m, const complex const* x, complex *y, 
  const int k_end, const unsigned int flags)
{ 
  /* Get transformation data. */
  dpt_data *data = &(set->dpt[m]);  
  /** */
  const int Nk = next_power_of_2(k_end);
  /** */
  const int tk = (int)(ceil(log((double)k_end)/log(2.0)));
  /** */
  int const k_start_tilde = K_START_TILDE(data->k_start,Nk);
  /** */
  int const k_end_tilde = K_END_TILDE(k_end,Nk);
   
  /** Level index \f$tau\f$ */
  int tau;
  /** Index of first block at current level */  
  int firstl;
  /** Index of last block at current level */
  int lastl;
  /** Block index \f$l\f$ */
  int l;
  /** Length of polynomial coefficient arrays at next level */
  int plength;
  /** Polynomial array length for stabilization */
  int plength_stab;
  /** Current matrix \f$U_{n,tau,l}\f$ */
  dpt_step *step;
  /** */
  fftw_plan plan;
  int length = k_end+1;
  fftw_r2r_kind kinds[2] = {FFTW_REDFT01,FFTW_REDFT01};
  
  /** */
  double mygamma;

  /** Loop counter */
  int k;
  
  complex *work_ptr;
  const complex *x_ptr;
  complex *y_ptr;
  
  /* Check if fast transform is activated. */
  if (set->flags & DPT_NO_FAST_TRANSFORM)
  { 
    return;
  }  
  
  /*fprintf(stdout,"fpt_trafo: k_start = %d\n",data->k_start);
  fprintf(stdout,"fpt_trafo: k_end = %d\n",k_end);
  fprintf(stdout,"fpt_trafo: Nk = %d\n",Nk);
  fprintf(stdout,"fpt_trafo: tk = %d\n",tk);
  fprintf(stdout,"fpt_trafo: k_start_tilde = %d\n",k_start_tilde);
  fprintf(stdout,"fpt_trafo: k_end_tilde = %d\n",k_end_tilde);*/

  if (flags & DPT_FUNCTION_VALUES)
  {
    plan = fftw_plan_many_r2r(1, &length, 2, (double*)set->work, NULL, 2, 1, 
      (double*)set->work, NULL, 2, 1, kinds, 0U);
  }
  
  /* Initialize working arrays. */
  memset(set->result,0U,2*Nk*sizeof(complex));

  /* The first step. */
  
  /* Set the first 2*data->k_start coefficients to zero. */
  memset(set->work,0U,2*data->k_start*sizeof(complex));

  work_ptr = &set->work[2*data->k_start];
  x_ptr = x;

  for (k = 0; k < k_end_tilde-data->k_start+1; k++) 
  {
    *work_ptr++ = *x_ptr++;
    *work_ptr++ = 0;
  }

  /* Set the last 2*(set->N-1-k_end_tilde) coefficients to zero. */
  memset(&set->work[2*(k_end_tilde+1)],0U,2*(Nk-1-k_end_tilde)*sizeof(complex));

  /* If k_end == Nk, use three-term recurrence to map last coefficient x_{Nk} to 
   * x_{Nk-1} and x_{Nk-2}. */
  if (k_end == Nk)
  { 
    set->work[2*(Nk-2)]   += data->gammaN[tk-2]*x[Nk-data->k_start];
    set->work[2*(Nk-1)]   += data->betaN[tk-2]*x[Nk-data->k_start];
    set->work[2*(Nk-1)+1]  = data->alphaN[tk-2]*x[Nk-data->k_start];
  }
  
  /* Compute the remaining steps. */
  plength = 4;
  for (tau = 1; tau < tk; tau++)
  {     
    /* Compute first l. */
    firstl = FIRST_L(k_start_tilde,plength);
    /* Compute last l. */
    lastl = LAST_L(k_end_tilde,plength);
    
    /* Compute the multiplication steps. */
    for (l = firstl; l <= lastl; l++)
    {      
      /* Copy vectors to multiply into working arrays zero-padded to twice the length. */
      memcpy(set->vec3,&(set->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(complex));
      memcpy(set->vec4,&(set->work[(plength/2)*(4*l+3)]),(plength/2)*sizeof(complex));     
      memset(&set->vec3[plength/2],0U,(plength/2)*sizeof(complex));
      memset(&set->vec4[plength/2],0U,(plength/2)*sizeof(complex));
      
      /* Copy coefficients into first half. */
      memcpy(&(set->work[(plength/2)*(4*l+2)]),&(set->work[(plength/2)*(4*l+1)]),(plength/2)*sizeof(complex));
      memset(&(set->work[(plength/2)*(4*l+1)]),0U,(plength/2)*sizeof(complex));
      memset(&(set->work[(plength/2)*(4*l+3)]),0U,(plength/2)*sizeof(complex));
      
      /* Get matrix U_{n,tau,l} */
      step = &(data->steps[tau][l]);
      
      /* Check if step is stable. */
      if (step->stable)
      {        
        /* Multiply third and fourth polynomial with matrix U. */       
        dpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0], 
          step->a21[0], step->a22[0], step->gamma[0], tau, set);

        if (step->gamma[0] != 0.0)
        {  
          for (k = 0; k < plength; k++)
          {
            set->work[plength*2*l+k] += set->vec3[k];
          }          
        }
        for (k = 0; k < plength; k++)
        {
          set->work[plength*(2*l+1)+k] += set->vec4[k];
        }  
      }
      else
      {        
        /* Stabilize. */

        /* The lengh of the polynomials */
        plength_stab = 1<<tk;
                
        /* Set rest of vectors explicitely to zero */
        memset(&set->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        memset(&set->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        
        /* Multiply third and fourth polynomial with matrix U. */
        if (set->flags & DPT_BANDWIDTH_WINDOW)
        {
          dpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0], 
            step->a21[0], step->a22[0], step->gamma[0], tk-1, set);
        }
        else
        {
          dpt_do_step(set->vec3, set->vec4, step->a11[tk-tau-1], 
            step->a12[set->t-tau-1], step->a21[tk-tau-1], 
            step->a22[set->t-tau-1], step->gamma[tk-tau-1], tk-1, set);
        }

        if (step->gamma[tk-tau-1] != 0.0)
        {  
          for (k = 0; k < plength_stab; k++)
          {
            set->result[k] += set->vec3[k];
          }          
        }
        for (k = 0; k < plength_stab; k++)
        {
          set->result[plength_stab+k] += set->vec4[k];
        }  
      }
    }
    /* Double length of polynomials. */
    plength = plength<<1;
  } 
  
  /* Add the resulting cascade coeffcients to the coeffcients accumulated from 
   * the stabilization steps. */
  for (k = 0; k < 2*Nk; k++)
  {
    set->result[k] += set->work[k];
  }  
 
  /* The last step. Compute the Chebyshev coeffcients c_k^n from the 
   * polynomials in front of P_0^n and P_1^n. */
  /*for (k = 0; k < 2*set->N; k++)
  {
    fprintf(stdout,"result[%d] = %1.16le +I*%1.16le\n",k,creal(set->result[k]),
      cimag(set->result[k]));
  }*/ 
    
  y[0] = data->gamma_m1*(set->result[0] + data->beta_0*set->result[Nk] + 
    data->alpha_0*set->result[Nk+1]*0.5);
  y[1] = data->gamma_m1*(set->result[1] + data->beta_0*set->result[Nk+1]+
    data->alpha_0*(set->result[Nk]+set->result[Nk+2]*0.5));
  y[k_end-1] = data->gamma_m1*(set->result[k_end-1] +
    data->beta_0*set->result[Nk+k_end-1] +
    data->alpha_0*set->result[Nk+k_end-2]*0.5);
  y[k_end] = data->gamma_m1*0.5*data->alpha_0*set->result[Nk+k_end-1];
  for (k = 2; k <= k_end-2; k++)
  {
    y[k] = data->gamma_m1*(set->result[k] + data->beta_0*set->result[Nk+k] +
      data->alpha_0*0.5*(set->result[Nk+k-1]+set->result[Nk+k+1]));
  } 

  if (flags & DPT_FUNCTION_VALUES)
  {
    y[0] *= 2.0;
    fftw_execute_r2r(plan,(double*)y,(double*)y);
    fftw_destroy_plan(plan);
    for (k = 0; k <= k_end; k++)
    {
      y[k] *= 0.5;
    }
  }  
}

void dpt_transposed(dpt_set set, const int m, const complex const* x, complex *y, 
  const int k_end, const unsigned int flags)
{
}

void fpt_transposed(dpt_set set, const int m, const complex const* x, complex *y, 
  const int k_end, const unsigned int flags)
{
}

void dpt_finalize(dpt_set set)
{
  int tau;
  int l;
  int m;
  dpt_data *data;
  int k_start_tilde;
  int N_tilde;  
  int tau_stab;
  int firstl, lastl;
  int plength;
  
  /* TODO Clean up DPT transform data structures. */
  for (m = 0; m <= set->M; m++)
  {
    /* Check if precomputed. */
    data = &set->dpt[m];
    if (data->steps != (dpt_step**)NULL)
    {
      free(data->alphaN);
      free(data->betaN);
      free(data->gammaN);
      
      /* Free precomputed data. */
      k_start_tilde = K_START_TILDE(data->k_start,set->N);
      N_tilde = N_TILDE(set->N);
      plength = 4;
      for (tau = 1; tau < set->t; tau++)
      {     
        /* Compute first l. */
        firstl = FIRST_L(k_start_tilde,plength);
        /* Compute last l. */
        lastl = LAST_L(N_tilde,plength);
        
        /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
        for (l = firstl; l <= lastl; l++)
        {                     
          if (set->flags & DPT_NO_STABILIZATION || 
            data->steps[tau][l].stable == true ||
            data->steps[tau][l].stable == true && set->flags & 
              DPT_BANDWIDTH_WINDOW)
          {
            /* Free components. */
            free(data->steps[tau][l].a11[0]);
            free(data->steps[tau][l].a12[0]);
            free(data->steps[tau][l].a21[0]);
            free(data->steps[tau][l].a22[0]);
            data->steps[tau][l].a11[0] = NULL;
            data->steps[tau][l].a12[0] = NULL;
            data->steps[tau][l].a21[0] = NULL;
            data->steps[tau][l].a22[0] = NULL;
          }
          else
          {  
            for (tau_stab = tau-1; tau_stab <= set->t-2; tau_stab++)
            {
              /* Free components. */
              free(data->steps[tau][l].a11[tau_stab-tau+1]);
              free(data->steps[tau][l].a12[tau_stab-tau+1]);
              free(data->steps[tau][l].a21[tau_stab-tau+1]);
              free(data->steps[tau][l].a22[tau_stab-tau+1]);
            }
          }
          /* Free components. */
          free(data->steps[tau][l].a11);
          free(data->steps[tau][l].a12);
          free(data->steps[tau][l].a21);
          free(data->steps[tau][l].a22);
          free(data->steps[tau][l].gamma);          
          data->steps[tau][l].a11 = NULL;
          data->steps[tau][l].a12 = NULL;
          data->steps[tau][l].a21 = NULL;
          data->steps[tau][l].a22 = NULL;
          data->steps[tau][l].gamma = NULL;        
        }
        /* Free pointers for current level. */
        free(data->steps[tau]);
        data->steps[tau] = NULL;
        /* Double length of polynomials. */
        plength = plength<<1;        
      } 
      /* Free steps. */
      free(data->steps);
      data->steps = NULL;           
    }
    
    if (set->flags & DPT_NO_SLOW_TRANSFORM)
    {
    }
    else
    {
      /* Check, if recurrence coefficients must be copied. */
      if (set->flags & DPT_PERSISTENT_DATA)
      {
      }
      else
      {
        free(data->alpha); 
        free(data->beta); 
        free(data->gamma); 
      }      
      data->alpha = NULL;
      data->beta = NULL;
      data->gamma = NULL;
    }
  }
  
  /* Delete array of DPT transform data. */
  free(set->dpt);
  set->dpt = NULL;
  
  /* Check if fast transform is activated. */
  if (set->flags & DPT_NO_FAST_TRANSFORM)
  {
  }
  else
  { 
    /* Delete arrays of Chebyshev nodes. */
    free(set->xc);
    set->xc = NULL;
    for (tau = 1; tau < set->t; tau++)
    {
      free(set->xcvecs[tau-1]);
      set->xcvecs[tau-1] = NULL;
    }
    free(set->xcvecs);
    set->xcvecs = NULL;
    
    /* Free auxilliary arrays. */
    free(set->work);  
    free(set->result);
    free(set->vec3);  
    free(set->vec4);  
    free(set->z);
    set->work = NULL;
    set->result = NULL;
    set->vec3 = NULL;
    set->vec4 = NULL;
    set->z = NULL;
    
    /* Free FFTW plans. */
    for(tau = 0; tau < set->t-1; tau++)
    {
      fftw_destroy_plan(set->plans_dct3[tau]);
      fftw_destroy_plan(set->plans_dct2[tau]);
      set->plans_dct3[tau] = NULL;
      set->plans_dct2[tau] = NULL;
    }  
     
    free(set->plans_dct3);
    free(set->plans_dct2);
    free(set->kinds);
    free(set->kindsr);
    free(set->lengths);
    set->plans_dct3 = NULL;
    set->plans_dct2 = NULL;
    set->kinds = NULL;
    set->kindsr = NULL;
    set->lengths = NULL;
  }
  
  if (set->flags & DPT_NO_SLOW_TRANSFORM)
  {
  }
  else
  { 
    free(set->xc_slow);
  }         
      
  /* Free DPT set structure. */
  free(set);
}
