#include "dpt.h"

#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "util.h" 

#define FIRST_L (int)floor(k_start_tilde/(double)plength)
#define LAST_L (int)ceil((N_tilde+1)/(double)plength)-1

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
  set->dpt = malloc(M*sizeof(dpt_data));
 
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
 
  /* Return the newly created DPT set. */ 
  return set;
}

void dpt_precompute(dpt_set set, const int m, double const* alpha, 
                    double const* beta, double const* gamma, int k_start,
                    double threshold)
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
  //set->dpt[m] = (dpt_data*) malloc(sizeof(dpt_data));

  /* Get pointer to DPT data. */
  data = &(set->dpt[m]);
  
  /* Save k_start. */
  data->k_start = k_start;
  
  /* Save recursion coefficients. */
  data->alphaN = alpha[set->N]; 
  data->betaN = beta[set->N]; 
  data->gammaN = gamma[set->N]; 
  data->alpha_0 = alpha[1];
  data->beta_0 = beta[1];
  data->gamma_m1 = gamma[0];
    
  int k_start_tilde = MAX(MIN(data->k_start,set->N-2),0);
  int N_tilde = set->N-1;
    
  // printf("k_start = %d, k_start_tilde = %d, N = %d, N_tilde = %d\n",
  // k_start,k_start_tilde,set->N,N_tilde);

  /* Allocate memory for the cascade with t = log_2(N) many levels. */
  data->steps = (dpt_step**) fftw_malloc(sizeof(struct dpt_step *) * set->t);
    
  /* For tau = 1,...t compute the matrices U_{n,tau,l}. */
  plength = 4;
  for (tau = 1; tau < set->t; tau++)
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
            data->steps[tau][l].stable = false;
          }
        }
      }
    }
    /** Increase polynomial degree to next power of two. */
    plength = plength << 1;
  }  
}

void dpt_trafo(dpt_set set, const int m, const int k_end, complex *x)
{ 
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
  double mygamma;

  /** Loop counter */
  int j;
  
  complex *work_ptr;
  complex *x_ptr;
  
  /* Get transformation data. */
  dpt_data *data = &(set->dpt[m]);
  
  int k_start_tilde = MAX(MIN(data->k_start,set->N-2),0);
  int N_tilde = set->N-1;
  
  /* Initialize working arrays. */
  memset(set->work,0U,2*set->N*sizeof(complex));
  memset(set->result,0U,2*set->N*sizeof(complex));

  /* Set first n Fourier coefficients explicitly to zero. */
  memset(x,0U,data->k_start*sizeof(complex));
  memset(&(x[k_end+1]),0U,(set->N-k_end)*sizeof(complex));
  
  work_ptr = set->work;
  x_ptr = x;
  /* First step */ 
  for (j = 0; j < set->N; j++) 
  {
    *work_ptr = *x_ptr++;
    work_ptr += 2;
  }
  
  /* Use three-term recurrence to map last coefficient a_N to a_{N-1} and 
   * a_{N-2}. */
  set->work[2*(set->N-2)]   += data->gammaN*x[set->N];
  set->work[2*(set->N-1)]   += data->betaN*x[set->N];
  set->work[2*(set->N-1)+1]  = data->alphaN*x[set->N];
  
  /* Compute the remaining steps. */
  plength = 4;
  for (tau = 1; tau < set->t; tau++)
  {    
    /* Compute first l. */
    firstl = FIRST_L;
    /* Compute last l. */
    lastl = LAST_L;
    
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
          for (j = 0; j < plength; j++)
          {
            set->work[plength*2*l+j] += set->vec3[j];
          }          
        }
        for (j = 0; j < plength; j++)
        {
          set->work[plength*(2*l+1)+j] += set->vec4[j];
        }  
      }
      else
      {        
        /* Stabilize. */

        /* The lengh of the polynomials */
        plength_stab = 1<<set->t;
                
        /* Set rest of vectors explicitely to zero */
        memset(&set->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        memset(&set->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(complex));
        
        /* Multiply third and fourth polynomial with matrix U. */
        if (set->flags & DPT_BANDWIDTH_WINDOW)
        {
          dpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0], 
            step->a21[0], step->a22[0], step->gamma[0], set->t-1, set);
        }
        else
        {
          dpt_do_step(set->vec3, set->vec4, step->a11[set->t-tau-1], 
            step->a12[set->t-tau-1], step->a21[set->t-tau-1], 
            step->a22[set->t-tau-1], step->gamma[set->t-tau-1], set->t-1, set);
        }

        if (step->gamma[set->t-tau-1] != 0.0)
        {  
          for (j = 0; j < plength_stab; j++)
          {
            set->work[j] += set->vec3[j];
          }          
        }
        for (j = 0; j < plength_stab; j++)
        {
          set->work[plength_stab+j] += set->vec4[j];
        }  
      }
    }
    /* Double length of polynomials. */
    plength = plength<<1;
  } 
  
  /* Add the resulting cascade coeffcients to the coeffcients accumulated from 
   * the stabilization steps. */
  for (j = 0; j < 2*set->N; j++)
  {
    set->result[j] += set->work[j];
  }  
 
  /* The last step. Compute the Chebyshev coeffcients c_k^n from the 
   * polynomials in front of P_0^n and P_1^n. */ 
  x[0] = data->gamma_m1*(set->result[0]+data->beta_0*set->result[set->N]+
    data->alpha_0*set->result[set->N+1]*0.5);
  x[1] = data->gamma_m1*(set->result[1]+set->result[set->N+1]+
    data->alpha_0*(set->result[set->N]+set->result[set->N+2]*0.5));
  x[set->N-1] = data->gamma_m1*(set->result[set->N-1]+
    data->beta_0*set->result[set->N+set->N-1]+
    data->alpha_0*set->result[set->N+set->N-2]*0.5);
    /* !!!!!!!!!!!!!!!!!!! */
  x[set->N] = data->gamma_m1*(set->result[set->N+set->N]-
    set->result[set->N+set->N-1]*0.5);
  for (j = 2; j < set->N-1; j++)
  {
    x[j] = data->gamma_m1*(set->result[j]+data->beta_0*set->result[j+set->N]+
      data->alpha_0*(set->result[j+set->N-1]+set->result[j+set->N+1])*0.5);
  } 
}

void dpt_finalize(dpt_set set)
{
  int tau;
  /* TODO Clean up DPT transform data structures. */
  
  /* Delete array of DPT transform data. */
  free(set->dpt);
  
  /* Delete arrays of Chebyshev nodes. */
  free(set->xc);
  for (tau = 1; tau < set->t; tau++)
  {
    free(set->xcvecs[tau-1]);
  }
  free(set->xcvecs);
  
  /* Free auxilliary arrays. */
  free(set->work);  
  free(set->result);
  free(set->vec3);  
  free(set->vec4);  
  free(set->z);
  
  /* Free FFTW plans. */
  for(tau = 0; tau < set->t-1; tau++)
  {
    fftw_destroy_plan(set->plans_dct3[tau]);
    fftw_destroy_plan(set->plans_dct2[tau]);
  }  
   
  free(set->plans_dct3);
  free(set->plans_dct2);
  free(set->kinds);
  free(set->kindsr);
  free(set->lengths);
    
  /* Free DPT set structure. */
  free(set);
}
