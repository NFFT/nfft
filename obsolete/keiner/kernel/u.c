#include "legendre.h"

struct U_type*** precomputeU(int t, double threshold, double *walpha, 
                             double *wbeta, double *wgamma)
{
  /** Maximum bandwidth */
  int M = 1<<t;
  /** Legendre index n */
  int n;
  /** Cascade level */
  int tau;
  /** Level index */
  int l;
  
  /** Length of polynomials for the current level in the cascade */
  int plength;
  /** Degree of polynomials for the current level in the cascade */
  int degree;
  /** First index l for current cascade level and current n */
  int firstl;
  /** Last index l for current cascade level and current n */
  int lastl;
  /** Number of matrices U for current cascade level and current n .*/
  int nsteps;
  
  /** 
    * Length of polynomials for the current level in the cascade for 
    * stabilization.
    */
  /** Cascade level for stabilization */
  int tau_stab;
  /** 
   * Length of polynomials for the current level in the cascade for 
   * stabilization 
   */
  int plength_stab;
  /** 
   * Degree of polynomials for the current level in the cascade for 
   * stabilization 
   */
  int degree_stab;
  
  /** Three-dimensional array of matrices U_{n,tau,l} */
  struct U_type*** U;
  /** Array containing function values of the (1,1)-component of U_k^n. */
  double *m1;
  /** Array containing function values of the (1,2)-component of U_k^n. */
  double *m2;
  /** Array containing function values of the (2,1)-component of U_k^n. */
  double *m3;
  /** Array containing function values of the (2,2)-component of U_k^n. */
  double *m4;
  /**  
    * Array for three-term recurrence coefficients of associated Legendre
    * functions. 
    */
  double *alpha;
  /**  
    * Array for three-term recurrence coefficients of associated Legendre
    * functions. 
    */
  double *beta;
  /**  
    * Array for three-term recurrence coefficients of associated Legendre
    * functions. 
    */
  double *gamma;
  /** Array of pointers to arrays containing the Chebyshev nodes */
  double **xvecs;
  /** Array for Chebychev-nodes. */
  double *xc;

  /** loop counter */
  int i;
  /** Used to indicate that stabilization is neccessary. */
  bool needstab = false;  
  
  
  /* Initialize array with Chebyshev coefficients for the polynomial x. This 
   * would be trivially an array containing a 1 as second entry with all other 
   * coefficients set to zero. In order to compensate for the multiplicative 
   * factor 2 introduced by the DCT-III, we set this coefficient to 0.5 here. */
  xc = (double *) calloc(1<<t,sizeof(double));
  xc[1] = 0.5;
  
  /* Allocate memory for array of pointers to node arrays. */
  xvecs = (double**) malloc((t-1)*sizeof(double*));
  /* For each polynomial length starting with 4, compute the Chebyshev nodes 
   * using a DCT-III. */
  plength = 4;
  for (tau = 1; tau < t; tau++)
  {
    /* Allocate memory for current array. */
    xvecs[tau-1] = (double*) malloc(plength*sizeof(double));
    /* Create plan for DCT-III. */
    fftw_plan plan = fftw_plan_r2r_1d(plength, xc, xvecs[tau-1], FFTW_REDFT01, 
                                      FFTW_PRESERVE_INPUT);
    /* Execute it. */
    fftw_execute(plan);
    /* Destroy the plan. */
    fftw_destroy_plan(plan);
    /* Increase length to next power of two. */
    plength = plength << 1;
  }
  
  
  /* Allocate memory for matrices U_k^n(\cdot,4l+1). */
  U = (struct U_type***) fftw_malloc(sizeof(struct U_type **) * (M+1));
  
  /* We have to precompute the matrices 
   *   U_{n,tau,l} := U_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1)
   * for n = 0,...,M; tau = 1,...,t and l = 0,...,2^{t-tau-1}-1. */ 
  
  /* For nleg = 0,...,M compute the matrices U_{n,tau,l}. */
  for (n = 0; n <= M; n++)
  {   
    /* Allocate memory for current matrix array. The cascade will have 
     * t = log_2(M) many levels. */
    U[n] = (struct U_type**) fftw_malloc(sizeof(struct U_type*) * t);
    
    /* For tau = 1,...t compute the matrices U_{n,tau,l}. */
    plength = 4;
    for (tau = 1; tau < t; tau++)
    {     
      /* Compute auxilliary values. */
	     degree = plength>>1;
      
      /* Compute first l. */
      firstl = 0; //pow2((int)log2(n)-(n==M?1:0)-tau-1);
      /* Compute last l. */
      lastl = (int)(((double)pow2(t))/plength) - 1;
      /* Compute number of matrices for this level. */
      nsteps = lastl - firstl + 1;
      
      /* Allocate memory for current matrix array. The level will contain
        * 2^{t-tau-1} many matrices. */
	     U[n][tau] = (struct U_type*) fftw_malloc(sizeof(struct U_type) * nsteps); 
      
      /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
      l = firstl;
  	   for (i = 0; i < nsteps; i++)
	     {        
        /* Allocate memory for the components of U_{n,tau,l}. */
	       m1 = (double*) fftw_malloc(sizeof(double)*plength);
	       m2 = (double*) fftw_malloc(sizeof(double)*plength);
	       m3 = (double*) fftw_malloc(sizeof(double)*plength);
	       m4 = (double*) fftw_malloc(sizeof(double)*plength);               
        
        /* Evaluate the associated Legendre polynomials at the 2^{tau+1} 
         * Chebyshev nodes. */
        
        /* Get the pointers to the three-term recurrence coeffcients. */
        alpha = &(walpha[ROW(n)+plength*l+1+1]);
        beta = &(wbeta[ROW(n)+plength*l+1+1]);
        gamma = &(wgamma[ROW(n)+plength*l+1+1]);
        /* Evaluate P_{2^{tau}-2}^n(\cdot,2^{tau+1}l+2). */
        needstab = eval_al_thresh(xvecs[tau-1], m1, plength, degree-2, alpha, beta, 
                                  gamma, threshold);
        if (needstab == false)
        {
          /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+2). */
          needstab = eval_al_thresh(xvecs[tau-1], m2, plength, degree-1, alpha, beta, 
                                    gamma, threshold);
          if (needstab == false)
          { 
            alpha--;
            beta--;
            gamma--;
            /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1). */
            needstab = eval_al_thresh(xvecs[tau-1], m3, plength, degree-1, alpha, 
                                      beta, gamma, threshold);
            if (needstab == false)
            { 
              /* Evaluate P_{2^{tau}}^n(\cdot,2^{tau+1}l+1). */
              needstab = eval_al_thresh(xvecs[tau-1], m4, plength, degree, alpha, 
                                        beta, gamma, threshold);
            }
          }
        }      
        
        /* Check if stabilization needed. */
        if (needstab == false)
        {  
          /* No stabilization needed. */
          U[n][tau][i].m1 = m1;
          U[n][tau][i].m2 = m2;
          U[n][tau][i].m3 = m3;
          U[n][tau][i].m4 = m4;
          U[n][tau][i].stable = true;
        }          
        else 
        {    
          //fprintf(stderr,"(%d,%d,%d)\n",n,tau,l);
          
          /* Stabilize. */
  		      degree_stab = degree*(2*l+1);          
          plength_stab = pow2(t);
          tau_stab = t-2;
          
          /* Old arrays are to small. */
          fftw_free(m1);
          fftw_free(m2);
          fftw_free(m3);
          fftw_free(m4);
          
          /* Allocate memory for arrays. */
          m1 = (double*) fftw_malloc(sizeof(double)*plength_stab);
          m2 = (double*) fftw_malloc(sizeof(double)*plength_stab);
          m3 = (double*) fftw_malloc(sizeof(double)*plength_stab);
          m4 = (double*) fftw_malloc(sizeof(double)*plength_stab);
          
          /* Get the pointers to the three-term recurrence coeffcients. */
          alpha = &(walpha[ROW(n)+2]);
          beta = &(wbeta[ROW(n)+2]);
          gamma = &(wgamma[ROW(n)+2]);         
          /* Evaluate P_{2^{tau}(2l+1)-2}^n(\cdot,2). */
          eval_al(xvecs[tau_stab], m1, plength_stab, degree_stab-2, alpha, 
                  beta, gamma);
          /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
          eval_al(xvecs[tau_stab], m2, plength_stab, degree_stab-1, alpha, 
                  beta, gamma); 
          alpha--;
          beta--;
          gamma--;
          /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
          eval_al(xvecs[tau_stab], m3, plength_stab, degree_stab-1, alpha, 
                  beta, gamma);
          /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
          eval_al(xvecs[tau_stab], m4, plength_stab, degree_stab+0, alpha, 
                  beta, gamma);
                    
          U[n][tau][i].m1 = m1;
          U[n][tau][i].m2 = m2;
          U[n][tau][i].m3 = m3;
          U[n][tau][i].m4 = m4;          
          U[n][tau][i].stable = false;                    
        }
        l++;
      }
      /** Increase polynomial degree to next power of two. */
      plength = plength << 1;
   	}
  }
  
  /* Free memory for Chebyshev nodes. */
  for (tau = 1; tau < t; tau++)
  {
    free(xvecs[tau-1]);
  }
  free(xvecs);
  free(xc);
    
  return U;
}


void multiplyU(complex  *a, complex *b, struct U_type u, int tau, int n, int l, 
               struct nfsft_transform_wisdom *tw, double gamma)
{ 
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  int length2 = length<<1;
  /** Auxilliary variable */
  complex z;
  /** Loop counter */
  int i;
  
  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;   
  
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)b,(double*)b);
  
  /* Check, if gamma_k^n is zero. This is the case when l <= n holds. */
  if (false/*l <= n*/)
  {
    /* Perform multiplication only for second row. */
    for (i = 0; i < length; i++)
    {
      b[i] *= u.m4[i];
      b[i] += u.m3[i] * a[i];
    }    
  }
  else 
  {
    /* Perform multiplication for both rows. */
    for (i = 0; i < length; i++)
    {
      z = u.m3[i] * a[i] + u.m4[i] * b[i];
      a[i] *= u.m1[i];
      a[i] += u.m2[i]*b[i];
      a[i] *= gamma;
      b[i]  = z;
    }
    
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)a,(double*)a);   
    /* Compensate for lack of normalization of DCT-II. */        
    for (i = 0; i < length; i++)
    {
      a[i] = a[i]/length2;
    }  
    /* Compensate for factors introduced by a raw DCT-II. */    
    a[0] *= 0.5;
  }
  
  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)b,(double*)b);  
  /* Compensate for lack of normalization of DCT-II. */        
  for (i = 0; i < length; i++)
  {
    b[i] = b[i]/length2;
  }      
  /* Compensate for factors introduced by a raw DCT-II. */      
  b[0] *= 0.5;  
}


void multiplyU_adjoint(complex  *a, complex *b, 
                       struct U_type u, int tau, int n, int l, 
                       struct nfsft_transform_wisdom *tw, double gamma)
{ 
  /** Used to store temporary valiues. */
  complex z;
  /** Counter for loops. */
  int i;
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  int length2 = length<<1;
  
  
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)b,(double*)b);

  /* Make copies. */
  memcpy(tw->a2,a,length*sizeof(complex));
  memcpy(tw->b2,b,length*sizeof(complex));
  
  /* Perform matrix multiplication. */
  for (i = 0; i < length; i++)
  {
    a[i] = a[i] * gamma * u.m1[i] + b[i] * u.m3[i];
    b[i] = tw->a2[i] * gamma * u.m2[i] + tw->b2[i] * u.m4[i];
  }    
  
  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)a,(double*)a);   
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)b,(double*)b);   

  /* Compensate for lack of normalization of DCT-II. */          
  for (i = 0; i < length; i++)
  {
    a[i] = a[i]/length2;
    b[i] = b[i]/length2;
  }    
}
