#include "legendre.h"
#include "api.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>

struct U_type**** precomputeU(int t, double threshold, double *walpha, 
                             double *wbeta, double *wgamma, bool window)
{
  const int N = 1<<t;
  /** Maximum bandwidth */
  int M = 1<<t;
  /** Legendre index n */
  int n;
  /** Cascade level */
  int tau;
  /** Level index */
  int l;
  
  /** Length of polynomials for the next level in the cascade */
  int plength;
  /** Degree of polynomials for the current level in the cascade */
  int degree;
  /** First index l for current cascade level and current n */
  int firstl;
  /** Last index l for current cascade level and current n */
  int lastl;
  /** Number of matrices U for current cascade level and current n .*/
  //int nsteps;
  
  /** Cascade level for stabilization */
  int tau_stab;
  /** 
   * Length of polynomials for the next level in the cascade for 
   * stabilization 
   */
  int plength_stab;
  /** 
   * Degree of polynomials for the current level in the cascade for 
   * stabilization 
   */
  int degree_stab;
  
  /** Four-dimensional array of matrices U_{n,tau,l} */
  struct U_type ****U;
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
  //int i;
  /** Used to indicate that stabilization is neccessary. */
  bool needstab = false; 
  
  int nstab = 0;
  
  int ntilde;
  int Mtilde = min(M,N-1);
	
#ifdef LOGFILE
  int j;
	char filename[100];
  FILE *logfile = fopen(LOGFILENAME,"w");
  if (logfile != NULL)
  {
    fclose(logfile);
  }  
#endif	
  
  fprintf(stderr,"Threshold = %lf\n",threshold);  
    
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
  U = (struct U_type****) fftw_malloc(sizeof(struct U_type ***) * (M+1));
  
  /* We have to precompute the matrices 
   *   U_{n,tau,l} := U_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1)
   * for n = 0,...,M; tau = 1,...,t and l = 0,...,2^{t-tau-1}-1. */ 
  
  /* For nleg = 0,...,M compute the matrices U_{n,tau,l}. */
  for (n = 0; n <= M; n++)
  {   
    ntilde = max(min(n,/*N*/(1<<ngpt(n))-2),0);
    //printf("ntilde: n = %d, 1<<(ngpt(n))-2 = %d, min(n,1<<(ngpt(n))-2) = %d\n",n,(1<<ngpt(n))-2,ntilde);
    //fprintf(stderr,"n = %5d\n",n);
  		//fflush(stderr);
    /* Allocate memory for current matrix array. The cascade will have 
     * t = log_2(M) many levels. */
    U[n] = (struct U_type***) fftw_malloc(sizeof(struct U_type**) * t);
    
    /* For tau = 1,...t compute the matrices U_{n,tau,l}. */
    plength = 4;
    for (tau = 1; tau < t; tau++)
    {     
      /* Compute auxilliary values. */
	     degree = plength>>1;
      
      /* Compute first l. */
		  //printf("firstl: ntile = %d, plength = %d, n = %d, dings = %d\n",ntilde,plength,n,1<<ngpt(n));
      firstl = FIRST_L;
      /* Compute last l. */
      lastl = LAST_L;
      /* Compute number of matrices for this level. */
      //nsteps = lastl - firstl + 1;
      
      /* Allocate memory for current matrix array. The level will contain
        * 2^{t-tau-1} many matrices. */
		  //printf("Alloc tau: n = %d, tau = %d, firstl = %d, lastl = %d, lastl+1 = %d\n",n,tau,firstl,lastl,lastl+1);
			//fflush(stdout);
	     U[n][tau] = (struct U_type**) fftw_malloc(sizeof(struct U_type*) * (lastl+1)); 
      
      /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
  	   for (l = firstl; l <= lastl; l++)
	     {        

        //fprintf(stderr,"n = %d [%d,%d], tau = %d [%d,%d], l = %d [%d,%d]\n",n,0,M,tau,1,t-1,l,firstl,lastl);

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
        /*if (1 == 0)
        {
          eval_al(xvecs[tau-1], m1, plength, degree-2, alpha, beta, gamma);
          eval_al(xvecs[tau-1], m2, plength, degree-1, alpha, beta, gamma);
          alpha--;
          beta--;
          gamma--;
          eval_al(xvecs[tau-1], m3, plength, degree-1, alpha, beta, gamma);
          eval_al(xvecs[tau-1], m4, plength, degree, alpha, beta, gamma);
        }*/
        
        needstab = eval_al_thresh(xvecs[tau-1], m1, plength, degree-2, alpha, 
                                  beta, gamma, threshold);
        if (needstab == false)
        {
          /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+2). */
          needstab = eval_al_thresh(xvecs[tau-1], m2, plength, degree-1, alpha, 
                                    beta, gamma, threshold);
          if (needstab == false)
          { 
            alpha--;
            beta--;
            gamma--;
            /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1). */
            needstab = eval_al_thresh(xvecs[tau-1], m3, plength, degree-1, 
                                      alpha, beta, gamma, threshold);
            if (needstab == false)
            { 
              /* Evaluate P_{2^{tau}}^n(\cdot,2^{tau+1}l+1). */
              needstab = eval_al_thresh(xvecs[tau-1], m4, plength, degree, 
                                        alpha, beta, gamma, threshold);
            }
          }
        }      
        
        /* Check if stabilization needed. */
        if (needstab == false)
        {  
#ifdef LOGFILE2
					/*logfile = fopen(LOGFILENAME,"a");
					if (logfile != NULL)
					{
						fprintf(logfile,"%d %d %d\n",n,tau,(1<<(tau+1))*l+(1<<tau));
						fprintf(logfile,"%d %d %d\n",n,tau,(1<<(tau+1))*l+(1<<tau)+1);
						fclose(logfile);
					}*/
					sprintf(filename,"%d__%d_%d_%d_1_u.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m1[j]);						  
						}
						fclose(logfile);
					}
					sprintf(filename,"%d__%d_%d_%d_2_u.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m2[j]);						  
						}
						fclose(logfile);
					}
					alpha--;
					beta--;
					gamma--;
					sprintf(filename,"%d__%d_%d_%d_3_u.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m3[j]);						  
						}
						fclose(logfile);
					}
					sprintf(filename,"%d__%d_%d_%d_4_u.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m4[j]);						  
						}
						fclose(logfile);
					}
#endif				  //printf("Alloc l: n = %d, tau = %d, l = %d\n",n,tau,l);
					//fflush(stdout);
          U[n][tau][l] = (struct U_type*) fftw_malloc(sizeof(struct U_type)); 
          /* No stabilization needed. */
          U[n][tau][l][0].m1 = m1;
          U[n][tau][l][0].m2 = m2;
          U[n][tau][l][0].m3 = m3;
          U[n][tau][l][0].m4 = m4;
          U[n][tau][l][0].stable = true;
        }          
        else 
        {    
#ifdef LOGFILE
					logfile = fopen(LOGFILENAME,"a");
					if (logfile != NULL)
					{
						fprintf(logfile,"%d %d %d\n",n,tau,(1<<(tau+1))*l+(1<<tau));
						fprintf(logfile,"%d %d %d\n",n,tau,(1<<(tau+1))*l+(1<<tau)+1);
						fclose(logfile);
					}
					eval_al(xvecs[tau-1], m1, plength, degree-2, alpha, beta, gamma);
					sprintf(filename,"%d__%d_%d_%d_1.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m1[j]);						  
						}
						fclose(logfile);
					}
					eval_al(xvecs[tau-1], m2, plength, degree-1, alpha, beta, gamma);
					sprintf(filename,"%d__%d_%d_%d_2.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m2[j]);						  
						}
						fclose(logfile);
					}
					alpha--;
					beta--;
					gamma--;
					eval_al(xvecs[tau-1], m3, plength, degree-1, alpha, beta, gamma);
					sprintf(filename,"%d__%d_%d_%d_3.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m3[j]);						  
						}
						fclose(logfile);
					}
					eval_al(xvecs[tau-1], m4, plength, degree, alpha, beta, gamma);
					sprintf(filename,"%d__%d_%d_%d_4.dat",tau,(1<<tau)-1,n,(1<<(tau+1))*l+1);
					logfile = fopen(filename,"w");
					if (logfile != NULL)
					{
					  for (j = 0; j < plength; j++)
						{
  						fprintf(logfile,"%.16E\n",m4[j]);						  
						}
						fclose(logfile);
					}
#endif

          nstab++;
          //fprintf(stderr,"(%d,%d,%d)\n",n,tau,l);
          
          /* Stabilize. */
  		      degree_stab = degree*(2*l+1);          
          
          /* Old arrays are to small. */
          fftw_free(m1);
          fftw_free(m2);
          fftw_free(m3);
          fftw_free(m4);

          if (window == true)
          {
            U[n][tau][l] = (struct U_type*) fftw_malloc(sizeof(struct U_type));             
            plength_stab = 1<<t;

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
            eval_al(xvecs[t-2], m1, plength_stab, degree_stab-2, alpha, 
                    beta, gamma);
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
            eval_al(xvecs[t-2], m2, plength_stab, degree_stab-1, alpha, 
                    beta, gamma);
            alpha--;
            beta--;
            gamma--;
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
            eval_al(xvecs[t-2], m3, plength_stab, degree_stab-1, alpha, 
                    beta, gamma);
            /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
            eval_al(xvecs[t-2], m4, plength_stab, degree_stab+0, alpha, 
                    beta, gamma);
            
            U[n][tau][l][0].m1 = m1;
            U[n][tau][l][0].m2 = m2;
            U[n][tau][l][0].m3 = m3;
            U[n][tau][l][0].m4 = m4;          
            U[n][tau][l][0].stable = false;                    
          }  
          else
          {
            U[n][tau][l] = (struct U_type*) fftw_malloc(sizeof(struct U_type)*(t-tau));             
            for (tau_stab = tau-1; tau_stab <= t-2; tau_stab++)
            {
              //tau_stab = t-2;
              plength_stab = 1<<(tau_stab+2);
              /* Allocate memory for arrays. */
              m1 = (double*) fftw_malloc(sizeof(double)*plength_stab);
              m2 = (double*) fftw_malloc(sizeof(double)*plength_stab);
              m3 = (double*) fftw_malloc(sizeof(double)*plength_stab);
              m4 = (double*) fftw_malloc(sizeof(double)*plength_stab);
							
							/*if (m1 == NULL || m2 == NULL || m3 == NULL || m4 == NULL)
							{
							  fprintf(stderr,"Precompute U: stabilized U -> malloc failed!");
								fflush(stderr);
							}*/
              
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
              
              U[n][tau][l][tau_stab-tau+1].m1 = m1;
              U[n][tau][l][tau_stab-tau+1].m2 = m2;
              U[n][tau][l][tau_stab-tau+1].m3 = m3;
              U[n][tau][l][tau_stab-tau+1].m4 = m4;          
              U[n][tau][l][tau_stab-tau+1].stable = false;                    
            }
          }
        }
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
    
  fprintf(stdout,"Stabilized %d times.\n",nstab);
  
  return U;
}


void forgetU(struct U_type**** U, int M, int t, bool window)
{ 
  const int N = 1<<t;
  /** Legendre index n */
  int n;
  /** Cascade level */
  int tau;
  /** Length of polynomials for the current level in the cascade */
  int plength;
  /** First index l for current cascade level and current n */
  int firstl;
  /** Last index l for current cascade level and current n */
  int lastl;
  /** Number of matrices U for current cascade level and current n .*/
  //int nsteps;
  int l;
  int tau_stab;
  int lb;
  int ntilde;
  int Mtilde = min(M,N-1);

  for (n = 0; n <= M; n++)
  {   
    ntilde = max(min(n,/*N*/(1<<ngpt(n))-2),0);
    plength = 4;
    for (tau = 1; tau < t; tau++)
    {
      /* Compute first l. */
      firstl = FIRST_L;//0; //1<<((int)log2(n)-(n==M?1:0)-tau-1);
      /* Compute last l. */
      lastl = LAST_L;//(int)(((double)(1<<t))/plength) - 1;
      /* Compute number of matrices for this level. */
      //nsteps = lastl - firstl + 1;
      /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
  	   for (l = firstl; l <= lastl; l++)
	     {        
        if (U[n][tau][l][0].stable == true)
        {  
	         fftw_free(U[n][tau][l][0].m1);
	         fftw_free(U[n][tau][l][0].m2);
	         fftw_free(U[n][tau][l][0].m3);
	         fftw_free(U[n][tau][l][0].m4);
        }
        else
        {
          if (window == true)
          {
            fftw_free(U[n][tau][l][0].m1);
            fftw_free(U[n][tau][l][0].m2);
            fftw_free(U[n][tau][l][0].m3);
            fftw_free(U[n][tau][l][0].m4);
          }
          else
          {
            for (tau_stab = tau-1; tau_stab <= t-2; tau_stab++)
            {
              fftw_free(U[n][tau][l][tau_stab-tau+1].m1);
              fftw_free(U[n][tau][l][tau_stab-tau+1].m2);
              fftw_free(U[n][tau][l][tau_stab-tau+1].m3);
              fftw_free(U[n][tau][l][tau_stab-tau+1].m4);
            }
          }
        }
        fftw_free(U[n][tau][l]);
      }
      fftw_free(U[n][tau]);
      /** Increase polynomial degree to next power of two. */
      plength = plength << 1;
    }
    fftw_free(U[n]);
  }
  fftw_free(U);
}


inline void multiplyU(complex  *a, complex *b, struct U_type u, int tau, int n, int l, 
               struct nfsft_wisdom *tw, double gamma)
{ 
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  double normalize = 1.0/(length<<1);
  
  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;   
  
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)b,(double*)b);
  
  /* Check, if gamma_k^n is zero. This is the case when l <= n holds. */
  if (gamma == 0.0)
  {
    /* Perform multiplication only for second row. */
    auvxpwy(normalize,b,b,u.m4,a,u.m3,length);
  }
  else 
  {
    /* Perform multiplication for both rows. */
    auvxpwy(normalize,tw->z,b,u.m4,a,u.m3,length);
    auvxpwy(normalize*gamma,a,a,u.m1,b,u.m2,length);
    memcpy(b,tw->z,length*sizeof(complex));    
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)a,(double*)a);   
    /* Compensate for factors introduced by a raw DCT-II. */    
    a[0] *= 0.5;
  }  

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)b,(double*)b);  
  /* Compensate for factors introduced by a raw DCT-II. */      
  b[0] *= 0.5;  
}


inline void multiplyU_adjoint(complex  *a, complex *b, 
                       struct U_type u, int tau, int n, int l, 
                       struct nfsft_wisdom *tw, double gamma)
{ 
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  double normalize = 1.0/(length<<1);
    
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)b,(double*)b);
  
  /* Perform matrix multiplication. */
  abuvxpwy(normalize,gamma,tw->z,a,u.m1,b,u.m3,length);
  abuvxpwy(normalize,gamma,b,a,u.m2,b,u.m4,length);
  memcpy(a,tw->z,length*sizeof(complex));
  
  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)a,(double*)a);   
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)b,(double*)b);   
}


struct U_type**** precomputeU_stab(int t, double threshold, double *walpha, 
                             double *wbeta, double *wgamma)
{
  const int N = 1<<t;
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
  //int nsteps;
  
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
  struct U_type ****U;
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
  //int i;
  /** Used to indicate that stabilization is neccessary. */
  bool needstab = false;  
  
  int ntilde;
  int Mtilde = min(M,N-1);
  
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
  U = (struct U_type****) fftw_malloc(sizeof(struct U_type ***) * (M+1));
  
  /* We have to precompute the matrices 
   *   U_{n,tau,l} := U_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1)
   * for n = 0,...,M; tau = 1,...,t and l = 0,...,2^{t-tau-1}-1. */ 
  
  m1 = (double*) fftw_malloc(sizeof(double)*M);
  m2 = (double*) fftw_malloc(sizeof(double)*M);
  m3 = (double*) fftw_malloc(sizeof(double)*M);
  m4 = (double*) fftw_malloc(sizeof(double)*M);                 
  
  /* For nleg = 0,...,M compute the matrices U_{n,tau,l}. */
  for (n = 0; n <= M; n++)
  {   
    ntilde = max(min(n,/*N*/(1<<ngpt(n))-2),0);
    fprintf(stderr,"n = %5d\n",n);
  		fflush(stderr);
    /* Allocate memory for current matrix array. The cascade will have 
     * t = log_2(M) many levels. */
    U[n] = (struct U_type***) fftw_malloc(sizeof(struct U_type**) * t);
    
    /* For tau = 1,...t compute the matrices U_{n,tau,l}. */
    plength = 4;
    for (tau = 1; tau < t; tau++)
    {     
      /* Compute auxilliary values. */
	     degree = plength>>1;
      
      /* Compute first l. */
      firstl = FIRST_L;//0; //1<<((int)log2(n)-(n==M?1:0)-tau-1);
      /* Compute last l. */
      lastl = LAST_L;//(int)(((double)(1<<t))/plength) - 1;
      /* Compute number of matrices for this level. */
      //nsteps = lastl - firstl + 1;
      
      /* Allocate memory for current matrix array. The level will contain
        * 2^{t-tau-1} many matrices. */
	     U[n][tau] = (struct U_type**) fftw_malloc(sizeof(struct U_type*) * (int)(((double)(1<<t))/plength)/*nsteps*/); 
      
      /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
  	   for (l = firstl; l <= lastl; l++)
	     {        

        //fprintf(stderr,"n = %d [%d,%d], tau = %d [%d,%d], l = %d [%d,%d]\n",n,0,M,tau,1,t-1,l,firstl,lastl);

        
        /* Allocate memory for the components of U_{n,tau,l}. */
	       /*m1 = (double*) fftw_malloc(sizeof(double)*plength);
	       m2 = (double*) fftw_malloc(sizeof(double)*plength);
	       m3 = (double*) fftw_malloc(sizeof(double)*plength);
	       m4 = (double*) fftw_malloc(sizeof(double)*plength);*/ 
        
        /* Evaluate the associated Legendre polynomials at the 2^{tau+1} 
         * Chebyshev nodes. */
        
        /* Get the pointers to the three-term recurrence coeffcients. */
        alpha = &(walpha[ROW(n)+plength*l+1+1]);
        beta = &(wbeta[ROW(n)+plength*l+1+1]);
        gamma = &(wgamma[ROW(n)+plength*l+1+1]);
        /* Evaluate P_{2^{tau}-2}^n(\cdot,2^{tau+1}l+2). */
        needstab = eval_al_thresh(xvecs[tau-1], m1, plength, degree-2, alpha, beta, gamma, threshold);
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
          U[n][tau][l] = (struct U_type*) fftw_malloc(sizeof(struct U_type)); 
          /* No stabilization needed. */
          /*U[n][tau][l][0].m1 = m1;
          U[n][tau][l][0].m2 = m2;
          U[n][tau][l][0].m3 = m3;
          U[n][tau][l][0].m4 = m4;*/
          U[n][tau][l][0].stable = true;
        }          
        else 
        {    
          //fprintf(stderr,"(%d,%d,%d)\n",n,tau,l);
          
          /* Stabilize. */
  		      //degree_stab = degree*(2*l+1);          
          
          /* Old arrays are to small. */
          /*fftw_free(m1);
          fftw_free(m2);
          fftw_free(m3);
          fftw_free(m4);*/

          U[n][tau][l] = (struct U_type*) fftw_malloc(sizeof(struct U_type)*(t-tau)); 

          for (tau_stab = tau-1; tau_stab <= t-2; tau_stab++)
          {
            //tau_stab = t-2;
            //plength_stab = 1<<(tau_stab+2);
            /* Allocate memory for arrays. */
            /*m1 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            m2 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            m3 = (double*) fftw_malloc(sizeof(double)*plength_stab);
            m4 = (double*) fftw_malloc(sizeof(double)*plength_stab);*/
            
            /* Get the pointers to the three-term recurrence coeffcients. */
            /*alpha = &(walpha[ROW(n)+2]);
            beta = &(wbeta[ROW(n)+2]);
            gamma = &(wgamma[ROW(n)+2]);         */
            /* Evaluate P_{2^{tau}(2l+1)-2}^n(\cdot,2). */
            /*eval_al(xvecs[tau_stab], m1, plength_stab, degree_stab-2, alpha, 
                    beta, gamma);*/
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
            /*eval_al(xvecs[tau_stab], m2, plength_stab, degree_stab-1, alpha, 
                    beta, gamma); 
            alpha--;
            beta--;
            gamma--;*/
            /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
            /*eval_al(xvecs[tau_stab], m3, plength_stab, degree_stab-1, alpha, 
                    beta, gamma);*/
            /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
            /*eval_al(xvecs[tau_stab], m4, plength_stab, degree_stab+0, alpha, 
                    beta, gamma);*/
                      
            /*U[n][tau][l][tau_stab-tau+1].m1 = m1;
            U[n][tau][l][tau_stab-tau+1].m2 = m2;
            U[n][tau][l][tau_stab-tau+1].m3 = m3;
            U[n][tau][l][tau_stab-tau+1].m4 = m4;*/
            U[n][tau][l][tau_stab-tau+1].stable = false;                    
          }
        }
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

void forgetU_stab(struct U_type**** U, int M, int t)
{ 
  const int N = 1<<t;
  /** Legendre index n */
  int n;
  /** Cascade level */
  int tau;
  /** Length of polynomials for the current level in the cascade */
  int plength;
  /** First index l for current cascade level and current n */
  int firstl;
  /** Last index l for current cascade level and current n */
  int lastl;
  /** Number of matrices U for current cascade level and current n .*/
  //int nsteps;
  int l;
  int tau_stab;
  int ntilde;
  int Mtilde = min(M,N-1);	
  
  for (n = 0; n <= M; n++)
  {   
    ntilde = max(min(n,/*N*/(1<<ngpt(n))-2),0);
    plength = 4;
    for (tau = 1; tau < t; tau++)
    {
      /* Compute first l. */
      firstl = FIRST_L;//0; //1<<((int)log2(n)-(n==M?1:0)-tau-1);
      /* Compute last l. */
      lastl = LAST_L;//(int)(((double)(1<<t))/plength) - 1;
        /* Compute number of matrices for this level. */
        //nsteps = lastl - firstl + 1;
        /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
        for (l = firstl; l <= lastl; l++)
        {        
          if (U[n][tau][l][0].stable == true)
          {  
            /*fftw_free(U[n][tau][l][0].m1);
            fftw_free(U[n][tau][l][0].m2);
            fftw_free(U[n][tau][l][0].m3);
            fftw_free(U[n][tau][l][0].m4);*/
          }
          else
          {
            for (tau_stab = tau-1; tau_stab <= t-2; tau_stab++)
            {
              /*fftw_free(U[n][tau][l][tau_stab-tau+1].m1);
              fftw_free(U[n][tau][l][tau_stab-tau+1].m2);
              fftw_free(U[n][tau][l][tau_stab-tau+1].m3);
              fftw_free(U[n][tau][l][tau_stab-tau+1].m4);*/
            }
          }
          fftw_free(U[n][tau][l]);
        }
        fftw_free(U[n][tau]);
        /** Increase polynomial degree to next power of two. */
        plength = plength << 1;
    }
    fftw_free(U[n]);
  }
  fftw_free(U);
}

#ifdef LOGFILE
inline void multiplyU_print(complex  *a, complex *b, struct U_type u, int tau, int n, int l, 
                      struct nfsft_wisdom *tw, double gamma, FILE *logfile, FILE *logfile2)
{ 
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  double normalize = 1.0/(length<<1);
  int j;
  
  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;   
  
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)b,(double*)b);
  
  for (j = 0; j < length; j++)
  {
    fprintf(logfile,"%.16E + %.16E * I\n",creal(a[j]),cimag(a[j]));
    fprintf(logfile2,"%.16E + %.16E * I\n",creal(b[j]),cimag(b[j]));
  }  
  
  /* Check, if gamma_k^n is zero. This is the case when l <= n holds. */
  if (gamma == 0.0)
  {
    /* Perform multiplication only for second row. */
    auvxpwy(normalize,b,b,u.m4,a,u.m3,length);
  }
  else 
  {
    /* Perform multiplication for both rows. */
    auvxpwy(normalize,tw->z,b,u.m4,a,u.m3,length);
    auvxpwy(normalize*gamma,a,a,u.m1,b,u.m2,length);
    memcpy(b,tw->z,length*sizeof(complex));    
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)a,(double*)a);   
    /* Compensate for factors introduced by a raw DCT-II. */    
    a[0] *= 0.5;
  }  
  
  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)b,(double*)b);  
  /* Compensate for factors introduced by a raw DCT-II. */      
  b[0] *= 0.5;  
}
#endif