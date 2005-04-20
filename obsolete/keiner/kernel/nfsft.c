#include "api.h"
#include "util.h"
#include "legendre.h"
#include "direct.h"

/* Global structure for precomputed values. */
static struct nfsft_wisdom wisdom = {false};
int nfft_size[2] = {0,0};
int fftw_size[2] = {0,0};

#define THRESHOLD 1000000

/**
 * Precomputation of the matrices U^n_k
 * 
 * \param level The number of levels
 */
struct U_type*** precomputeU(int M, double threshold, int t)
{ 
  int N;
  /** N minus 1 */
  int N1;
  int i;
  int l;
  int Ng;
  int m_stab;
  int firstl, lastl, plength, nsteps;
  /** Counter for M. */
  int n;
  /** Array that contains the (1,1)-components of matrices U. */
  double *m1;
  /** Array that contains the (1,1)-components of matrices U. */
  double *m2;
  /** Array that contains the (1,1)-components of matrices U. */
  double *m3;
  /** Array that contains the (1,1)-components of matrices U. */
  double *m4;
  /**  
    * Array for three-term recurrence coefficients of associated Legendre-
    * functions. 
    */
  double *alpha;
  /**  
    * Array for three-term recurrence coefficients of associated Legendre-
    * functions. 
    */
  double *beta;
  /**  
    * Array for three-term recurrence coefficients of associated Legendre-
    * functions. 
    */
  double *gamma;
  /** Array for Chebychev-nodes. */
  int tau;
  struct U_type*** U;
  double **xvecs;
  double *xc;
    
  /* Initialize arrays with Chebyshev coefficients for the polynomial x. This 
   * would be trivially an array containing a 1 as second entry with all other 
   * coefficients set to zero. In order to compensate for the multiplicative 
   * factor 2 introduced by the DCT-III, we set this coefficient to 0.5 here. */
  //printf("%d\n\n",1<<t); 
  xc = (double *) calloc(1<<t,sizeof(double));
  xc[1] = 0.5;
  
  /* Allocate memory for array of pointers. */
  xvecs = (double**) malloc((t-1)*sizeof(double*));
  /* For each polynomial length, compute the Chebyshev nodes using a DCT-III. */
  plength = 4;
  for (tau = 1; tau < t; tau++)
  {
    /* Allocate memory for vector. */
    xvecs[tau-1] = (double*) malloc(plength*sizeof(double));
    /* Create plan for DCT-III. */
    fftw_plan plan = fftw_plan_r2r_1d(plength, xc, xvecs[tau-1], FFTW_REDFT01, FFTW_PRESERVE_INPUT);
    /* Execute it. */
    fftw_execute(plan);
    /* Destroy the plan. */
    fftw_destroy_plan(plan);
    /* Increase length. */
    plength = plength << 1;
  }
  
  /** Used to indicate that stabilization is neccessary. */
  bool needstab = false;
  
  /* Allocate memory for matrices U^n_k(\cdot,4l+1). */
  U = (struct U_type***) fftw_malloc(sizeof(struct U_type **) * (M+1));
  
  for (n = 0; n <= M; n++)
  {   
    /* Allocate memory for current matrix array. */
    U[n] = (struct U_type**) fftw_malloc(sizeof(struct U_type*)*t);

    /* Compute U for all levels. */
    plength = 4;
    for (tau = 1; tau < t; tau++) /* m=2^{tau}-1 */
    {     
      /* Compute first l. */
      firstl = 0;//pow2((int)log2(n)-(n==M?1:0)-tau-1);
      /* Compute last l. */
      lastl = (int)(((double)pow2(t))/plength) - 1;//int)ceil(((double)M)/plength) - 1;
#ifdef DEBUG
      printf("N = %d, plength = %d\n",N,plength);
#endif
      /* Compute number of steps. */
      nsteps = lastl - firstl + 1;
      
      /* Allocate memory for matrix entry. */
	     U[n][tau] = (struct U_type*) fftw_malloc(sizeof(struct U_type) * nsteps); 
      
      /* Compute auxilliary values. */
	     N = plength>>1;
	     N1 = N-1;
      
      l = firstl;
  	   for (i = 0; i < nsteps; i++)
	     {
        /* Compute U_{2^tau-1}^n(x,pelgth*l+1) */
        
        /* Allocate memory for arrays. */
	       m1 = (double*) fftw_malloc(sizeof(double)*(plength+1));
	       m2 = (double*) fftw_malloc(sizeof(double)*(plength+1));
	       m3 = (double*) fftw_malloc(sizeof(double)*(plength+1));
	       m4 = (double*) fftw_malloc(sizeof(double)*(plength+1));               
        
        /* Compute the corresponding polynomials evaluated at the Chebyshev nodes. */
        alpha = &(wisdom.alpha[ROW(n)+plength*l+1+1]);
        beta = &(wisdom.beta[ROW(n)+plength*l+1+1]);
        gamma = &(wisdom.gamma[ROW(n)+plength*l+1+1]);
        needstab = eval_al_thresh(xvecs[tau-1], m1, plength, N1-1, alpha, beta, gamma, threshold);
        if (needstab == false)
        {
          needstab = eval_al_thresh(xvecs[tau-1], m2, plength, N1, alpha, beta, gamma, threshold);
          if (needstab == false)
          { 
            alpha--;
            beta--;
            gamma--;
            needstab = eval_al_thresh(xvecs[tau-1], m3, plength, N1, alpha, beta, gamma, threshold);
            if (needstab == false)
            { 
              needstab = eval_al_thresh(xvecs[tau-1], m4, plength, N1+1, alpha, beta, gamma, threshold);
            }
          }
        }      
        
        if (needstab == false)
        {  
          /* No stabilization needed. */
          U[n][tau][i].m1 = m1;
          U[n][tau][i].m2 = m2;
          U[n][tau][i].m3 = m3;
          U[n][tau][i].m4 = m4;
          U[n][tau][i].stable = true;
          //printf("U: n = %d, tau = %d, i = %d\n",n,tau,i);
        }          
        else 
        {    
          /* Stabilization needed. */
  		      m_stab = N1+plength*l+1-1; /* =j*dim+pow2(tau)-1; (alte Version)*/
          
		        if (plength == pow2(t)) 
		        {
		          Ng = plength; 
		        }
		        else
		        {
		          Ng = pow2(1 + (int) (log((double)(l*plength+N)))/log(2.0)); 
		          /* Ng = Zweierpotenz um schnelle Polynmomultiplikation auszuf¸hren */
          }
          
          /* Old arrays are to small. */
          fftw_free(m1);
          fftw_free(m2);
          fftw_free(m3);
          fftw_free(m4);
          
          /* Allocate memory for arrays. */
          m1 = (double*) fftw_malloc(sizeof(double)*(Ng+1));
          m2 = (double*) fftw_malloc(sizeof(double)*(Ng+1));
          m3 = (double*) fftw_malloc(sizeof(double)*(Ng+1));
          m4 = (double*) fftw_malloc(sizeof(double)*(Ng+1));
                    
          /* Compute the corresponding polynomials evaluated at the Chebyshev nodes. */
          alpha = &(wisdom.alpha[ROW(n)+2]);
          beta = &(wisdom.beta[ROW(n)+2]);
          gamma = &(wisdom.gamma[ROW(n)+2]);         
          eval_al(xvecs[(int)(log((double)Ng)/log(2.0))-2], m1, Ng, m_stab-1, alpha, beta, gamma);
          eval_al(xvecs[(int)(log((double)Ng)/log(2.0))-2], m2, Ng, m_stab, alpha, beta, gamma); 
          alpha--;
          beta--;
          gamma--;
          eval_al(xvecs[(int)(log((double)Ng)/log(2.0))-2], m3, Ng, m_stab, alpha, beta, gamma);
          eval_al(xvecs[(int)(log((double)Ng)/log(2.0))-2], m4, Ng, m_stab+1, alpha, beta, gamma);
          
          U[n][tau][i].m1=m1;
          U[n][tau][i].m2=m2;
          U[n][tau][i].m3=m3;
          U[n][tau][i].m4=m4;          
          U[n][tau][i].stable=false;                    
        }
        l++;
      }
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

/**
 * Fast matrix Multiplication with matrices U_k^n.
 *
 * \param a The first array of Chebyshev-coefficients
 * \param b The second array of Chebyshev-coefficients
 * \param u The 2x2-matrix U_k^n(.,l)
 * \param tau The parameter tau with k = 2^tau-1
 * \param n The parameter n
 * \param l The parameter l
 */
void multiplyU(complex  *a, complex *b, struct U_type u, int tau, int n, int l, 
               struct nfsft_transform_wisdom *tw, double gamma)
{ 
  /** Used to store temporary valiues. */
  complex z;
  /** Counter for loops. */
  int i;
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  int length2 = length<<1;

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
  
//#define DEBUG
#ifdef DEBUG
  printf("Before transform:\n");
  myvprc(a,length,"a");
  myvprc(tw->a2,length,"a2");
  myvprc(b,length,"b");
  myvprc(tw->b2,length,"b2");
#endif
  
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(tw->plans_dct3[tau-1],(double*)b,(double*)b);
  /* Make copies. */
#ifdef DEBUG
  printf("length = %d\n",length*sizeof(complex));
#endif
  memcpy(tw->a2,a,length*sizeof(complex));
  memcpy(tw->b2,b,length*sizeof(complex));
  
#ifdef DEBUG
  printf("After transform:\n");
  myvprc(a,length,"a");
  myvprc(tw->a2,length,"a2");
  myvprc(b,length,"b");
  myvprc(tw->b2,length,"b2");
#endif
  
  //printf("reached\n");
  for (i = 0; i < length; i++)
  {
    //printf("i = %d\n",i);
    //printf("u.m1[%d] = %d\n",i,u.m2);
    a[i] *= gamma * u.m1[i];
    tw->a2[i] *= gamma * u.m2[i];
    b[i] *= u.m3[i];
    tw->b2[i] *= u.m4[i];
  }    
  
  for (i = 0; i < length; i++)
  {
    a[i] += b[i];
    tw->a2[i] += tw->b2[i];
  }    
  
  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)a,(double*)a);   
  fftw_execute_r2r(tw->plans_dct2[tau-1],(double*)tw->a2,(double*)tw->a2);   
  /* Compensate for lack of normalization of DCT-II. */          
  for (i = 0; i < length; i++)
  {
    a[i] = a[i]/length2;
    tw->a2[i] = tw->a2[i]/length2;
  }  
  
  memcpy(&a[length/2],tw->a2,(length/2)*sizeof(complex));
}

/**
 * Fast Legendre Function Transform
 *
 * \param 
 */

//#define DEBUG

void flft(struct U_type ***U, int M, int t, int n, complex *daten, 
         struct nfsft_transform_wisdom *tw)
{
  int plength, nsteps, firstl, lastl, tau;
  int N, N1, Ng, tg;
  double alpha, beta, gamma, gammaconst;
  int i,k,l,j;
  int coeff_index;
  struct U_type act_U;
  int rindex = ROWK(n);
  double *ngamma = &(wisdom.gamma[ROWK(n)]);

  /* Calculate auxilliary values. */
  N = pow2(t);
  N1 = N + 1;
  
  /* Initialize working arrays. */
  memset(tw->work,0U,(N1<<1)*sizeof(complex));
  memset(tw->ergeb,0U,(N1<<1)*sizeof(complex));
  /* Set first n Fourier-coefficients explicitely to zero. They should be 
   * zero right from the beginning anyway. */
  memset(daten,0U,n*sizeof(complex));
  
#ifdef DEBUG
  printf("FLFT: Before first step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  /* First step */ 
  for (k = 0; k < N; k++) 
  {
    tw->work[2*k] = daten[k];
  }
  
#ifdef DEBUG
  printf("FLFT: After copy\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  /* Use three-term recurrence to map last coefficient a_N to a_{N-1} and 
    * a_{N-2}. */
  coeff_index = N - n;
  alpha = wisdom.alpha[rindex+coeff_index];
  beta = wisdom.beta[rindex+coeff_index];
  gamma = wisdom.gamma[rindex+coeff_index];
  
#ifdef DEBUG
  printf("alpha_{N-1 = %d}^{|n| = %d} = %f\n",N-1,n,alpha);
  printf("beta_{N-1 = %d}^{|n| = %d} = %f\n",N-1,n,beta);
  printf("gamma_{N-1 = %d}^{|n| = %d} = %f\n",N-1,n,gamma);
  printf("\n");
#endif
  
  tw->work[2*(N-2)] += gamma*daten[N];
  tw->work[2*(N-1)] = daten[N-1] + beta*daten[N];
  tw->work[2*(N-1)+1] = alpha*daten[N];
  
#ifdef DEBUG
  printf("FLFT: After map\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  
  /* Compute the remaining steps. */
  plength = 4;
  for (tau = 1; tau < t; tau++)
  {    
#ifdef DEBUG
    printf("FLFT: \tau = %d\n",tau);
    printf("FLFT: plength = %d\n",plength);
#endif    
    /* Compute first l. */
    firstl = pow2((int)(log((double)n)/log(2.0))-(n==N?1:0)-tau-1);
    /* Compute last l. */
    lastl = (int)ceil(((double)M)/plength) - 1;
    /* Compute number of steps. */
    nsteps = lastl - firstl + 1;
    
    /* Compute the multiplication steps. */
    l = firstl;
    for (i = 0; i < nsteps; i++)
    {  
      /* Initialize second half of coefficient arrays with zeros. */
      memset(&tw->vec1[plength/2],0U,plength/2*sizeof(complex));
      memset(&tw->vec2[plength/2],0U,plength/2*sizeof(complex));
      memset(&tw->vec3[plength/2],0U,plength/2*sizeof(complex));
      memset(&tw->vec4[plength/2],0U,plength/2*sizeof(complex));
      
      /* Copy coefficients into first half. */
      memcpy(tw->vec1,&(tw->work[plength*(4*l+0)/2]),plength/2*sizeof(complex));
      memcpy(tw->vec2,&(tw->work[plength*(4*l+1)/2]),plength/2*sizeof(complex));
      memcpy(tw->vec3,&(tw->work[plength*(4*l+2)/2]),plength/2*sizeof(complex));
      memcpy(tw->vec4,&(tw->work[plength*(4*l+3)/2]),plength/2*sizeof(complex));     
      
      /* Get matrix U_(2^tau-1)^n() */
      act_U = U[n][tau][firstl+i];
      
      /* Check if step is stable. */
      if (act_U.stable)
      {
        gammaconst = ngamma[plength*l+1-n+1];            
        /* Multiply third and fourth polynomial with matrix U. */
        multiplyU(tw->vec3, tw->vec4, act_U, tau, n, plength*l+1, tw, gammaconst);        
        for (j = 0; j < plength; j++)
        {
          tw->work[plength*2*l+j] = tw->vec1[j] += tw->vec3[j];
          tw->work[plength*(2*l+1)+j] = tw->vec2[j] += tw->vec4[j];
        }          
      }
      else
      {	
        /* Stabilize. */
        if (plength == N) 
        {
          Ng = plength;
          tg = tau;
        }  
        else 
        {
          tg = (int) (log((double)(l*plength+pow2(tau)))/log(2.0));
          Ng = pow2(tg+1); 
        }
        /* Ng = Zweierpotenz um schnelle Polynomultiplikation auszuf¸hren */
        /* Vektor auf die volle Laenge verlaengern */
        memset(&tw->vec3[plength/2],0U,(Ng-plength/2)*sizeof(complex));
        memset(&tw->vec4[plength/2],0U,(Ng-plength/2)*sizeof(complex));
        
        multiplyU(tw->vec3, tw->vec4, act_U, tg, n, 1, tw, 0.0);
        
        for (j = 0; j < Ng; j++)
        {
          tw->ergeb[N+j] += tw->vec4[j];
        }
        
        /* Don't change result. */
        memcpy(&(tw->work[plength*2*l]),tw->vec1,plength*sizeof(complex));
        memcpy(&(tw->work[plength*(2*l+1)]),tw->vec2,plength*sizeof(complex));        
      }
      l++;
    }
    plength = plength<<1;
  } 
  
#ifdef DEBUG
  printf("FLFT: Before final copy\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  // Copy result.
  for (i = 0; i < 2*N1; i++)
  {
    tw->ergeb[i] += tw->work[i];
  }  
#ifdef DEBUG
  printf("FLFT: Before last step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  /* The final step */ 
  gamma = wisdom.gamma_m1[n];
  
#ifdef DEBUG
  printf("FLFT: gamma_0^{|n| = %d} = %f\n\n",n,gamma);
#endif
  
  /* funktioniert nur fuer diese alphas !!*/
  if (n%2 == 0)
  {
    if (n == 0)
    {
      daten[0] = gamma*(tw->ergeb[0]+tw->ergeb[N+1]*0.5);
      daten[1] = gamma*(tw->ergeb[1]+(tw->ergeb[N]+tw->ergeb[N+2]*0.5));
      daten[N-1] = gamma*(tw->ergeb[N-1]+tw->ergeb[N+N-2]*0.5);
      daten[N] = gamma*(tw->ergeb[N+N-1]*0.5);
      for (i = 2; i < N-1; i++)
      {
        daten[i] = gamma*(tw->ergeb[i]+(tw->ergeb[i+N-1]+tw->ergeb[i+N+1])*0.5);
      } 
    }
    else
    {
      daten[0] = gamma*(tw->ergeb[0]+tw->ergeb[N]-tw->ergeb[N+1]*0.5);
      daten[1] = gamma*(tw->ergeb[1]+tw->ergeb[N+1]-(tw->ergeb[N]+tw->ergeb[N+2]*0.5));
      daten[N-1] = gamma*(tw->ergeb[N-1]+tw->ergeb[N+N-1]-tw->ergeb[N+N-2]*0.5);
      daten[N] = gamma*(tw->ergeb[N+N]-tw->ergeb[N+N-1]*0.5);
      for (i = 2; i < N-1; i++)
      {
        daten[i] = gamma*(tw->ergeb[i]+tw->ergeb[i+N]-(tw->ergeb[i+N-1]+tw->ergeb[i+N+1])*0.5);
      } 
    }
  }
  else
  {
    daten[N] = 0.0;
    for (i = 0; i < N; i++)
    {
      daten[i] = gamma*(tw->ergeb[i]+tw->ergeb[i+N]);
    }
  }
  
#ifdef DEBUG
  printf("FLFT: After last step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
}

/**
* Fast Legendre Function Transform
 *
 * \param 
 */
void flft_adjoint(struct U_type ***U, int M, int t, int n, complex *daten, 
          struct nfsft_transform_wisdom *tw)
{
  int plength, nsteps, firstl, lastl, tau;
  int N, N1, Ng, tg;
  double alpha, beta, gamma, gammaconst;
  int i,k,l,j;
  int coeff_index;
  struct U_type act_U;
  int rindex = ROWK(n);
  double *ngamma = &(wisdom.gamma[ROWK(n)]);
  
  /* Calculate auxilliary values. */
  N = pow2(t);
  N1 = N + 1;
  
  /* Initialize working arrays. */
  memset(tw->work,0U,(N1<<1)*sizeof(complex));
  memset(tw->ergeb,0U,(N1<<1)*sizeof(complex));
  /* Set first n Fourier-coefficients explicitely to zero. They should be 
    * zero right from the beginning anyway. */
  //memset(daten,0U,n*sizeof(complex));
  
//#define DEBUG
#ifdef DEBUG
  printf("FLFT-adjoint: Before 'last' step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  /* The final step */ 
  gamma = wisdom.gamma_m1[n];
  
  /* funktioniert nur fuer diese alphas !!*/
  /* Multiplication with T_0^T */
  
#ifdef DEBUG
  printf("FLFT: gamma_0^{|n| = %d} = %f\n\n",n,gamma);
#endif
  
  /* First half consists always of coefficient vector multiplied by I_{N+1}, 
   * i.e. a copy of this vector. */
  for (i = 0; i <= N; i++)
  {
    tw->work[i] = gamma*daten[i]; 
  }
  
  /* Distinguish by n for the second half. */
  if (n%2 == 0)
  {
    if (n == 0)
    {
      /* Second half is T_{N+1}^T */
      tw->work[N+1+0] = gamma * daten[1];
      for (i = 1; i < N; i++)
      {
        tw->work[N+1+i] = gamma*0.5*(daten[i-1] + daten[i+1]);
      } 
      tw->work[N+1+N] = 0.5*gamma*daten[N-1];     
    }
    else
    {
      /* Second half is I_{N+1} - T_{N+1}^T */
      tw->work[N+1+0] = gamma * (daten[0] - daten[1]);
      for (i = 1; i < N; i++)
      {
        tw->work[N+1+i] = gamma * (daten[i] - 0.5*(daten[i-1] + daten[i+1]));
      } 
      tw->work[N+1+N] = gamma * (daten[N] - 0.5*daten[N-1]);     
    }
  }
  else
  {
    /* Second half is I_{N+1} */    
    for (i = 0; i <= N; i++)
    {
      tw->work[N+1+i] = gamma*daten[i];
    }
  }  
  
#ifdef DEBUG
  printf("FLFT-adjoint: after 'last' step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  memmove(&tw->work[N],&tw->work[N+1],N*sizeof(complex));
  memset(&tw->work[2*N],0U,2*sizeof(complex));
  
//#define DEBUG
#ifdef DEBUG
  for (k = 0; k < 2*N; k++)
  {
    printf("%17.16f\n%17.16f\n",creal(tw->work[k]),cimag(tw->work[k]));
  }
#endif
#undef DEBUG  
  
  complex *old = malloc(2*N*sizeof(complex));
  memcpy(old,tw->work,2*N*sizeof(complex));
  
#ifdef DEBUG
  printf("FLFT-adjoint: After 'last' step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  /* Compute the remaining steps. */
  plength = N;
  
  if (true/*n == 0*/)
  {
    
  for (tau = t-1; tau >= 1; tau--)
  {    
#ifdef DEBUG
    printf("tau = %d\n",tau);
#endif
    /* Compute first l. */
    //firstl = pow2((int)log2(n)-(n==N?1:0)-tau-1);
    firstl = 0;    
    /* Compute last l. */
    lastl = (int)(((double)N)/plength) - 1;
    
    /* Compute number of steps. */
    nsteps = lastl - firstl + 1;
    
    /* Compute the multiplication steps. */
    l = firstl;
    for (i = 0; i < nsteps; i++)
    {  
#ifdef DEBUG
      printf("i = %d\n",i);
#endif
      /* Initialize second half of coefficient arrays with zeros. */
      memcpy(tw->vec3,&(tw->work[plength*(4*l+0)/2]),plength*sizeof(complex));
      memcpy(tw->vec4,&(tw->work[plength*(4*l+2)/2]),plength*sizeof(complex));     

//#define DEBUG
#ifdef DEBUG
      if (n == 0)
      {
        printf("-----\n");
        for (j = 0; j < plength; j++)
        {
          printf("%f\n",creal(tw->vec3[j]));
        }  
        printf("--\n");
        for (j = 0; j < plength; j++)
        {
          printf("%f\n",creal(tw->vec4[j]));
        }  
        printf("-----\n");        
      }
#endif      
#undef DEBUG

      //memcpy(&tw->work[plength*(4*l+0)/2],&(tw->work[plength*(4*l+0)/2]),(plength/2)*sizeof(complex));
      memcpy(&tw->work[plength*(4*l+1)/2],&(tw->work[plength*(4*l+2)/2]),(plength/2)*sizeof(complex));
      
      
#ifdef DEBUG
      printf("Before U:\n");
      myvprc(tw->vec3,plength,"vec3:");
      myvprc(tw->vec4,plength,"vec4:");
      myvprc(tw->work,2*N1,"work:");
      printf("\n");
#endif
      
      /* Get matrix U_(2^tau-1)^n() */
      act_U = U[n][tau][i];
      
      /* Check if step is stable. */
      if (act_U.stable)
      {
        /* Multiply third and fourth polynomial with matrix U. */
        gammaconst = ngamma[plength*l+1-n+1];        
        multiplyU_adjoint(tw->vec3, tw->vec4, act_U, tau, n, plength*l+1, tw, gammaconst);
        
//#define DEBUG
#ifdef DEBUG
        printf("-----\n");
        for (j = 0; j < plength; j++)
        {
          printf("%f\n",gammaconst*act_U.m1[j]);
        }  
        for (j = 0; j < plength; j++)
        {
          printf("%f\n",gammaconst*act_U.m2[j]);
        }  
        for (j = 0; j < plength; j++)
        {
          printf("%f\n",act_U.m3[j]);
        }  
        for (j = 0; j < plength; j++)
        {
          printf("%f\n",act_U.m4[j]);
        }  
        printf("-----\n");
#endif
#undef DEBUG
        
#ifdef DEBUG
        printf("After U:\n");
        myvprc(daten,N1,"Daten:");
        myvprc(tw->ergeb,2*N1,"ergeb:");
        myvprc(tw->work,2*N1,"work:");
        myvprc(tw->vec1,plength,"vec1:");
        myvprc(tw->vec2,plength,"vec2:");
        myvprc(tw->vec3,plength,"vec3:");
        myvprc(tw->vec4,plength,"vec4:");
        printf("\n");
#endif

//#define
#ifdef DEBUG
          printf("gamma = %f\n",gammaconst);
#endif
          for (j = 0; j < plength; j++)
          {
            tw->work[plength*(4*l+2)/2+j] = tw->vec3[j];
            //tw->work[plength+j] = gammaconst * tw->vec3[j] + tw->vec4[j];
          }

//#define DEBUG
#ifdef DEBUG
          if (n == 0)
          {
            printf("-----\n");
            for (j = 0; j < 2*plength; j++)
            {
              printf("%f + %f i\n",creal(tw->work[plength*(4*l+0)/2+j]),cimag(tw->work[plength*(4*l+0)/2+j]));
            }  
            printf("-----\n");        
          }
#endif      
#undef DEBUG          
          
#ifdef DEBUG
          printf("After copy:\n");
          myvprc(tw->vec1,plength,"vec1:");
          myvprc(tw->vec2,plength,"vec2:");
          myvprc(tw->vec3,plength,"vec3:");
          myvprc(tw->vec4,plength,"vec4:");
          myvprc(tw->work,(N1<<1),"work:");
          myvprc(tw->ergeb,(N1<<1),"ergeb:");
          myvprc(daten,N1,"daten:");
          printf("\n");
#endif
        //}        
      }
      else
      {	
        //printf("must stabilize!!!!!\n");
        /* Stabilize. */
        //printf("Warning: Not implemented yet!\n");
        if (plength == N) 
        {
          Ng = plength;
          tg = tau;
        }  
        else 
        {
          tg = (int) (log((double)(l*plength+pow2(tau)))/log(2.0));
          Ng = pow2(tg+1); 
        }
        /* Ng = Zweierpotenz um schnelle Polynomultiplikation auszuf¸hren */
        /* Vektor auf die volle Laenge verlaengern */
        //printf("before\n");
        memcpy(tw->vec3,old,N*sizeof(complex));
        memcpy(tw->vec4,&old[N],N*sizeof(complex));
        //printf("after\n");
        //memset(&tw->vec3[plength/2],0U,(Ng-plength/2)*sizeof(complex));
        //memset(&tw->vec4[plength/2],0U,(Ng-plength/2)*sizeof(complex));
        
        //printf("n = %d, tau = %d, i = %d\n",n,tau,i);
        //multiplyU_adjoint(tw->vec3, tw->vec4, act_U, tau, n, plength*l+1, tw, gammaconst);
          multiplyU_adjoint(tw->vec3, tw->vec4, act_U,  tg, n,           1, tw, 0.0);
        //printf("next\n");

        for (j = 0; j < plength; j++)
        {
          tw->work[plength*(4*l+2)/2+j] = tw->vec3[j];
          //tw->work[plength+j] = gammaconst * tw->vec3[j] + tw->vec4[j];
        }
        
        /*for (j = 0; j < Ng; j++)
        {
          tw->ergeb[N+j] += tw->vec4[j];
        }*/
        
        /* Don't change result. */
        /*memcpy(&(tw->work[plength*2*l]),tw->vec1,plength*sizeof(complex));
        memcpy(&(tw->work[plength*(2*l+1)]),tw->vec2,plength*sizeof(complex));*/        
	     }
      l++;
    }
    plength = plength>>1;
  } 
    
  }
    
#ifdef DEBUG
  printf("FLFT-adjoint: Before 'first' step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  /* First step */ 
  memset(daten,0U,N1*sizeof(complex));
  for (k = 0; k < N; k++) 
  {
    daten[k] = tw->work[2*k];
  }
  
  coeff_index = N - n;
  alpha = wisdom.alpha[rindex+coeff_index];
  beta = wisdom.beta[rindex+coeff_index];
  gamma = wisdom.gamma[rindex+coeff_index];

//#define DEBUG
#ifdef DEBUG
  printf("alpha_{N-1 = %d}^{|n| = %d} = %f\n",N-1,n,alpha);
  printf("beta_{N-1 = %d}^{|n| = %d} = %f\n",N-1,n,beta);
  printf("gamma_{N-1 = %d}^{|n| = %d} = %f\n",N-1,n,gamma);
  printf("\n");
#endif
#undef DEBUG
  
  daten[N] = gamma*tw->work[2*N-4] + beta*tw->work[2*N-2] + alpha*tw->work[2*N-1];
  
  //memset(daten,0U,n*sizeof(complex));
  
#ifdef DEBUG
  printf("FLFT-adjoint: After 'first' step\n");
  myvprc(daten,N1,"Daten:");
  myvprc(tw->ergeb,2*N1,"ergeb:");
  myvprc(tw->work,2*N1,"work:");
  printf("\n");
#endif
  
  free(old);

//#define DEBUG  
#ifdef DEBUG
  for (k = 0; k < N+1; k++)
  {
    printf("%17.16f\n%17.16f\n",creal(daten[k]),cimag(daten[k]));
  }
#endif
#undef DEBUG  
  
//#define DEBUG  
#ifdef DEBUG
  for (k = 0; k < N+1; k++)
  {
    printf("%17.16f + %17.16f i\n",creal(daten[k]),cimag(daten[k]));
  }
#endif
#undef DEBUG  
}

/**
 * Converts Chebyshev coefficients to Fourier coefficients.
 */
void cheb2exp(complex *f_hat, complex **DATA,int M,int N)
{
  int k,n;
  complex *data;
  complex last, act;
  complex *f_hat_p;
  complex *f_hat_n;
  int l1,l2,u1,u2;
  int dim, dimh;
  int rowz;
  int colz;
  dimh = N + 1;
  dim = dimh<<1;
  rowz = (dimh)*dim;
  colz = dimh;

  if (M%2 == 0)
  {
    l1 = -M;
    u1 = M;
    l2 = -M+1;
    u2 = M-1;
  }
  else
  {
    l2 = -M;
    u2 = M;
    l1 = -M+1;
    u1 = M-1;
  }
  
  memset(f_hat,0U,dim*dim*sizeof(complex));
  
  /* Process even terms. */
  for (n = l1; n <= u1; n += 2)
  {
    data = DATA[n+M];
    
    f_hat[rowz+n*dim+colz] = *(data++);
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *f_hat_p++ = (*(data))/2.0;
      *f_hat_n-- = (*(data++))/2.0;
    }
  }
  
  /* Process odd terms. */
  for (n = l2; n <= u2; n += 2)
  {
    data = DATA[n+M];
    
    f_hat[rowz+n*dim+colz] = *(data++);
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *f_hat_p++ = (*(data))/2.0;
      *f_hat_n-- = (*(data++))/2.0;
    }
    
    /* Incorporate sine term. */
    last = f_hat[rowz+n*dim+colz-N];
    f_hat[rowz+n*dim+colz-N] = I * f_hat[rowz+n*dim+colz-N+1]/2.0;
    for (k = -N+1; k <= N-1; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = I * (f_hat[rowz+n*dim+colz+k+1] - last)/2.0;       
      last = act;
    }
    f_hat[rowz+n*dim+colz+N] = - I * last/2.0;
    /*last = f_hat[rowz+n*dim+colz-N-1];
    f_hat[rowz+n*dim+colz-N-1] = I * f_hat[rowz+n*dim+colz-N]/2.0;
    for (k = -N; k <= N; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = I * (f_hat[rowz+n*dim+colz+k+1] - last)/2.0;       
      last = act;
    }
    f_hat[rowz+n*dim+colz+N+1] = - I * last/2.0;*/
  }  
}


/**
* Converts Chebyshev coefficients to Fourier coefficients.
 */
void cheb2exp_adjoint(complex *f_hat, complex **DATA,int M,int N)
{
  int k,n;
  complex *data;
  complex last, act;
  complex *f_hat_p;
  complex *f_hat_n;
  int l1,l2,u1,u2;
  int dim, dimh;
  int rowz;
  int colz;
  
  dimh = N + 1;
  dim = dimh<<1;
  rowz = (dimh)*dim;
  colz = dimh;
  
  if (M%2 == 0)
  {
    l1 = -M;
    u1 = M;
    l2 = -M+1;
    u2 = M-1;
  }
  else
  {
    l2 = -M;
    u2 = M;
    l1 = -M+1;
    u1 = M-1;
  }
  
  /* Process even terms. */
  for (n = l1; n <= u1; n += 2)
  {
//#define DEBUG
#ifdef DEBUG
    printf("n = %d\n",n);
    for (k = -N; k<= N; k++)
    {
      printf("%f + %f i\n",creal(f_hat[rowz+n*dim+colz+k]),cimag(f_hat[rowz+n*dim+colz+k]));
    }
    printf("------\n");
#endif
    
    data = DATA[n+M];
	
	/*printf("data = %p\n",data);
	printf("N = %d\n",N);
	printf("length = %d\n",(N+1)*sizeof(complex));*/
    memset(data, 0U, (N+1)*sizeof(complex));
    
    *data++ = f_hat[rowz+n*dim+colz];
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *data++ = (*f_hat_p++ + *f_hat_n--)/2.0;
    }
  }
  
  /* Process odd terms. */
  for (n = l2; n <= u2; n += 2)
  {    
#ifdef DEBUG
    printf("n = %d\n",n);
    for (k = -N; k<= N; k++)
    {
      printf("%f + %f i\n",creal(f_hat[rowz+n*dim+colz+k]),cimag(f_hat[rowz+n*dim+colz+k]));
    }
    printf("------\n");
#endif
    
    /* Incorporate sine term. */
    last = f_hat[rowz+n*dim+colz-N];
    f_hat[rowz+n*dim+colz-N] = I * f_hat[rowz+n*dim+colz-N+1]/2.0;
    for (k = -N+1; k <= N-1; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = I * (f_hat[rowz+n*dim+colz+k+1] - last)/2.0;       
      last = act;
    }
    f_hat[rowz+n*dim+colz+N] = -I * last/2.0;
    /*last = f_hat[rowz+n*dim+colz-N-1];
    f_hat[rowz+n*dim+colz-N-1] = I * f_hat[rowz+n*dim+colz-N]/2.0;
    for (k = -N; k <= N; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = I * (f_hat[rowz+n*dim+colz+k+1] - last)/2.0;       
      last = act;
    }
    f_hat[rowz+n*dim+colz+N+1] = -I * last/2.0;*/
    
//#define DEBUG
#ifdef DEBUG
    printf("n = %d\n",n);
    for (k = -N; k<= N; k++)
    {
      printf("%f + %f i\n",creal(f_hat[rowz+n*dim+colz+k]),cimag(f_hat[rowz+n*dim+colz+k]));
    }
    printf("------\n");
#endif
#undef DEBUG

    data = DATA[n+M];
    memset(data, 0U, (N+1)*sizeof(complex));
    
    *data++ = f_hat[rowz+n*dim+colz];
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *data++ = (*f_hat_p++ + *f_hat_n--)/2.0;
    }
  }  
}

// --- Wisdom stuff ---

struct nfsft_transform_wisdom* init_transform_wisdom(int N, double threshold, 
                                                     int t)
{  
  int i;
  int ti;
  int N1;
  struct nfsft_transform_wisdom *tw;
  
  N1 = N+1;
  tw = malloc(sizeof(struct nfsft_transform_wisdom));
  wisdom.transform_wisdoms[t] = tw;
  tw->N = N;
  tw->t = t;
  
  tw->work = (complex*) calloc(3*N1,sizeof(complex));   
  tw->ergeb = (complex*) calloc(3*N1,sizeof(complex)); /* hier werden schrittweise die Cheb.-koeffizienten aufgebaut */
  tw->vec1 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->vec2 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->vec3 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->vec4 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->a2 = (complex*) fftw_malloc(sizeof(complex)*N1);
  tw->b2 = (complex*) fftw_malloc(sizeof(complex)*N1);
  
  tw->kinds = malloc(2*sizeof(fftw_r2r_kind));
  tw->kinds[0] = FFTW_REDFT01;
  tw->kinds[1] = FFTW_REDFT01;
  tw->kindsr = malloc(2*sizeof(fftw_r2r_kind));
  tw->kindsr[0] = FFTW_REDFT10;
  tw->kindsr[1] = FFTW_REDFT10;
  
  tw->plans_dct3 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(t-1));
  tw->plans_dct2 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(t-1));
  tw->lengths = (int *) malloc((t-1)*sizeof(int));
  for(i=0,ti=4;i<t-1;i++,ti<<=1)
  {
    tw->lengths[i] = ti;
    tw->plans_dct3[i] = fftw_plan_many_r2r(1, &tw->lengths[i], 2, (double*)tw->vec1, NULL, 2, 1,
                                           (double*)tw->vec1, NULL, 2, 1, tw->kinds, 0);
    tw->plans_dct2[i] = fftw_plan_many_r2r(1, &tw->lengths[i], 2, (double*)tw->vec1, NULL, 2, 1,
                                           (double*)tw->vec1, NULL, 2, 1, tw->kindsr, 0);
  }

  tw->U = precomputeU(N, threshold, t);
    
  return tw;
}

void forget_transform_wisdom(struct nfsft_transform_wisdom *tw)
{
  int ts,mu;
  int i,i_j,Ni,nleg;
  
  for (nleg = 0; nleg <= tw->N; nleg++)
  {
    if (nleg > 1)
    { 
      mu = (int) (log((double)nleg)/log(2.0));
    }  
    
    for (ts = 1; ts < tw->t; ts++)
	   {
	     if (nleg > 1) 
      {
        i_j = pow2(mu-ts-1);
      }  
	     else 
      {
        i_j = 0;
      }  
      
	     Ni = (int) ceil((double)tw->N/(double)pow2(ts+1)) - i_j;
      
      for (i = 0; i < Ni; i++)
	     {
	       fftw_free(tw->U[nleg][ts][i].m1);
	       fftw_free(tw->U[nleg][ts][i].m2);
	       fftw_free(tw->U[nleg][ts][i].m3);
	       fftw_free(tw->U[nleg][ts][i].m4);
      }
      
      fftw_free(tw->U[nleg][ts]);
    }
    
    fftw_free(tw->U[nleg]);
  }
  
  fftw_free(tw->U);
  
  fftw_free(tw->work);
  fftw_free(tw->ergeb);
  fftw_free(tw->vec1);
  fftw_free(tw->vec2);
  fftw_free(tw->vec3);
  fftw_free(tw->vec4);
  fftw_free(tw->a2);
  fftw_free(tw->b2);
  
  for(i = 0; i < tw->t-1; i++)
  {
    fftw_destroy_plan(tw->plans_dct3[i]);
    fftw_destroy_plan(tw->plans_dct2[i]);
  }  
  free(tw->plans_dct3);
  free(tw->plans_dct2);
  free(tw->kinds);
  free(tw->kindsr);
  free(tw->lengths);
}

inline void init_wisdom()
{ 
  if (wisdom.initialized == false)
  {
    wisdom.initialized = true;
    wisdom.transform_wisdoms = (struct nfsft_transform_wisdom**) 
      calloc(BWEXP_MAX+1,sizeof(struct nfsft_transform_wisdom*));
    
    alpha_al_all(wisdom.alpha,BW_MAX);
    beta_al_all(wisdom.beta,BW_MAX);
    gamma_al_all(wisdom.gamma,BW_MAX);
    gamma_al_m1_all(wisdom.gamma_m1,BW_MAX);
  }  
}

void nfsft_forget_wisdom()
{
  static int i;
  for (i = 0; i <= BWEXP_MAX; i++)
  {
    if (wisdom.transform_wisdoms[i] != 0)
    {  
      forget_transform_wisdom(wisdom.transform_wisdoms[i]);
    }
  }
  free(wisdom.transform_wisdoms);
  wisdom.initialized = false;
}

/*void export_transform_wisdom(struct nfsft_transform_wisdom *tw, FILE *f)
{
  int ts,mu;
  int i,i_j,Ni,nleg;
  int N,N2;
  
  for (nleg = 0; nleg <= tw->N; nleg++)
  {
    if (nleg > 1)
    { 
      mu = (int) (log2(nleg));
    }  
    
    for (ts = 1; ts < tw->t; ts++)
	   {
	     if (nleg > 1) 
      {
        i_j = pow2(mu-ts-1);
      }  
	     else 
      {
        i_j = 0;
      }  
      
	     Ni = (int) ceil((double)tw->N/(double)pow2(ts+1)) - i_j;
      
      N = pow2(ts);
	     N2 = 2*N;
      
      for (i = 0; i < Ni; i++)
	     {
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
        fwrite(tw->U[nleg][ts][i].m1,sizeof(double),N2+1,f);
      }
    }
  }
  
  fftw_free(tw->work);
  fftw_free(tw->ergeb);
  fftw_free(tw->vec1);
  fftw_free(tw->vec2);
  fftw_free(tw->vec3);
  fftw_free(tw->vec4);
  
  for(i = 0; i < tw->t-1; i++)
  {
    fftw_destroy_plan(tw->plans_dct3[i]);
    fftw_destroy_plan(tw->plans_dct2[i]);
  }  
  free(tw->plans_dct3);
  free(tw->plans_dct2);
  free(tw->kinds);
  free(tw->kindsr);
  free(tw->lengths);
}

void nfsft_export_wisdom(const char* filename)
{
  int i;
  int nul;
  int bwexp_max = BWEXP_MAX;
  nul = 0;
  FILE *f = fopen(filename,"w");
  fwrite(&bwexp_max,sizeof(int),1,f);
  for (i = 0; i <= BWEXP_MAX; i++)
  {
    if (wisdom.transform_wisdoms[i] != 0)
    {
      export_transform_wisdom(wisdom.transform_wisdoms[i],f);
    }
    else
    {
      fwrite(&nul,sizeof(int),1,f);      
    }
  }  
  fclose(f);
}

void nfsft_import_wisdom(const char* filename)
{
  FILE *f = fopen(filename,"r");
  int max;
  fread(&max,sizeof(int),1,f);
  printf("BWEXP_MAX=%d\n",max);
  fclose(f);
}*/

// --- Public interface ---

void nfsft_compute_wisdom(int M)
{
  int N;
  int t;
  t = (int) ceil(log((double)M)/log(2.0));
  /* Calculate N as next greater power of 2 of the bandwidth M. */
  N = pow2(t);  
  init_wisdom();
  if (M > 1)
  {  
    /*printf("M = %d\n",M);
    printf("r = %f\n",log((double)M)/log(2.0));
    printf("t = %d\n",t);*/
    init_transform_wisdom(N,THRESHOLD,t);
  }
}

nfsft_plan nfsft_create_plan(int type, int d, int m, double *angles, 
                             fftw_complex **f_hat, fftw_complex *f, 
                             nfsft_flags flags)
{
  nfsft_plan plan = malloc(sizeof(struct nfsft_plan_s));
  
  plan->D = d;
  plan->M = m;
  plan->angles = angles;
  plan->f_hat = f_hat;
  plan->f = f;
  plan->kind = type;
  plan->threshold = THRESHOLD;
  
  return(plan);
}

void nfsft_destroy_plan(nfsft_plan plan)
{
  free(plan);
}

void ndsft_execute(nfsft_plan plan)
{
  int t,N,n;
  struct nfsft_transform_wisdom *tw;
  
  init_wisdom();
  
  t = (int) ceil(log((double)plan->M)/log(2.0));
  /* Calculate N as next greater power of 2 of the bandwidth M. */
  N = 1<<t;
  
  if (plan->M == 0)
  {
    for (n = 0; n < plan->D; n++)
    {
      plan->f[n] = plan->f_hat[0][0];
    }  
  }
  else
  {
    /* Compute direct transform. */
    ndsft(plan->angles,plan->D,plan->f_hat, plan->f, plan->M, N, 
      tw, &wisdom);
  }  
}

void nfsft_execute(nfsft_plan plan)
{
  /** Exponent of next greater power of 2 realtive to M */
  int t = plan->M==0?0:(int) ceil(log((double)plan->M)/log(2.0));
  /** Next greater power of 2 relative to M */
  int N = plan->M==0?0:1<<t;
  /** Counter for loops*/
  int i, n;
  /** */
  struct nfsft_transform_wisdom *tw;
  /** */
  nfft_plan myplan;

  /* Init global structure. */
  init_wisdom();

  if (plan->M < 3)
  {
    ndsft_execute(plan);
  }
  else
  {
    /* Ensure that precomputation has been done. */
    if (wisdom.transform_wisdoms[t] == 0)
    {
      tw = init_transform_wisdom(N,plan->threshold,t);
    }  
    else
    {
      tw = wisdom.transform_wisdoms[t];
    }  
    
    if (plan->kind == NFSFT_BACKWARD)
    {
#ifdef DEBUG
      printf("Executing NFSFT_BACKWARD.\n");
#endif
      /* Calculate FLFT. */
      for(n = -plan->M, i = 0; n <= plan->M; n++, i++) 
      {
        flft(tw->U, plan->M, t, abs(n), plan->f_hat[i], tw);
      }
//#define DEBUG
#ifdef DEBUG
      printf("M = %d\n",plan->M);
      for (n = -plan->M; n <= plan->M; n++)
      {
        printf("-- n = %d --\n",n);
        for(i = 0; i < N+1; i++)
        {
          printf("%f + %fi\n",creal(plan->f_hat[n+plan->M][i]),cimag(plan->f_hat[n+plan->M][i]));
        }  //myvprc(myplan.f_hat,2*(N+1)*2*(N+1),"f_hat after adjoint nfft:");
        printf("----\n");
        //myvprc(plan->f_hat[n+plan->M],N+1,"f_hat after adjoint cheb2exp:");        
      }
#endif
#undef DEBUG
      
      /* Calculate NFFT. */
      nfft_size[0] = 2*(N+1);
      nfft_size[1] = 2*(N+1);
      fftw_size[0] = 4*N;
      fftw_size[1] = 4*N;
      nfft_init_specific(&myplan, 2, nfft_size, plan->D, fftw_size, 
                         6, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT | FFT_OUT_OF_PLACE, 
                         FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
      /* Assign angle array. */
      myplan.x = plan->angles;
      if (myplan.nfft_flags & PRE_PSI) 
      {  
        nfft_precompute_psi(&myplan); 
      } 
      /* Assign result array. */
      myplan.f = plan->f;
      
      /* Convert Chebyshev coefficients to Fourier coefficients. */
      cheb2exp(myplan.f_hat, plan->f_hat, plan->M, N); 
      
      /* Execute NFFT. */
      nfft_trafo(&myplan);
      nfft_finalize(&myplan);
    }  
    else if (plan->kind == NFSFT_ADJOINT)
    {
#ifdef DEBUG
      printf("Executing NFSFT_ADJOINT.\n");
#endif
      
      /* Calculate adjoint NFFT. */
      nfft_size[0] = 2*(N+1);
      nfft_size[1] = 2*(N+1);
      fftw_size[0] = 4*N;
      fftw_size[1] = 4*N;
      nfft_init_specific(&myplan, 2, nfft_size, plan->D, fftw_size, 
                         6, PRE_PHI_HUT| PRE_PSI| MALLOC_F_HAT | FFT_OUT_OF_PLACE, 
                         FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
      /* Assign angle array. */
      myplan.x = plan->angles;
      /* Assign result array. */
      myplan.f = plan->f;
      if (myplan.nfft_flags & PRE_PSI) 
      {  
        nfft_precompute_psi(&myplan); 
      } 
      
      /* Execute adjoint NFFT. */
      nfft_adjoint(&myplan);

//#define DEBUG
#ifdef DEBUG
      for(n = 0; n < 4*(N+1)*(N+1); n++)
      {
        printf("%f + %fi\n",creal(myplan.f_hat[n]),cimag(myplan.f_hat[n]));
      }  //myvprc(myplan.f_hat,2*(N+1)*2*(N+1),"f_hat after adjoint nfft:");
#endif
      
      /* Convert Chebyshev coefficients to Fourier coefficients. */
      cheb2exp_adjoint(myplan.f_hat, plan->f_hat, plan->M, N); 

//#define DEBUG
#ifdef DEBUG      
      int k;
      for (n = -plan->M; n <= plan->M; n++)
      {  
        for (k = 0; k <= N; k++)
        {
          printf("%17.16f\n%17.16f\n",creal(plan->f_hat[n+plan->M][k]),cimag(plan->f_hat[n+plan->M][k]));
        }
      }
#endif
#undef DEBUG
      
//#define DEBUG
#ifdef DEBUG
      printf("M = %d\n",plan->M);
      for (n = -plan->M; n <= plan->M; n++)
      {
        printf("-- n = %d --\n",n);
        for(i = 0; i < N+1; i++)
        {
          printf("%f + %fi\n",creal(plan->f_hat[n+plan->M][i]),cimag(plan->f_hat[n+plan->M][i]));
        }  //myvprc(myplan.f_hat,2*(N+1)*2*(N+1),"f_hat after adjoint nfft:");
        printf("----\n");
        //myvprc(plan->f_hat[n+plan->M],N+1,"f_hat after adjoint cheb2exp:");        
      }
#endif
#undef DEBUG
      
      /* Calculate FLFT. */
      for(n = -plan->M, i = 0; n <= plan->M; n++, i++) 
      {
        flft_adjoint(tw->U,plan->M,t,abs(n),plan->f_hat[i],tw);
      }    
      
      nfft_finalize(&myplan);
    }
    else if (plan->kind == NFSFT_FORWARD)
    {
    }
    else
    {
#ifdef DEBUG
      printf("Wrong transform type.\n");
#endif
    }
  }    
}

double norm(complex* x, int length)
{
  double r = 0.0;
  int i;
  for (i = 0; i < length; i++)
  {
    r += cabs(x[i]);
  }  
  return r;
}

void nfsft_test()
{
  const int N = 128;
  int t = (int)ceil(log((double)N)/log(2.0));
  int j;
  
  double m1[4] = { 2.0             ,  2.0             ,  2.0             ,  2.0             };
  double m2[4] = { 3.69551813004515,  1.53073372946036, -1.53073372946036, -3.69551813004515};
  double m3[4] = {-0.92387953251129, -0.38268343236509,  0.38268343236509,  0.92387953251129};
  double m4[4] = { 1.70710678118655,  0.29289321881345,  0.29289321881345,  1.70710678118655};
  complex x[4] = {1.0, 1.0, 0.0, 0.0};
  complex y[4] = {3.0,-2.0, 0.0, 0.0};
  struct U_type u;
  int tau = 1;
  int l = 1;
  int n = 0;
  struct nfsft_transform_wisdom *tw;
  init_wisdom();
  tw = init_transform_wisdom(16,1000000,4);
  u.m1 = m1;
  u.m2 = m2;
  u.m3 = m3;
  u.m4 = m4;
  //multiplyU(x,y,u,tau,n,l,tw,1.3);
  multiplyU_adjoint(x,y,u,tau,n,l,tw,1.3);
  myvprc(x,4,"x");
  myvprc(y,4,"y");
  /*int i,j;
  complex *f_hat_orig;
  complex *f, *f_hat_iter;
  complex *r, *z, *p, *v;
  double alpha, beta, old_norm;
  
  f_hat_orig = malloc((N+1)*sizeof(complex));
  f = malloc((N+1)*sizeof(complex));
  f_hat_i = malloc((N+1)*sizeof(complex));
  r = malloc((N+1)*sizeof(complex));
  z = malloc((N+1)*sizeof(complex));
  p = malloc((N+1)*sizeof(complex));
  v = malloc((N+1)*sizeof(complex));
  for (i = 0; i < (N+1); i++)
  {
    f_hat_orig[i] = (double)(i+1);
    f_hat_i[i] = 0.0;
  }*/  
    
  /*printf("------------------------\nN = %d, t = %d, n = %d\n-------------------\n",N,t,n);
  
  memcpy(f,f_hat_orig,(N+1)*sizeof(complex));  
  myvprc(f,(N+1),"f");
  printf("Applying FLFT:\n");
  flft(tw->U,N,t,abs(n),f,tw);
  myvprc(f,(N+1),"f");
  printf("Applying adjoint FLFT:\n");
  flft_adjoint(tw->U,N,t,abs(n),f,tw);
  myvprc(f,(N+1),"f");
  return;
  
  memcpy(f,f_hat_orig,(N+1)*sizeof(complex));
  printf("Before FLFT:\n");
  myvprc(f_hat_orig,(N+1),"f_hat_orig");
  myvprc(f,(N+1),"f");
  printf("\n");*/
  
  /* Calculate FLFT. */
  /*flft(tw->U,N,t,abs(n),f,tw);

  printf("After FLFT:\n");
  myvprc(f_hat_orig,(N+1),"f_hat_orig");
  myvprc(f,(N+1),"f");
  printf("\n");*/
  
  /* Landweber */
  /*printf("Initialization:");
  memcpy(r,f,(N+1)*sizeof(complex));
  memcpy(z,r,(N+1)*sizeof(complex));  
  flft_adjoint(tw->U,N,t,abs(n),z,tw);
  myvprc(f_hat_i,(N+1),"f_hat_i");
  myvprc(r,(N+1),"r");
  myvprc(z,(N+1),"z");
  return;
/*  printf("\n");
  for (i = 1; i < 200; i++)
  {
    printf("Iteration %d:\n",i);
    for (j = 0; j < (N+1); j++)
    {
      f_hat_i[j] += z[j];
    }
    myvprc(f_hat_i,(N+1),"f_hat_i");
    memcpy(r,f_hat_i,(N+1)*sizeof(complex));
    flft(tw->U,N,t,abs(n),r,tw);
    for (j = 0; j < (N+1); j++)
    {
      r[j] = f[j] - r[j];
    }
    myvprc(r,(N+1),"r");
    memcpy(z,r,(N+1)*sizeof(complex));
    flft_adjoint(tw->U,N,t,abs(n),z,tw);
    myvprc(z,(N+1),"z");
    printf("\n");
  }*/
  
  // CGNE
  /*printf("Initialization:");
  memcpy(r,f,(N+1)*sizeof(complex));
  memcpy(z,r,(N+1)*sizeof(complex));  
  flft_adjoint(tw->U,N,t,abs(n),z,tw);
  memcpy(p,z,(N+1)*sizeof(complex));  
  myvprc(f_hat_i,(N+1),"f_hat_i");
  myvprc(r,(N+1),"r");
  myvprc(z,(N+1),"z");
  myvprc(p,(N+1),"p");
  myvprc(v,(N+1),"v");

  for (i = 1; i < 200; i++)
  {*/
    // v = AWz
    //memcpy(v,p,(N+1)*sizeof(complex));  
    //flft(tw->U,N,t,abs(n),v,tw);    
    // \alpha_l = ...
    //alpha = norm(z,N+1)/norm(v,N+1);
    // f_hat = f_Hat + \alpha*p
    // r = r - alpha*v;
    /*for (j = 0; j < N+1; j++)
    {
      f_hat_i[j] += alpha*p[j];
      r[j] -= alpha*v[j]; 
    } */ 
    // z = A^H W r
    /*old_norm = norm(z,N+1);
    memcpy(z,r,(N+1)*sizeof(complex));  
    flft_adjoint(tw->U,N,t,abs(n),z,tw);*/
    // \beta = ...
    //beta = norm(z,N+1)/old_norm;
    // p = \beta*p + z
    /*for (j = 0; j < N+1; j++)
    {
      p[j] = beta*p[j] + z[j]; 
    }  
  }*/

  /*myvprc(f_hat_orig,(N+1),"f_hat_orig");
  myvprc(f_hat_i,(N+1),"f_hat_i");
  myvprc(f,(N+1),"f");
  myvprc(r,(N+1),"r");
  
  free(f);
  free(f_hat_orig);
  free(f_hat_i);
  free(r);
  free(z);
  free(p);
  free(v);*/
  /*double x[5] = {-1.0,0.5,0.0,0.5,1.0};
  double y[5];
  int n = 1;
  for (j = n; j < 10; j++)
  {
    eval_al(x,y,5,j,&(wisdom.alpha[ROW(n)]),&(wisdom.beta[ROW(n)]),&(wisdom.gamma[ROW(n)]));
    printf("P_%d^%d(1) =",j,n);
    myvpr(y,5,"");
  }*/  
}
