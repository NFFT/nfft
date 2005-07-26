#include "direct.h"


void ndsft(int D, double *angles, complex *f, int M, complex **f_hat, 
           struct nfsft_wisdom *wisdom)
{
  /** Node index */
  int d;
	 /** Legendre index k */
  int l;
	 /** Legendre index n */
	 int n;
	 /** nleg = |n| */
	 int nleg;
  
	 /**
   * Pointer to array containing three-term recurrence coefficients \alpha for 
	  * associated Legendre functions. 
	  */
  double *alpha;
	 /**
   * Pointer to array containing three-term recurrence coefficients \gamma for 
	  * associated Legendre functions. 
	  */
	 double *gamma;
  
  /** Index used in Clenshaw algorithm. */
  int index;
 	/** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *a;
 	/** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *temp;
  
  /** The final result for a single node. */
  complex f_d;
  
  /** Used to store the angles theta_d */
  double *theta;
  /** Used to store the angles phi_d */
  double *phi;
    
	 /* Allocate memory for auxilliary arrays. */
  temp = (complex*) malloc (sizeof(complex) * (M+1));
  theta = (double*) malloc (sizeof(double) * D);    
  phi = (double*) malloc (sizeof(double) * D);    
  
  /* Scale angles phi_j from [-1,1] to [-pi,pi] and angles \theta_j from [0,1] 
   * to and [0,pi], respectively. */ 
  for (d = 0; d < D; d++)
  {
    theta[d] = 2.0*PI*angles[2*d+1];
    phi[d] = -2.0*PI*angles[2*d];
  }
  
	 /* Distinguish by bandwidth M. */
  if (M == 0)
  {
    /* Constant function */
    for (d = 0; d < D; d++)
   	{
      f[d] = f_hat[0][0];
    }
  }
  else
  { 
    /* Apply cosine to angles theta_j. */
    for (d = 0; d < D; d++)
    {
	     theta[d] = cos(theta[d]);
    }
    
    /* Evaluate 
     *      \sum_{k=0}^M \sum_{n=-k}^k a_k^n P_k^{|n|}(\cos\theta_j) e^{i n \phi_j} 
		   *    = \sum_{n=-M}^M \sum_{k=|n|}^M a_k^n P_k^{|n|}(\cos\theta_j) e^{i n \phi_j}. */
    for (d = 0; d < D; d++)
	   {
			   /* Initialize result for current node. */
      f_d = 0.0;

      /* For n = -M,...,M, evaluate 
			    *   b_n := \sum_{k=|n|}^M a_k^n P_k^{|n|}(\cos\vtheta_j)
			    * using Clenshaw's algorithm. */
      for (n = -M; n <= M; n++)
      {
			     /* Get Fourier coefficients vector. */
        a = f_hat[n+M];
        
				    /* Take absolute value of n. */
        nleg = abs(n);
        
        /* Get three-term recurrence coefficients vectors. */
        alpha = &(wisdom->alpha[ROWK(nleg)]);
        gamma = &(wisdom->gamma[ROWK(nleg)]);
	       
        /* Make copy of array a. */ 
        memcpy(temp,a,(M+1)*sizeof(complex));
        
        /* Clenshaw's algorithm */        
        for (l = M; l > nleg + 1; l--)
        {
          index = l - nleg;
          temp[l-1] += temp[l] * alpha[index] * theta[d]; 
          temp[l-2] += temp[l] * gamma[index];
        }
        
        /* Compute final step if neccesary. */
        if (nleg < M)
        {  
          temp[0+nleg] += temp[1+nleg] * wisdom->alpha[ROWK(nleg)+1] * theta[d];
        }
        
        /* Write final result b_n of multiplication by normalization constant to 
         * array b  = (b_{-M},...,b_M). */
        f_d += temp[0+nleg] * wisdom->gamma[ROW(nleg)] *
          pow(1- theta[d] * theta[d], 0.5*nleg) * cexp(I*n*phi[d]);
	     }
            
			   /* Write result to vector f. */
      f[d] = f_d;
    }
  }  
  
 	/* Free auxilliary arrays. */
  free(phi);
  free(theta);
  free(temp);
}


void adjoint_ndsft(int D, double *angles, complex *f, int M, complex **f_hat,
                   struct nfsft_wisdom *wisdom)
{
  /** Node index */
  int d;
	 /** Legendre index k */
  int l;
	 /** Legendre index n */
	 int n;
	 /** nleg = |n| */
	 int nleg;
  
	 /**
   * Pointer to array containing three-term recurrence coefficients \alpha for 
	  * associated Legendre functions. 
	  */
  double *alpha;
	 /**
    * Pointer to array containing three-term recurrence coefficients \gamma for 
	  * associated Legendre functions. 
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
  temp = (complex*) malloc (sizeof(complex) * (M+1));
  theta = (double*) malloc (sizeof(double) * D);    
  phi = (double*) malloc (sizeof(double) * D);    
  
  /* Scale angles phi_j from [-1,1] to [-pi,pi] and angles \theta_j from [0,1] 
   * to and [0,pi], respectively. */ 
  for (d = 0; d < D; d++)
  {
    theta[d] = 2.0*PI*angles[2*d+1];
    phi[d] = -2.0*PI*angles[2*d];
  }

  for (n = -M; n <= M; n++)
  {
    for (l = abs(n); l <= M; l++)
    {
      f_hat[n+M][l] = 0.0;
    }  
  }
  
	 /* Distinguish by bandwidth M. */
  if (M == 0)
  {
    /* Constant function */
    for (d = 0; d < D; d++)
   	{
      f_hat[0][0] += f[d];
    }
  }
  else
  { 
    /* Apply cosine to angles theta_j. */
    for (d = 0; d < D; d++)
    {
	     theta[d] = cos(theta[d]);
    }
    
    for (d = 0; d < D; d++)
	   {
      for (n = -M; n <= M; n++)
      {
				    /* Take absolute value of n. */
        nleg = abs(n);

        /* Get three-term recurrence coefficients vectors. */
        alpha = &(wisdom->alpha[ROWK(nleg)]);
        gamma = &(wisdom->gamma[ROWK(nleg)]);

        temp[nleg] = f[d] * cexp(-I*n*phi[d]) * wisdom->gamma[ROW(nleg)]/*gamma_m1[nleg]*/ *
          pow(1- theta[d] * theta[d], 0.5*nleg);
                
        /* Compute final step if neccesary. */
        if (nleg < M)
        {  
          temp[nleg+1] = temp[nleg] * wisdom->alpha[ROWK(nleg)+1] * theta[d];
        }

        /* Clenshaw's algorithm */        
        for (l = nleg+2; l <= M; l++)
        {
          index = l - nleg;
          temp[l] = alpha[index] * theta[d] * temp[l-1] + gamma[index] * temp[l-2];
        }
        
        /* Copy result */
        for (l = nleg; l <= M; l++)
        {
          f_hat[n+M][l] += temp[l];
        }  
	     }      
    }
  }  
  
 	/* Free auxilliary arrays. */
  free(phi);
  free(theta);
  free(temp);
}
