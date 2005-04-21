#include "direct.h"

void ndsft(int D, double *angles, complex *f, int M, int N, 
  complex **f_hat, struct nfsft_transform_wisdom *tw, 
	struct nfsft_wisdom *wisdom)
{
  /** Node index */
  int j;
	/** Fourier-index k */
  int l;
	/** Fourier-index n */
	int n;
	/** nleg = |n| */
	int nleg;

	/**
	 * Pointer to array containing three-term recurrence coefficients \alpha for 
	 * associated Legendre functions. 
	 */
  double *alpha_pre_n;
	/**
	 * Pointer to array containing three-term recurrence coefficients \gamma for 
	 * associated Legendre functions. 
	 */
	double *gamma_pre_n;

  /** Index used in Clenshaw algorithm. */
  int index;
	/** Pointer to auxilliary array used in Clenshaw algorithm. */
  complex *a;

	/** Auxilliary arrays. */
  complex *h;
  complex *b;  

  complex result_j;
  double sin_n_phi,cos_n_phi;

  
	/* Allocate memory for auxilliary arrays. */
  h = (complex*) malloc (sizeof(complex) * (N+1));
  b = (complex*) malloc (sizeof(complex) * (2*M+1));    
  
  /* 
	 * Convert angles from [-1,1] and [0,1] to [-pi,pi] and [0,pi] 
	 * respectively. 
	 */ 
  for(j = 0; j < D; j++)
  {
    angles[2*j+1] = 2.0*PI*angles[2*j+1];
    angles[2*j] = -2.0*PI*angles[2*j];
  }
  
	/* Distinguish by bandwidth M. */
  if (M == 0)
  {
    /* Constant function. */
    for (j = 0; j < D; j++)
   	{
      f[j] = f_hat[0][0];
	  }
  }
  else
  { 
    /* Apply cosine function to angles theta. */
    for (j = 0; j < D; j++)
    {
	    angles[2*j+1] = cos(angles[2*j+1]);
    }
    
    /* 
		 * Evaluate 
     *    \sum_{k=0}^M \sum_{n=-k}^k a_k^n P_k^{|n|}(\cos\theta_j) e^{i n \phi_j} 
		 *  = \sum_{n=-M}^M \sum_{k=|n|}^M a_k^n P_k^{|n|}(\cos\theta_j) e^{i n \phi_j}.
		 */
    for (j = 0; j < D; j++)
	  {
      /* 
			 * For n = -M,...,M, evaluate 
			 *   b_n := \sum_{k=|n|}^M a_k^n P_k^{|n|}(\cos\theta_j)
			 * by Clenshaw algorithm. 
			 */
      for (n = -M; n <= M; n++)
      {
			  /* Get corresponding Fourier-coefficients vector. */
        a = f_hat[n+M];
				
				/* Take absolute value of n. */
        nleg = abs(n);
        
				/* Get corresponding three-term recurrence coefficients vectors. */
	      alpha_pre_n = &(wisdom->alpha[ROWK(nleg)]);
	      gamma_pre_n = &(wisdom->gamma[ROWK(nleg)]);
	       
				/* Make copy of array a. */ 
        memcpy(h,a,(N+1)*sizeof(complex));

        /* Clenshaw algorithm */        
        for (l = M; l > nleg + 1; l--)
		    {
		      index = l - nleg;
          h[l-1] += h[l] * alpha_pre_n[index] * angles[2*j+1]; 
          h[l-2] += h[l] * gamma_pre_n[index];
        }
        
				/* Compute final step if neccesary. */
        if (nleg < M)
        {  
 	        h[0+nleg] += h[1+nleg] * wisdom->alpha[ROWK(nleg)+1] * angles[2*j+1];          
        }
        
				/* 
				 * Write final result b_n of multiplication by normalization constant to 
				 * array b  = (b_{-M},...,b_M). 
				 */
	       b[n+M] = h[0+nleg] * wisdom->gamma_m1[nleg] * pow(1 - angles[2*j+1] * 
				   angles[2*j+1], nleg/2.0);
	     }
      
			/* Initialize result for current node. */
      result_j = 0.0;
      
			/* Compute
			 *   f_j := \sum_{n=-M}^M b_n e^{i n \phi_j}.
			 */
      for (n = -M; n <= M; n++)
      {
        sin_n_phi=sin(n*angles[2*j]); 
        cos_n_phi=cos(n*angles[2*j]);
        
				/** \todo Use complex arithmetic here. */
        result_j += (creal(b[n+M])*cos_n_phi-cimag(b[n+M])*sin_n_phi) + 
				  I * (creal(b[n+M])*sin_n_phi+cimag(b[n+M])*cos_n_phi);
      }
      
			/* Write result to vector f. */
      f[j] = result_j;
    }
  }  
	
	/* Free auxilliary arrays. */
  free(b);  
  free(h);
}


void adjoint_ndsft(int D, double *angles, complex *f, int M, int N, 
  complex **f_hat, struct nfsft_transform_wisdom *tw, 
	struct nfsft_wisdom *wisdom)
{
  /** Fourier-index k */
	int k;
	/** Fourier-index n */
	int n;
  int j, l, nleg;
  //int index;
  double *alpha_pre_n,*gamma_pre_n;
  //complex result_j;
  //double sin_n_phi,cos_n_phi;
  //int l,N1;
  //complex *h;
  //complex *b;  
  //int d;
  //double *theta;
  
  complex *a;
	complex *b;
	complex *temp;
	complex *y;
	complex result;
	
  a = (complex*) malloc (sizeof(complex) * D);    
  b = (complex*) malloc (sizeof(complex) * D);
  temp = (complex*) malloc (sizeof(complex) * D);
  y = (complex*) malloc (sizeof(complex) * D);
	
  /* 
	 * Convert angles from [-1,1] and [0,1] to [-pi,pi] and [0,pi] 
	 * respectively. 
	 */ 
  for(j = 0; j < D; j++)
  {
    angles[2*j+1] = 2.0 * PI * angles[2*j+1];
    angles[2*j+1] = cos(angles[2*j+1]);
    angles[2*j] = -2.0 * PI * angles[2*j];
  }

	for (k = 0; k < M; k++)
	{
	  for (n = -k; n < k; n++)
		{
		  /* Evaluate
			 *   \sum_{j=0}^{D-1} f_j P_k^{|n|}(\cos\theta_j) e^{-i n \phi_j}
			 */
			
			/* Evaluate b_j := (f_j P_k^{|n|}(\cos\theta_j))_{j=0}^D. */
						
			if (k == 0)
			{
			  for (j = 0; j < D; j++)
		  	{
			    y[j] = f[j];
			  }
			}
			else
			{
				/* Take absolute value of n. */
        nleg = abs(n);
        
				/* Get corresponding three-term recurrence coefficients vectors. */
	      alpha_pre_n = &(wisdom->alpha[ROWK(nleg)]);
	      gamma_pre_n = &(wisdom->gamma[ROWK(nleg)]);
	       
        /* Clenshaw algorithm */        
	  		
				/* Initialize vector a. */
  			for (j = 0; j < D; j++)
		  	{
			    a[j] = f[j];
				  b[j] = 0.0;
			  }  

        for (l = k; l > 1; l--)
		    {
  				/* Make copy of array a. */ 
          memcpy(temp, a, D*sizeof(complex));
  			  for (j = 0; j < D; j++)
	  	  	{
            a[j] = b[j] + temp[j]*(alpha_pre_n[l]*angles[2*j+1]);		        
            b[j] = temp[j]*gamma_pre_n[l];
			    }
        }	
			  for (j = 0; j < D; j++)
  	  	{
          y[j] = a[j]*(alpha_pre_n[0]*angles[2*j+1]) + b[j];                  
  	    }						  
			}
			
		  /* Evaluate
			 *   \sum_{j=0}^{D-1} b_j e^{-i n \phi_j}
			 */
			result = 0.0;
      for (j = 0; j < D; j++)
  	  {
        result += y[j] * cexp(-I*n*angles[2*j]);
			}						  
			f_hat[n+M][k] = result;
		}
	}

  /* Free auxilliary data arrays. */
  free(y);  
  free(temp);
  free(b);  
  free(a);  
}  
