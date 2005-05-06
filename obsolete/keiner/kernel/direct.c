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
 	/** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *b;  

  /** The final result for a single node. */
  complex f_d;
  
  /** Used to store the angles theta_d */
  double *theta;
  /** Used to store the angles phi_d */
  double *phi;
  
  
	 /* Allocate memory for auxilliary arrays. */
  temp = (complex*) malloc (sizeof(complex) * (M+1));
  b = (complex*) malloc (sizeof(complex) * (2*M+1));    
  theta = (double*) malloc (sizeof(double) * D);    
  phi = (double*) malloc (sizeof(double) * D);    
  
  /* 
 	 * Scale angles phi_j from [-1,1] to [-pi,pi] and angles \theta_j from [0,1] 
   * to and [0,pi], respectively. 
 	 */ 
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
    
    /* 
  		 * Evaluate 
     *      \sum_{k=0}^M \sum_{n=-k}^k a_k^n P_k^{|n|}(\cos\theta_j) e^{i n \phi_j} 
		   *    = \sum_{n=-M}^M \sum_{k=|n|}^M a_k^n P_k^{|n|}(\cos\theta_j) e^{i n \phi_j}.
		   */
    for (d = 0; d < D; d++)
	   {
      /* 
			    * For n = -M,...,M, evaluate 
			    *   b_n := \sum_{k=|n|}^M a_k^n P_k^{|n|}(\cos\theta_j)
			    * using Clenshaw's algorithm. 
			    */
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
				    * array b  = (b_{-M},...,b_M). 
				    */
       b[n+M] = temp[0+nleg] * wisdom->gamma_m1[nleg] * pow(1- theta[d] * 
				      theta[d], nleg/2.0);
	     }
      
			   /* Initialize result for current node. */
      f_d = 0.0;
      
			   /* Compute
			    *   f_j := \sum_{n=-M}^M b_n e^{i n \phi_j}.
			    */
      for (n = -M; n <= M; n++)
      {        
        f_d += b[n+M]*cexp(I*n*phi[d]);
      }
      
			   /* Write result to vector f. */
      f[d] = f_d;
    }
  }  
	
 	/* Free auxilliary arrays. */
  free(phi);
  free(theta);
  free(b);  
  free(temp);
}


void adjoint_ndsft(int D, double *angles, complex *f, int M, complex **f_hat, 
  struct nfsft_wisdom *wisdom)
{
  /** Legendre index k */
 	int k;
 	/** Legendre index n */
 	int n;

  /** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *a;
  /** Pointer to auxilliary array for Clenshaw algorithm. */
 	complex *temp;
  /** Pointer to auxilliary array for Clenshaw algorithm. */
  complex *b;  

  /** Used to store the angles phi_d */ 
  double *phi;
  /** Used to store the angles theta_d */
  double *theta;
  /** Used to store the values sin(theta_d) */
 	double *sint;
  /** Used to store the values cos(theta_d) */
 	double *cost;
 
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

  int j, l, nleg;
 	complex *y;
 	complex result;
 	double gamma_nleg;
	
  a = (complex*) malloc (sizeof(complex) * D);    
  b = (complex*) malloc (sizeof(complex) * D);
  temp = (complex*) malloc (sizeof(complex) * D);
  y = (complex*) malloc (sizeof(complex) * D);
  phi = (double*) malloc (sizeof(double) * D);
  sint = (double*) malloc (sizeof(double) * D);
  cost = (double*) malloc (sizeof(double) * D);
	
  /* 
 	 * Convert angles from [-1,1] and [0,1] to [-pi,pi] and [0,pi] 
	  * respectively. 
	  */ 
  for(j = 0; j < D; j++)
  {
    phi[j] = - 2.0 * PI * angles[2*j];
    sint[j] = - 2.0 * PI * angles[2*j+1];
    cost[j] = - 2.0 * PI * angles[2*j+1];
    sint[j] = sin(sint[j]);
    cost[j] = cos(cost[j]);
  }

	 for (k = 0; k <= M; k++)
	 {
	   for (n = -k; n <= k; n++)
		  {
		    /* Evaluate
		 	   *   \sum_{j=0}^{D-1} f_j P_k^{|n|}(\cos\theta_j) e^{-i n \phi_j}
			    */
			
			   /* Evaluate b_j := (f_j P_k^{|n|}(\cos\theta_j))_{j=0}^{D-1}. */
						
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
	       alpha = &(wisdom->alpha[ROW(nleg)]);
	       gamma = &(wisdom->gamma[ROW(nleg)]);
	       
        /* Clenshaw's algorithm */
	  		
 				   /* Initialize vector a. */
   			  for (j = 0; j < D; j++)
		   	  {
			       a[j] = 1.0;
				      b[j] = 0.0;
			     }  

        for (l = k; l >= nleg + 1; l--)
		      {
  				    /* Make copy of array a. */ 
          memcpy(temp, a, D * sizeof(complex));
  			     for (j = 0; j < D; j++)
	  	  	   {
            a[j] = b[j] + temp[j] * (alpha[l] * cost[j]);		        
            b[j] = temp[j] * gamma[l];
 			      }
        }	
				
				    gamma_nleg = wisdom->gamma_m1[nleg];
        
        for (j = 0; j < D; j++)
        {
          y[j] = a[j];
        }			

        /*for (j = 0; j < D; j++)
  	  	  {
          y[j] *= f[j] * gamma_nleg * pow(-1.0,nleg) * pow(1.0 - cost[j] * cost[j], nleg/2.0);
	       }*/
        
        for (l = 0; l < nleg; l++)
        {  
          for (j = 0; j < D; j++)
    	  	  {
            y[j] *= (-1.0)*sint[j];
  	       }						  
        }
        for (j = 0; j < D; j++)
        {
          y[j] *= gamma_nleg * f[j];
        }						  
			   }
			
		    /* Evaluate
			    *   \sum_{j=0}^{D-1} b_j e^{-i n \phi_j}
			    */
			   result = 0.0;
      for (j = 0; j < D; j++)
  	   {
        result += y[j] * cexp(-I*n*phi[j]);
			   }						  
			   f_hat[n+M][k] = result;
		  }
	 }

  /* Free auxilliary data arrays. */
  free(cost);
 	free(sint);
  free(phi);
  free(y);  
  free(temp);
  free(b);  
  free(a);  
}  
