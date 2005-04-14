#include "direct.h"

void direct_trafo(double *angles, int D, complex **f_hat, complex *f, 
  int M, int N, struct nfsft_transform_wisdom *tw, struct nfsft_wisdom *wisdom)
{
  int j, nleg,n;
  int index;
  double *alpha_pre_n,*gamma_pre_n;
  complex *a;
  complex result_j;
  double sin_n_phi,cos_n_phi;
  int l,N1;
  complex *h;
  complex *b;  
  int d;
  
  N1 = N + 1;
  
  h = (complex*) fftw_malloc (sizeof(complex) * (N1));
  b = (complex*) fftw_malloc (sizeof(complex) * (2*M+1));    
  
  /* Convert angles. */ 
  for(d = 0; d < D; d++)
  {
    angles[2*d+1] = 2.0*PI*angles[2*d+1];
    angles[2*d] = -2.0*PI*angles[2*d];
  }
  
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
    /* Convert theta angles. */
    for (j = 0; j < D; j++)
	   {
	     angles[2*j+1] = cos(angles[2*j+1]);
	   }
    
    /* Evaluate at every node. */
    for (j = 0; j < D; j++)
	   {
	     /* Clenshaw */
	     for (n = -M; n <= M; n++)
	     {
        a = f_hat[n+M];
        
	       nleg = abs(n);
        
	       alpha_pre_n = &(wisdom->alpha[ROWK(nleg)]);
	       gamma_pre_n = &(wisdom->gamma[ROWK(nleg)]);
	       
        memcpy(h,a,N1*sizeof(complex));
        
	       /* Clenshaw */
	       for (l = M; l > nleg + 1; l--)
		      {
		        index = l - nleg;
          h[l-1] += h[l]*alpha_pre_n[index]*angles[2*j+1]; 
          h[l-2] += h[l]*gamma_pre_n[index];
        }
        
        if (nleg < M)
        {  
 	        h[0+nleg] += h[1+nleg] * wisdom->alpha[ROWK(nleg)+1] * angles[2*j+1];          
        }
        
	       b[n+M] = h[0+nleg]*wisdom->gamma_m1[nleg]*pow(1-angles[2*j+1]*angles[2*j+1],nleg/2.0);
	     }
      
      /* fourier */
      
      result_j = 0.0;
      
      for (n = -M; n <= M; n++)
      {
        sin_n_phi=sin(n*angles[2*j]); 
        cos_n_phi=cos(n*angles[2*j]);
        
        result_j += (creal(b[n+M])*cos_n_phi-cimag(b[n+M])*sin_n_phi) + I*(creal(b[n+M])*sin_n_phi+cimag(b[n+M])*cos_n_phi);
      }
      
      f[j] = result_j;
    }
  }  
  fftw_free(b);  
  fftw_free(h);
}

void direct_trafo_adjoint(double *angles, int D, complex **f_hat, complex *f, 
                  int M, int N, struct nfsft_transform_wisdom *tw, struct nfsft_wisdom *wisdom)
{
  int j, nleg,n,k;
  int index;
  double *alpha_pre_n,*gamma_pre_n;
  complex *a;
  complex result_j;
  double sin_n_phi,cos_n_phi;
  int l,N1;
  complex *h;
  complex *b;  
  int d;
  double *theta;
  
  N1 = N + 1;
  
  theta = (complex*) fftw_malloc (sizeof(complex) * (D));
  a = (complex*) fftw_malloc (sizeof(complex) * (D));    
  b = (complex*) fftw_malloc (sizeof(complex) * (D));    
  
  /* Convert angles. */ 
  for(d = 0; d < D; d++)
  {
    angles[2*d+1] = 2.0*PI*angles[2*d+1];
    angles[2*d] = -2.0*PI*angles[2*d];
  }
  
  if (M == 0)
  {
    /* Constant function. */
    for (j = 0; j < D; j++)
   	{
      f_hat[0][0] = f[j];
	   }
  }
  else
  { 
    /* Convert theta angles. */
    for (j = 0; j < D; j++)
	   {
	     angles[2*j+1] = cos(angles[2*j+1]);
      theta[d] = angles[2*d+1];
	   }
    
    for (k = 0; k <= M; k++)
    {
      for (n = -k; n <= k; n++)
      {
        /* Clenshaw to evaluate P_k^{|n|}(cos(\theta)) */
        for (j = 0; j < D; j++)
        {
          a[j] = 1.0;
          b[j] = 0.0;
        }  
      }
    }
    
    /* Evaluate at every node. */
    for (j = 0; j < D; j++)
	   {
	     /* Clenshaw */
	     for (n = -M; n <= M; n++)
	     {
        a = f_hat[n+M];
        
	       nleg = abs(n);
        
	       alpha_pre_n = &(wisdom->alpha[ROWK(nleg)]);
	       gamma_pre_n = &(wisdom->gamma[ROWK(nleg)]);
	       
        memcpy(h,a,N1*sizeof(complex));
        
	       /* Clenshaw */
	       for (l = M; l > nleg + 1; l--)
		      {
		        index = l - nleg;
          h[l-1] += h[l]*alpha_pre_n[index]*angles[2*j+1]; 
          h[l-2] += h[l]*gamma_pre_n[index];
        }
        
        if (nleg < M)
        {  
 	        h[0+nleg] += h[1+nleg] * wisdom->alpha[ROWK(nleg)+1] * angles[2*j+1];          
        }
        
	       b[n+M] = h[0+nleg]*wisdom->gamma_m1[nleg]*pow(1-angles[2*j+1]*angles[2*j+1],nleg/2.0);
	     }
      
      /* fourier */
      
      result_j = 0.0;
      
      for (n = -M; n <= M; n++)
      {
        sin_n_phi=sin(n*angles[2*j]); 
        cos_n_phi=cos(n*angles[2*j]);
        
        result_j += (creal(b[n+M])*cos_n_phi-cimag(b[n+M])*sin_n_phi) + I*(creal(b[n+M])*sin_n_phi+cimag(b[n+M])*cos_n_phi);
      }
      
      f[j] = result_j;
    }
  }  
  fftw_free(a);  
  fftw_free(b);  
  fftw_free(theta);
}  