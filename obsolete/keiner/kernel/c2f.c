#include "c2f.h"

inline void cheb2exp(complex *f_hat, complex **cheb, int M, int N)
{
  int k;             /** Degree \f$k\f$                           */
  int n;             /** Order \f$k\f$                            */
  complex *cheb_n;   /** Chebyshev coefficients for fixed \f$n\f$ */
  complex last;      /** Stores temporary values                  */
  complex act;       /** Stores temporary values                  */
  complex *f_hat_p;  /** Auxilliary pointer                       */
  complex *f_hat_n;  /** Auxilliary pointer                       */
  int l1,l2,u1,u2;
  int dim, dimh;
  int rowz;
  int colz;
  dimh = N + 1;
  dim = dimh<<1;
  rowz = (dimh)*dim;
  colz = dimh;
  
  /* Determine lower and upper bounds for loops processing even and odd terms. */
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
  
  /* initialize target with zeros */
  memset(f_hat,0U,dim*dim*sizeof(complex));
  
  /* Process even terms. */
  for (n = l1; n <= u1; n += 2)
  {
    /* Get array for current n. */
    cheb_n = cheb[n+M];
    
    f_hat[rowz+n*dim+colz] = *(cheb_n++);
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *f_hat_p++ = 0.5*(*(cheb_n));
      *f_hat_n-- = 0.5*(*(cheb_n++));
    }
  }
  
  /* Process odd terms. */
  for (n = l2; n <= u2; n += 2)
  {
    cheb_n = cheb[n+M];
    
    f_hat[rowz+n*dim+colz] = *(cheb_n++);
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *f_hat_p++ = 0.5*(*(cheb_n));
      *f_hat_n-- = 0.5*(*(cheb_n++));
    }
    
    /* Incorporate sine term. */
    last = f_hat[rowz+n*dim+colz-N];
    f_hat[rowz+n*dim+colz-N] = 0.5 * I * f_hat[rowz+n*dim+colz-N+1];
    for (k = -N+1; k <= N-1; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = 0.5 * I * (f_hat[rowz+n*dim+colz+k+1] - last);
      last = act;
    }
    f_hat[rowz+n*dim+colz+N] = - 0.5 * I * last;
  }  
}

inline void cheb2exp_adjoint(complex *f_hat, complex **cheb, int M, int N)
{
  int k;             /** Degree \f$k\f$                           */
  int n;             /** Order \f$k\f$                            */
  complex *cheb_n;   /** Chebyshev coefficients for fixed \f$n\f$ */
  complex last;      /** Stores temporary values                  */
  complex act;       /** Stores temporary values                  */
  complex *f_hat_p;  /** Auxilliary pointer                       */
  complex *f_hat_n;  /** Auxilliary pointer                       */
  int l1,l2,u1,u2;
  int dim, dimh;
  int rowz;
  int colz;
  
  dimh = N + 1;
  dim = dimh<<1;
  rowz = (dimh)*dim;
  colz = dimh;
  
  /* Determine lower and upper bounds for loops processing even and odd terms. */
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
    /* Get array for current n. */
    cheb_n = cheb[n+M];
    
    memset(cheb_n, 0U, (N+1)*sizeof(complex));
    
    *cheb_n++ = f_hat[rowz+n*dim+colz];
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *cheb_n++ = 0.5 * (*f_hat_p++ + *f_hat_n--);
    }
  }
  
  /* Process odd terms. */
  for (n = l2; n <= u2; n += 2)
  {       
    /* Incorporate sine term. */
    last = f_hat[rowz+n*dim+colz-N];
    f_hat[rowz+n*dim+colz-N] = 0.5 * I * f_hat[rowz+n*dim+colz-N+1];
    for (k = -N+1; k <= N-1; k++)
    {
      act = f_hat[rowz+n*dim+colz+k];
      f_hat[rowz+n*dim+colz+k] = 0.5 * I * (f_hat[rowz+n*dim+colz+k+1] - last);
      last = act;
    }
    f_hat[rowz+n*dim+colz+N] = -0.5 * I * last;
    
    cheb_n = cheb[n+M];
    memset(cheb_n, 0U, (N+1)*sizeof(complex));
    
    *cheb_n++ = f_hat[rowz+n*dim+colz];
    
    f_hat_p = &f_hat[rowz+n*dim+colz+1];
    f_hat_n = &f_hat[rowz+n*dim+colz-1];
    
    for(k = 1; k <= N; k++)
    {
      *cheb_n++ = 0.5 * (*f_hat_p++ + *f_hat_n--);
    }
  }  
}
