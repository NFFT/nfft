#include "config.h"

#include <fftw3.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

#include "util.h"

inline int ngpt(int n)
{
  int e = 0;
  int p = 1;
  
  while (n > p)
  {
    e++;
    p = p<<1;
  }
  return e;  
}

inline void vpr_c_hat (complex **x, int n, char *text)
{
  int k,m,d;
  
  fprintf(stderr,"%s, adr=%p\n", text, x);
  
  d = 0;
  for (k=0; k<=n; k++)
  {
    for (m = -k; m <= k; m++)
    {  
      if (d%2==0) fprintf(stderr,"%d,%d.", k,m);
      fprintf(stderr,"(%20.16e,", creal(x[k][m+n]));
      fprintf(stderr,"%20.16e),", cimag(x[k][m+n]));
      if (d%2==1) fprintf(stderr,"\n");
      fflush(stdout);
      d++;
    }
  }
  
  if (n%5!=0)
    fprintf(stderr,"\n");
}

double error_complex_inf_r_old(complex *x0, complex *x, int n)
{
  double max_oldimum=0.0, max_olddiff=0.0, xda;
  int k;
  double res;
  double x0a,xa;
  
  for (k = 0; k < n; k++) 
  {
    x0a = cabs(x0[k]);
    xa = cabs(x[k]);
    max_oldimum = MAX(MAX(x0a,xa),max_oldimum);
    xda = cabs(x0[k]-x[k]);
    max_olddiff = MAX(max_olddiff,xda);
  }
  
  res = (max_oldimum<2*DBL_EPSILON ? max_olddiff : max_olddiff / max_oldimum);
  
  return res;
}

double error_complex_inf_old(complex *x0, complex *x, int n)
{
  double maximum1=0.0;
  double maximum2=0.0;
  int k;
  
  for (k = 0; k < n; k++) 
  {
    maximum1 = MAX(maximum1,cabs(x0[k]-x[k]));
    maximum2 = MAX(maximum2,cabs(x0[k]));
  }
    
  return maximum1/maximum2;
}

double norm_complex_inf_old(complex *x, int n)
{
  double max_oldimum=0.0;
  int k;
  
  for (k = 0; k < n; k++) 
  {
    max_oldimum = MAX(max_oldimum,cabs(x[k]));
  }
  
  return max_oldimum;
}

double norm_complex_2_old(complex *x, int n)
{
  double r = 0.0;
  int k;
  
  for (k = 0; k < n; k++) 
  {
    r += cabs(x[k])*cabs(x[k]);
  }
  
  return sqrt(r);
}

double norm_complex_1_old(complex *x, int n)
{
  double r = 0.0;
  int k;
  
  for (k = 0; k < n; k++) 
  {
    r += cabs(x[k]);
  }
  
  return r;
}

double error_complex_1_old(complex *a, complex *b, int n)
{
  double l1a,l1d;
  int k;
  
  l1a = 0.0;
  l1d = 0.0;
  
  for (k = 0; k < n; k++) 
  {
    l1d += cabs(a[k]-b[k]);
    l1a += cabs(a[k]);
  }
    
  return l1d/l1a;
}

double error_complex_2_old(complex *a, complex *b, int n)
{
  double l1a,l1d;
  int k;
  
  l1a = 0.0;
  l1d = 0.0;
  
  for (k = 0; k < n; k++) 
  {
    l1d += cabs(a[k]-b[k])*cabs(a[k]-b[k]);
    l1a += cabs(a[k])*cabs(a[k]);
  }
    
  return sqrt(l1d/l1a);
}

double error_complex3_old(complex *a, complex *b, int n)
{
  double m;
  int k;
  
  m = 0.0;
  
  for (k = 0; k < n; k++) 
  {
    m = MAX(m,cabs(a[k]-b[k])/cabs(a[k]));
  }
  
  return m;
}

double err_f_hat_infty_old(complex **f_hat, complex **f_hat2, int M)
{
  double err;
  double absmax_old;
  int k,n;
  
  err = 0.0;
  absmax_old = 0.0;
  
  for (k = 0; k <= M; k++) 
  {
    for (n = -k; n <= k; n++) 
    {
      absmax_old = MAX(absmax_old,cabs(f_hat[n+M][k]));
      err = MAX(err,cabs(f_hat[n+M][k]-f_hat2[n+M][k]));
    }
  }  
  return err/absmax_old;
}

double err_f_hat_1_old(complex **f_hat, complex **f_hat2, int M)
{
  double err = 0;
  double temp = 0;
  int k,n;
  
  for (k = 0; k <= M; k++) 
  {
    for (n = -k; n <= k; n++) 
    {
      err += cabs(f_hat[n+M][k]-f_hat2[n+M][k]);
      temp += cabs(f_hat[n+M][k]);
    }
  }  
  return err/temp;
}

double norm_nfft_1_old(complex *f_hat, int M)
{
  double r = 0.0;
  int k,n;
	complex* t = f_hat;
  
  for (k = 0; k < M; k++) 
  {
    for (n = 0; n < M; n++) 
    {
		  //fprintf(stderr,"r = %.4E\n",r);
      r += cabs(*t++);
    }
  }  
  return r;
}

double norm_f_hat_1_old(complex **f_hat, int M)
{
  double r = 0.0;
  int k,n;
  
  for (k = 0; k <= M; k++) 
  {
    for (n = -k; n <= k; n++) 
    {
      r += cabs(f_hat[n+M][k]);
    }
  }  
  return r;
}

double norm_f_hat_1s_old(complex **f_hat, int M)
{
  double r = 0.0;
  int k,n;
  
  for (k = 1; k <= M; k++) 
  {
    for (n = -k; n <= k; n++) 
    {
      r += cabs(f_hat[n+M][k])/sqrt((double)k);
    }
  }  
  r += cabs(f_hat[0][0]);
  return r;
}

double err_f_hat_2_old(complex **f_hat, complex **f_hat2, int M)
{
  double err = 0;
  double temp = 0;
  int k,n;
  
  for (k = 0; k <= M; k++) 
  {
    for (n = -k; n <= k; n++) 
    {
      err += cabs(f_hat[n+M][k]-f_hat2[n+M][k])*cabs(f_hat[n+M][k]-f_hat2[n+M][k]);
      temp += cabs(f_hat[n+M][k])*cabs(f_hat[n+M][k]);
    }
  }  
  return sqrt(err/temp);
}

/** Copies \f$x \leftarrow y\f$.
*/
void copyc_hat_old(complex **x,complex **y, int M)
{
  int n,k;
  
  for (k = 0; k <= M; k++)
  {
    for (n = -k; n <= k; n++)
    {
      x[n+M][k] = y[n+M][k];
    }
  }
}


/** Copies \f$x \leftarrow w\odot y\f$.
*/
void copyc_w_hat_old(complex **x, double **w, complex **y, int M)
{
  int n,k;
  
  for (k = 0; k <= M; k++)
  {
    for (n = -k; n <= k; n++)
    {
      x[n+M][k] = w[n+M][k] * y[n+M][k];
    }
  }
}

/** Computes the inner/dot product \f$x^H x\f$.
*/
double dotproductc_hat_old(complex** x, int M)
{
  int n,k;
  double result=0.0;
  
  for(k = 0; k <= M; k++)
  {
    for(n = -k; n <= k; n++)
    {
      result += (x[n+M][k]*x[n+M][k]);
    }
  }
  
  return result;
}

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double dotproductc_w_hat_old(complex **x, double **w, int M)
{
  int n,k;
  double result=0.0;
  
  for(k = 0; k <= M; k++)
  {
    for(n = -k; n <= k; n++)
    {
      result += w[n+M][k] * (x[n+M][k]*x[n+M][k]);
    }
  }
  
  return result;
}

/** Updates \f$x \leftarrow x + a y\f$.
*/
inline void updatec_xpay_hat_old(complex **x,double a, complex **y, int M)
{
  int n,k;
  
  for(k = 0; k <= M; k++)
  {
    for(n = -k; n <= k; n++)
    {
      x[n+M][k] += a*y[n+M][k];
    }
  }
}

/** Updates \f$x \leftarrow x + a w\odot y\f$.
*/
inline void updatec_xpawy_hat_old(complex **x, double a, double **w, complex **y,
                   int M)
{
  int n,k;
  
  for(k = 0; k <= M; k++)
  {
    for(n = -k; n <= k; n++)
    {
      x[n+M][k] += a*w[n+M][k]*y[n+M][k];
    }
  }
}

/** Updates \f$x \leftarrow a x + y\f$.
*/
inline void updatec_axpy_hat_old(complex **x,double a, complex **y, int M)
{
  int n,k;
  
  for(k = 0; k <= M; k++)
  {
    for(n = -k; n <= k; n++)
    {
      x[n+M][k]=a*x[n+M][k]+y[n+M][k];
    }
  }
}

/** Updates \f$x \leftarrow a x + b y\f$.
*/
inline void updatec_axpby_hat_old(complex **x, double a, complex **y, double b, int M)
{
  int n,k;
  
  for(k = 0; k <= M; k++)
  {
    for(n = -k; n <= k; n++)
    {
      x[n+M][k] = a*x[n+M][k] + b*y[n+M][k];
    }
  }
}

/** Updates \f$x \leftarrow a x + y\f$.
*/
inline void updatec_axpy_2_old(complex* x,double a, complex* y, int n)
{
  int l;
  
  for(l=0;l<n;l++)
  {
    x[l] = a*x[l] + y[l];
  }
}

/** Updates \f$x \leftarrow x + a y\f$.
*/
inline void updatec_xpay_2_old(complex* x,double a, complex* y, int n)
{
  int l;
  
  for(l=0;l<n;l++)
  {
    x[l] += a*y[l];
  }
}

/** Updates \f$x \leftarrow a x + b y\f$.
*/
inline void updatec_axpby_2(complex* x, double a, complex* y, double b,
                   int n)
{
  int l;
  
  for(l=0;l<n;l++)
  {
    x[l]=a*x[l]+b*y[l];
  }
}

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
inline void updatec_xpawy_2_old(complex* x,double a, double* w, complex* y,
                   int n)
{
  int l;
  
  for(l=0;l<n;l++)
  {
    x[l] += a*w[l]*y[l];
  }
}

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
inline void uvxpwy_old(complex* u, complex* x, double* v, complex* y, double* w, int n)
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
    *u++ = (*v++) * (*x++) + (*w++) * (*y++);
  }
}

/** Updates \f$x \leftarrow x + a w\odot y\f$.
*/
inline void auvxpwy_old(double a, complex* u, complex* x, double* v, complex* y, double* w, int n)
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

inline void abuvxpwy_old(double a, double b, complex* u, complex* x, double* v, complex* y, double* w, int n)
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
    *u++ = a * (b * (*v++) * (*x++) + (*w++) * (*y++));
  }
}

inline void normalize_f_hat_old(complex **f_hat, int M)
{
  int k,n;
  for (k = 0; k <= M; k++)
  {
    for (n = -k; n <= k; n++)
    {
      f_hat[n+M][k] *= sqrt((2*k+1)/2.0);
    }
  }
}

inline int min_old(const int a, const int b)
{
  return a<b?a:b;
}

inline int max_old(const int a, const int b)
{
  return a>b?a:b;
}

inline long int max_oldl(const long int a, const long int b)
{
  return a>b?a:b;
}

inline int abs(const int a)
{
  return (a>0)?(a):(-a);
}

