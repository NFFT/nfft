#include "config.h"

#ifdef HAVE_FFTW3_H
# include <fftw3.h>
#else
# error Need fftw3.h
#endif

#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

#include "util.h"

#define max(a,b) a<b?b:a

inline int ngpt(int n)
{
  return (int) ceil(log((double)n)/log(2.0));
}

/** For timing */ 
inline double mysecond()
{
  static struct rusage temp;
  double foo, foo1;
  
  getrusage(RUSAGE_SELF,&temp)    ;
  foo     = temp.ru_utime.tv_sec  ;       /* seconds */
  foo1    = temp.ru_utime.tv_usec ;       /* uSecs */
  return  foo  + (foo1/1000000.0)         ;       /* milliseconds */
}

inline void myvpr (double* x, int n, char * text)
{
  int k;
  
  printf ("%s, adr=%p\n", text, x);
  
  for (k=0; k<n; k++)
  {
    if (k%5==0) printf("%d.", k);
    printf("%14.7e,", x[k]);
    if (k%5==4) printf("\n");
    fflush(stdout);
  }
  
  if (n%5!=0)
  {
    printf("\n");
  }
}


inline void myvprs (double* x, int n, char * text, int stride)
{
  int k;
  
  printf ("%s, adr=%p\n", text, x);
  
  for (k=0; k<n; k++)
  {
    if (k%5==0) printf("%d.", k);
    printf("%14.7e,", x[stride*k]);
    if (k%5==4) printf("\n");
    fflush(stdout);
  }
  
  if (n%5!=0)
  {
    printf("\n");
  }
}

inline void myvprc (complex *x, int n, char *text)
{
  int k;
  
  printf ("%s, adr=%p\n", text, x);
  
  for (k=0; k<n; k++)
  {
    if (k%5==0) printf("%d.", k);
    printf("(%14.7e,", creal(x[k]));
    printf("%14.7e),", cimag(x[k]));
    if (k%5==4) printf("\n");
    fflush(stdout);
  }
  
  if (n%5!=0)
    printf("\n");
}

inline void myvprci (complex *x, int n, char *text)
{
  int k;
  double *myx;
  myx = (double*) x;
  
  printf ("%s, adr=%p\n", text, x);
  
  for (k=0; k<2*n; k++)
  {
    if (k%5==0) printf("%d.", k);
    printf("%14.7e,", myx[k]);
    if (k%5==4) printf("\n");
    fflush(stdout);
  }
  
  if (n%5!=0)
    printf("\n");
}

inline int pow2(const int t)
{  
  return (1L << t); 
}


double error_complex_inf(complex *x0, complex *x, int n)
{
  double maximum=0.0, maxdiff=0.0, xda;
  int k;
  double res;
  double x0a,xa;
  
  for (k = 0; k < n; k++) 
  {
    x0a = cabs(x0[k]);
    xa = cabs(x[k]);
    maximum = max(max(x0a,xa),maximum);
    xda = cabs(x0[k]-x[k]);
    maxdiff = max(maxdiff,xda);
  }
  
  res = (maximum<2*DBL_EPSILON ? maxdiff : maxdiff / maximum);
  
  return res;
}

double error_complex_1(complex *a, complex *b, int n)
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

double error_complex3(complex *a, complex *b, int n)
{
  double m;
  int k;
  
  m = 0.0;
  
  for (k = 0; k < n; k++) 
  {
    m = max(m,cabs(a[k]-b[k])/cabs(a[k]));
  }
  
  return m;
}

double error_f_hat(complex **f_hat, complex **f_hat2, int M)
{
  double err;
  double absmax;
  int k,n;
  
  err = 0.0;
  absmax = 0.0;
  
  for (k = 0; k <= M; k++) 
  {
    for (n = -k; n <= k; n++) 
    {
      absmax = max(absmax,cabs(f_hat[n+M][k]));
      err = max(err,cabs(f_hat[n+M][k]-f_hat2[n+M][k]));
    }
  }
  
  return err/absmax;
}


