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

