#include "util.h"
#include "options.h"
#include "window_defines.h"
#include "stdlib.h"
#include "math.h"

#include "nfft3.h"

#define PHI_periodic(x) ((x>0.5)?(PHI(x-1.0,0)):((x<-0.5)?PHI(x+1.0,0):PHI(x,0)))

void exponential (double k, double x)
{
  int l,n,N;
  double xx;
  nfft_plan  *ths;
  complex result_exact,result_approx;

  ths = (nfft_plan*) malloc(sizeof(nfft_plan));

  N=2*ceil(k);
  
  nfft_init_1d(ths,N,1);

  n=ths->n[0];

  result_approx= 0.0;
  for (xx=-x;xx<x;xx=xx+0.001)
  {
    for (l = -n/2;l<n/2;l++)
    {
      result_approx += PHI_periodic(xx -(((double)l)/((double)n)))  * cexp( -2.0*PI*I*k*l/((double)n));
    }

    result_approx *= 1.0/(PHI_HUT(k,0));

    result_exact = cexp(-2.0*PI*I*k*xx );

    printf("%e %e %e %e\n",xx,creal(result_exact),creal(result_approx) ,cabs(result_exact-result_approx));
  }
  nfft_finalize(ths);
  free(ths);
}


int main(int argc, char **argv)
{
  // usage: k x 
  exponential(atof(argv[1]),atof(argv[2]));
  return 1;
}

