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
  double kk;
  nfft_plan  *ths;
  /* ergeb ist die exakte Berechnung
     ergeb2 die Approximation 
   */
  complex result_exact,result_approx;

  ths = (nfft_plan*) malloc(sizeof(nfft_plan));

  N=2*ceil(k);
  
  /* Nur fuer PHI und PHI_HUT es ist n=2* N/2 = N */
  nfft_init_1d(ths,N,1);

  n=ths->n[0];

  result_approx= 0.0;
  for (kk=-k;kk<k;kk=kk+0.1)
  {
  for (l = -n/2;l<n/2;l++)
  {
    /* Periodisierung beachten */
     result_approx += PHI_periodic(x -(((double)l)/((double)n)))  * cexp( -2.0*PI*I*kk*l/((double)n));
  }

  result_approx *= 1.0/(PHI_HUT(kk,0));

  result_exact = cexp(-2.0*PI*I*kk*x );

  /*  printf("%e %e \n",creal(ergeb),cimag(ergeb));
      printf("%e %e \n",creal(ergeb2),cimag(ergeb2)); */
  /*printf("%e %e %e \n",kk, creal(ergeb)-creal(ergeb2),cimag(ergeb)-cimag(ergeb2));*/
  printf("%e %e  \n",kk, cabs(result_exact-result_approx));
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

