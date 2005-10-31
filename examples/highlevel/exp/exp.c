#include "util.h"
#include "options.h"
#include "window_defines.h"
#include "stdlib.h"
#include "math.h"

#include "nfft3.h"

/*
 * window_funct_plan is a plan to use the window functions
 * independent of the nfft
 */

typedef struct window_funct_plan_ {
  int d;
	int m;
	int n[1];
	double sigma[1];
	double *b;
  double *spline_coeffs;                /**< input for de Boor algorithm, if   
                                             B_SPLINE or SINC_2m is defined   */ 
} window_funct_plan;

void window_funct_init(window_funct_plan* ths, int m, int n, double sigma) {
	ths->d=1;
	ths->m=m;
	ths->n[0]=n;
	ths->sigma[0]=sigma;
  WINDOW_HELP_INIT
}

void exponential (double k)
{
  int l,n,N;
  double xx;
  window_funct_plan  *ths;
  complex result_exact,result_approx;

  ths = (window_funct_plan*) malloc(sizeof(window_funct_plan));

	N=ceil(k);
  n=2*N;
  
  window_funct_init(ths,2, n,n/((double)N));

  result_approx= 0.0;
  for (xx=-0.5;xx<0.5;xx=xx+0.001)
  {
    for (l = -n/2;l<n/2;l++)
    {
      result_approx += PHI(xx -(((double)l)/((double)n)),0)  * cexp( -2.0*PI*I*k*l/((double)n));
    }

    result_approx *= 1.0/(PHI_HUT(k,0));

    result_exact = cexp(-2.0*PI*I*k*xx );

    printf("%e %e %e %e %e %e\n",xx,creal(result_exact),creal(result_approx),
        cimag(result_exact),cimag(result_approx) ,cabs(result_exact-result_approx));
  }
  WINDOW_HELP_FINALIZE
	free(ths);
}


int main(int argc, char **argv)
{
  // usage: k 
  exponential(atof(argv[1]));
  return 1;
}

