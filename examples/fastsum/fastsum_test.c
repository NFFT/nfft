#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fastsum.h"
#include "kernels.h"

/**********************************************************************
 **********************************************************************/
int main(int argc, char **argv)
{
  int j,k,t;                                         /**< indices                 */
  int d;                                             /**< number of dimensions    */
  int N;                                             /**< number of source nodes  */
  int M;                                             /**< number of target nodes  */
  int n;                                             /**< expansion degree        */
  int p;                                             /**< degree of smoothness    */
  char *s;                                           /**< name of kernel          */
  complex (*kernel)(double , int , const double *);  /**< kernel function         */
  double c;                                          /**< parameter for kernel    */
  fastsum_plan my_fastsum_plan;                      /**< plan for fast summation */
  complex *direct;                                   /**< array for direct computation */
  double time;                                       /**< for time measurement    */
  double error=0.0;                                  /**< for error computation   */
  double *temp;                                      /**< for temporary things    */

  if (argc!=8)
  {
    printf("\nfastsum_test d N M n p kernel c\n\n");
    printf("  d       dimension                 \n");
    printf("  N       number of source nodes    \n");
    printf("  M       number of target nodes    \n");
    printf("  n       expansion degree          \n");
    printf("  p       degree of smoothness      \n");
    printf("  kernel  kernel function  (e.g., gaussian)\n");
    printf("  c       kernel parameter          \n\n");
    exit(-1);
  }
  else
  {
    d=atoi(argv[1]);
    N=atoi(argv[2]); c=1.0/pow((double)N,1.0/(double)d);
    M=atoi(argv[3]);
    n=atoi(argv[4]);
    p=atoi(argv[5]);
    s=argv[6];
    c=atof(argv[7]);
    if (strcmp(s,"gaussian")==0)
      kernel = gaussian;
    else if (strcmp(s,"multiquadric")==0)
      kernel = multiquadric;
    else if (strcmp(s,"inverse_multiquadric")==0)
      kernel = inverse_multiquadric;
    else if (strcmp(s,"logarithm")==0)
      kernel = logarithm;
    else if (strcmp(s,"thinplate_spline")==0)
      kernel = thinplate_spline;
    else if (strcmp(s,"one_over_square")==0)
      kernel = one_over_square;
    else if (strcmp(s,"one_over_modulus")==0)
      kernel = one_over_modulus;
    else if (strcmp(s,"one_over_x")==0)
      kernel = one_over_x;
    else if (strcmp(s,"inverse_multiquadric3")==0)
      kernel = inverse_multiquadric3;
    else if (strcmp(s,"sinc_kernel")==0)
      kernel = sinc_kernel;
    else if (strcmp(s,"cosc")==0)
      kernel = cosc;
    else if (strcmp(s,"cot")==0)
      kernel = cot;
    else
    {
      s="multiquadric";
      kernel = multiquadric;
    }
  }
  printf("d=%d, N=%d, M=%d, n=%d, p=%d, kernel=%s, c=%g \n",d,N,M,n,p,s,c);

  temp = (double *)malloc(N*(sizeof(double)));

  /** init two dimensional fastsum plan */
  fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, 0, n, p);

  /** init source knots in a d-ball with radius 0.25-eps_b/2 */
  for (k=0; k<N; k++)
  {
    double r=(0.25-my_fastsum_plan.eps_B/2.0)*pow((double)rand()/(double)RAND_MAX,1.0/d);
    my_fastsum_plan.x[k*d+0] = r;
    for (j=1; j<d; j++)
    {
      double phi=2.0*PI*(double)rand()/(double)RAND_MAX;
      my_fastsum_plan.x[k*d+j] = r;
      for (t=0; t<j; t++)
      {
        my_fastsum_plan.x[k*d+t] *= cos(phi);
      }
      my_fastsum_plan.x[k*d+j] *= sin(phi);
    }

    my_fastsum_plan.alpha[k] = (double)rand()/(double)RAND_MAX + I*(double)rand()/(double)RAND_MAX;
  }

  /** init target knots in a d-ball with radius 0.25-eps_b/2 */
  for (k=0; k<M; k++)
  {
    double r=(0.25-my_fastsum_plan.eps_B/2.0)*pow((double)rand()/(double)RAND_MAX,1.0/d);
    my_fastsum_plan.y[k*d+0] = r;
    for (j=1; j<d; j++)
    {
      double phi=2.0*PI*(double)rand()/(double)RAND_MAX;
      my_fastsum_plan.y[k*d+j] = r;
      for (t=0; t<j; t++)
      {
        my_fastsum_plan.y[k*d+t] *= cos(phi);
      }
      my_fastsum_plan.y[k*d+j] *= sin(phi);
    }
  }

  /** direct computation */
  printf("direct computation: "); fflush(NULL);
  time=second();
  fastsum_exact(&my_fastsum_plan);
  time=second()-time;
  printf("%fsec\n",time);

  /** copy result */
  direct = (complex *)malloc(my_fastsum_plan.M_total*(sizeof(complex)));
  for (j=0; j<my_fastsum_plan.M_total; j++)
    direct[j]=my_fastsum_plan.f[j];

  /** precomputation */
  printf("pre-computation:    "); fflush(NULL);
  time=second();
  fastsum_precompute(&my_fastsum_plan);
  time=second()-time;
  printf("%fsec\n",time);

  /** fast computation */
  printf("fast computation:   "); fflush(NULL);
  time=second();
  fastsum_trafo(&my_fastsum_plan);
  time=second()-time;
  printf("%fsec\n",time);

  /** compute max error */
  error=0.0;
  for (j=0; j<my_fastsum_plan.M_total; j++)
  {
    if (cabs(direct[j]-my_fastsum_plan.f[j])/cabs(direct[j])>error)
      error=cabs(direct[j]-my_fastsum_plan.f[j])/cabs(direct[j]);
  }
  printf("max relative error: %e\n",error);

  /** finalise the plan */
  fastsum_finalize(&my_fastsum_plan);

  return 0;
}
