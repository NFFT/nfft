/*! \file fastsum_matlab.c
 *  \brief Simple test program for the fast NFFT-based summation algorithm, called by fastsum.m.
 *
 *  \author Markus Fenn
 *  \date 2006
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "fastsum.h"
#include "kernels.h"

int main(int argc, char **argv)
{
  int j,k,t;                                         /**< indices                 */
  int d;                                             /**< number of dimensions    */
  int N;                                             /**< number of source nodes  */
  int M;                                             /**< number of target nodes  */
  int n;                                             /**< expansion degree        */
  int m;                                             /**< cut-off parameter       */
  int p;                                             /**< degree of smoothness    */
  char *s;                                           /**< name of kernel          */
  complex (*kernel)(double , int , const double *);  /**< kernel function         */
  double c;                                          /**< parameter for kernel    */
  fastsum_plan my_fastsum_plan;                      /**< plan for fast summation */
  complex *direct;                                   /**< array for direct computation */
  double time;                                       /**< for time measurement    */
  double error=0.0;                                  /**< for error computation   */
  double eps_I;                                      /**< inner boundary          */
  double eps_B;                                      /**< outer boundary          */
  FILE *fid1, *fid2;
  double temp;

  if (argc!=11)
  {
    printf("\nfastsum_test d N M n m p kernel c\n\n");
    printf("  d       dimension                 \n");
    printf("  N       number of source nodes    \n");
    printf("  M       number of target nodes    \n");
    printf("  n       expansion degree          \n");
    printf("  m       cut-off parameter         \n");
    printf("  p       degree of smoothness      \n");
    printf("  kernel  kernel function  (e.g., gaussian)\n");
    printf("  c       kernel parameter          \n");
    printf("  eps_I   inner boundary            \n");
    printf("  eps_B   outer boundary            \n\n");
    exit(-1);
  }
  else
  {
    d=atoi(argv[1]);
    N=atoi(argv[2]); c=1.0/pow((double)N,1.0/(double)d);
    M=atoi(argv[3]);
    n=atoi(argv[4]);
    m=atoi(argv[5]);
    p=atoi(argv[6]);
    s=argv[7];
    c=atof(argv[8]);
    eps_I=atof(argv[9]);
    eps_B=atof(argv[10]);
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
  printf("d=%d, N=%d, M=%d, n=%d, m=%d, p=%d, kernel=%s, c=%g, eps_I=%g, eps_B=%g \n",d,N,M,n,m,p,s,c,eps_I,eps_B);

  /** init two dimensional fastsum plan */
  fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, 0, n, m, p, eps_I, eps_B);
  /*fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, EXACT_NEARFIELD, n, m, p);*/

  /** load source knots and coefficients */
  fid1=fopen("x.dat","r");
  fid2=fopen("alpha.dat","r");
  for (k=0; k<N; k++)
  {
    for (t=1; t<d; t++)
    {
      fscanf(fid1,"%le",&my_fastsum_plan.x[k*d+t]);
    }
    fscanf(fid2,"%le",&temp); my_fastsum_plan.alpha[k] = temp;
    fscanf(fid2,"%le",&temp); my_fastsum_plan.alpha[k] += temp*I;
  }
  fclose(fid1);
  fclose(fid2);

  /** load target knots */
  fid1=fopen("y.dat","r");
  for (j=0; j<M; j++)
  {
    for (t=1; t<d; t++)
    {
      fscanf(fid1,"%le",&my_fastsum_plan.y[j*d+t]);
    }
  }
  fclose(fid1);

  /** direct computation */
  printf("direct computation: "); fflush(NULL);
  time=nfft_second();
  fastsum_exact(&my_fastsum_plan);
  time=nfft_second()-time;
  printf("%fsec\n",time);

  /** copy result */
  direct = (complex *)malloc(my_fastsum_plan.M_total*(sizeof(complex)));
  for (j=0; j<my_fastsum_plan.M_total; j++)
    direct[j]=my_fastsum_plan.f[j];

  /** precomputation */
  printf("pre-computation:    "); fflush(NULL);
  time=nfft_second();
  fastsum_precompute(&my_fastsum_plan);
  time=nfft_second()-time;
  printf("%fsec\n",time);

  /** fast computation */
  printf("fast computation:   "); fflush(NULL);
  time=nfft_second();
  fastsum_trafo(&my_fastsum_plan);
  time=nfft_second()-time;
  printf("%fsec\n",time);

  /** compute max error */
  error=0.0;
  for (j=0; j<my_fastsum_plan.M_total; j++)
  {
   if (cabs(direct[j]-my_fastsum_plan.f[j])/cabs(direct[j])>error)
      error=cabs(direct[j]-my_fastsum_plan.f[j])/cabs(direct[j]);
  }
  printf("max relative error: %e\n",error);

  /** write result to file */
  fid1=fopen("f.dat","w+");
  fid2=fopen("f_direct.dat","w+");
  if (fid1==NULL)
  {
    printf("Fehler!\n");
    exit(-1);
  }
  for (j=0; j<M; j++)
  {
    temp=creal(my_fastsum_plan.f[j]);
    fprintf(fid1,"  % .16e",temp);
    temp=cimag(my_fastsum_plan.f[j]);
    fprintf(fid1,"  % .16e\n",temp);

    temp=creal(direct[j]);
    fprintf(fid2,"  % .16e",temp);
    temp=cimag(direct[j]);
    fprintf(fid2,"  % .16e\n",temp);
  }
  fclose(fid1);
  fclose(fid2);

  /** finalise the plan */
  fastsum_finalize(&my_fastsum_plan);

  return 0;
}
