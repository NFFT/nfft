/* 
   nfsft - Spherical Fourier transform of regularly and irregularly sampled data

   Copyright (C) 2004 Jens Keiner

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

#include <termios.h>
#include <grp.h>
#include <pwd.h>
*/
#ifdef HAVE_CONFIG_H
#  include "config.h"
#else
#  error Need config.h
#endif

#ifdef STDC_HEADERS
#  include <stdlib.h>
#  include <stdio.h>
#  include <float.h>
#  include <math.h>
#else
#  error Need ANSI-C headers
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#else
#  error Need sys/types.h
#endif

/* Auxilliary headers */
#include <complex.h>
#include <fftw3.h>
#include "nfsft.h"
#include "util.h"

#define max(a,b) a<b?b:a

/** Arrays for complex Fourier-coefficients. */
complex *f_hat;
complex *f_hat2;
complex **F_HAT;
complex **F_HAT2;


/**
*
 */
void get_input(int* D, int* M)
{ 
  //printf("\nComputation of f(theta_d,phi_d) for d=0,...,D-1\n");
  
  //printf("Number of nodes D = ");
  //fflush(stdout); 
  //scanf("%d", D);
  *D = 4096;//2048;
  
  //printf("Bandwidth M = "); 
  //fflush(stdout); 
  //scanf("%d", M);
  *M = 256;//256;
           //*D = *M + 1;//2048;
    
    //printf("Threshold for basis-transformation = "); 
    //fflush(stdout); 
    //scanf("%ld", threshold);
    //*threshold = 100000.0;
    
    //printf("Cut off gaussian bell for nfft = "); 
    //fflush(stdout); 
    //scanf("%d", kernelWidth);
    //*kernelWidth = 4;
    
    //printf("\n");
}

void free_data(int M)
{
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


/**
 * The main program.
 *
 * \arg argc The number of arguments
 * \arg argv An array containing the arguments as C-strings
 *
 * \return 
 */
int main (int argc, char **argv)
{  
  /** The bandwidth */
  int M;
  /** The number of nodes */
  int D; 
  /** Next greater power of 2 relative to M */
  int N;
  /** Array of angles phi of data points. */
  //double *nodes_phi;
  /** Array of angles theta of data points. */
  //double *nodes_theta;
  double *angles;
  /** Plan for fast spherical fourier transform. */
  nfsft_plan plan;
  /** Plan for fast spherical fourier transform. */
  nfsft_plan plan2;
  /** Fix random seed for testing purposes. */
  unsigned short seed[3] = {1,2,3};
  complex *result;
  complex *result2;
  double ctime;
  int i;  
  int n,k,nleg;
  int mode;
  
  //nfsft_test();
  //return EXIT_SUCCESS;
  
//#define DEBUG
  fscanf(stdin,"%d\n",&mode);

  /* Initialize nodes. */
  fscanf(stdin,"%d\n",&D);
  angles = (double*) calloc(2*D,sizeof(double));
  for (i = 0; i < D; i++)
  {
    fscanf(stdin,"%lf\n%lf\n",&angles[2*i],&angles[2*i+1]);    
  } 
#ifdef DEBUG
  myvpr(angles,2*D,"angles");
#endif
  
  /* Initialize Fourier coefficients. */
  fscanf(stdin,"%d\n",&M);
#ifdef DEBUG
  printf("M = %d\n",M);
#endif
  /* Calculate N as next greater power of 2 of the bandwidth M. */
  N = M==0?0:1<<(int)ceil(log2((double)M));
#ifdef DEBUG
  printf("N = %d\n",N);
#endif
  
  /* Allocate memory for data. */
  F_HAT = (complex**) calloc(2*M+1,sizeof(complex*));
  F_HAT2 = (complex**) calloc(2*M+1,sizeof(complex*));
  for (n = -M; n <= M; n++)
  {
    f_hat = (complex*) calloc(N+1,sizeof(complex));
    f_hat2 = (complex*) calloc(N+1,sizeof(complex));    
    F_HAT[n+M] = f_hat;
    F_HAT2[n+M] = f_hat2;    
  }

  result = (complex*) calloc(D,sizeof(complex));
  result2 = (complex*) calloc(D,sizeof(complex));   
  
  if (mode == 0)
  {      
    /* Read Fourier coefficients from standard input. */    
    double c;
    for (k = 0; k <= M; k++)
    {
      for (n = -k; n <= k; n++)
      {
        fscanf(stdin,"%lf\n",&c);
        F_HAT[n+M][k] = c;
      }  
    }  
    for (n = -M; n <= M; n++)
    {
      memcpy(F_HAT2[n+M],F_HAT[n+M],(N+1)*sizeof(complex));
    }  

    /* Create plan for fast spherical Fourier transform.*/
    plan = nfsft_create_plan(NFSFT_BACKWARD,D,M,angles,F_HAT,
      result,0U);
    /* Initialize */
    ctime = mysecond();
    init_wisdom(M);
    ctime = mysecond() - ctime;
#ifdef DEBUG
    printf("Time for initialization: %f sec.\n",ctime);
#endif
    /* Precompute. */
    ctime = mysecond();
    nfsft_compute_wisdom(M);
    ctime = mysecond() - ctime;
#ifdef DEBUG
    printf("Time for precomputation: %f sec.\n",ctime);
#endif
    /* Execute the plan. */
    ctime = mysecond();
    nfsft_execute(plan);
    ctime = mysecond() - ctime;
#ifdef DEBUG
    printf("Time for fast algorithm: %f sec.\n",ctime);
#endif
    /* Create plan for fast spherical Fourier transform.*/
    plan2 = nfsft_create_plan(NFSFT_BACKWARD,D,M,angles,
      F_HAT2,result2,0);
    /* Execute the plan. */
    ctime = mysecond();
    ndsft_execute(plan2);
    ctime = mysecond() - ctime;
//#define DEBUG
#ifdef DEBUG
    printf("Time for slow algorithm: %f sec.\n",ctime);
    printf("relative error:%e\n",error_complex_inf(result2,result,D));  
    printf("relative error:%e\n",error_complex_1(result2,result,D));  
    //myvprc(result,D,"result");
#endif
    for (k = 0; k < D; k++)
    {
      printf("%17.16f\n%17.16f\n",creal(result[k]),cimag(result[k]));
    }
  }
  else
  {
    /* Read function values from standard input. */    
    double c;
    for (k = 0; k < D; k++)
    {
      fscanf(stdin,"%lf\n",&c);
      result[k] = c;
      result2[k] = c;
    }  
    
    /* Create plan for fast spherical Fourier transform.*/
    plan = nfsft_create_plan(NFSFT_ADJOINT,D,M,angles,F_HAT,
                             result,0U);
    /* Initialize */
    ctime = mysecond();
    init_wisdom(M);
    ctime = mysecond() - ctime;
#ifdef DEBUG
    printf("Time for initialization: %f sec.\n",ctime);
#endif
    /* Precompute. */
    ctime = mysecond();
    nfsft_compute_wisdom(M);
    ctime = mysecond() - ctime;
#ifdef DEBUG
    printf("Time for precomputation: %f sec.\n",ctime);
#endif
    /* Execute the plan. */
    ctime = mysecond();
    nfsft_execute(plan);
    ctime = mysecond() - ctime;
#ifdef DEBUG
    printf("Time for fast algorithm: %f sec.\n",ctime);
#endif
    /* Create plan for fast spherical Fourier transform.*/
    /*plan2 = nfsft_create_plan(NFSFT_ADJOINT,D,M,angles,
                              F_HAT2,result2,0);*/
    /* Execute the plan. */
    /*ctime = mysecond();
    ndsft_execute(plan2);
    ctime = mysecond() - ctime;
    printf("Time for slow algorithm: %f sec.\n",ctime);
    printf("relative error:%e\n",error_complex_inf(result2,result,D));  
    printf("relative error:%e\n",error_complex_1(result2,result,D));  */
#define DEBUG
#ifdef DEBUG
    for (k = 0; k <= M; k++)
    {  
      for (n = -k; n <= k; n++)
      {
        printf("%17.16f\n%17.16f\n",creal(F_HAT[n+M][k]),cimag(F_HAT[n+M][k]));
      }
      //myvprc(F_HAT[k+M],N+1,"f_hat");
    }
#endif
#undef DEBUG
  }
  
  /* Destroy the plan. */
  nfsft_destroy_plan(plan);
  /* Destroy the plan. */
  nfsft_destroy_plan(plan2);
  
  /* Forget wisdom. */
  nfsft_forget_wisdom();
  
  fftw_free(result);
  fftw_free(result2);
    
  /* Free data vectors. */
  for (n=-M;n<=M;n++)
  {
    fftw_free(F_HAT[n+M]);
    fftw_free(F_HAT2[n+M]);
  }
  
  fftw_free(F_HAT);
  fftw_free(F_HAT2);
  
  /* Free nodes. */
  free(angles);
  
  return EXIT_SUCCESS;
}
