/* 
   accuracy - Fast spherical convolution example for the NFSFT

   Copyright (C) 2005 Jens Keiner

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
#  include <math.h>
#else
#  error Need ANSI-C headers
#endif

/* Auxilliary headers */
#include <complex.h>
#include "nfsft.h"
#include "util.h"
#include "../../nfft/utils.h"
#include "time.h"

#define H_MIN 0.6
#define H_MAX 0.9
#define H_STRIDE 0.2

#define M_MIN 32
#define M_MAX 512
#define M_STRIDE 32

#define L 128

#define D 1024

#define QUADRATIC_KERNEL(k) k==0?1.0:((2*k+1)/((double)(k*(k+1))*(k*(k+1))))
#define ABEL_POISSON_KERNEL(k,h) (2*k+1)*pow(h,k)*((1-h)*(1-h)/(1+h))
#define GAUSS_WEIERSTRASS(k,p) exp(-k*(k+1)*p)*(2*k+1)/(4*PI) 


inline double ip(double phi1,double theta1,double phi2,double theta2)
{
  return cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
}

inline double poisson(double phi1,double theta1,double phi2,double theta2,double h)
{
  double t = ip(phi1,theta1,phi2,theta2);
 	return pow((1-h)/sqrt(1-2*h*t+h*h),3);
}

/**
 * The main program.
 *
 * \param argc The number of arguments
 * \param argv An array containing the arguments as C-strings
 *
 * \return Exit code 
 */
int main (int argc, char **argv)
{  
  double h_min;
  double h_max;
  double h_stride;
  int m_min;
  int m_max;
  int m_stride;
  int l;
  int d;
  int n_max;

 	complex *b;
  complex **f_hat;
 	complex *a;
 	double *xi;
 	double *nu;
 	complex *f;
  complex *fd;
 	nfsft_plan plan, plan_adjoint;
 	int m,n,j,k;
  double h;
	
  if (argc == 1)
  {
    h_min = H_MIN;
    h_max = H_MAX; 
    h_stride = H_STRIDE; 
    m_min = M_MIN; 
    m_max = M_MAX; 
    m_stride = M_STRIDE; 
    l = L;
    d = D;
  }
  else if (argc == 9)
  {
    sscanf(argv[1],"%lf",&h_min);
    sscanf(argv[2],"%lf",&h_max);
    sscanf(argv[3],"%lf",&h_stride);
    sscanf(argv[4],"%d",&m_min);
    sscanf(argv[5],"%d",&m_max);
    sscanf(argv[6],"%d",&m_stride);
    sscanf(argv[7],"%d",&l);
    sscanf(argv[8],"%d",&d);
  }
  else
  {
    fprintf(stderr,"Convolution - Convolution test for NFSFT\n");
    fprintf(stderr,"Usage: convolution H_MIN H_MAX H_STRIDE M_MIN M_MAX M_STRIDE L D\n");
    return -1;
  }  

  fprintf(stderr,"%lf\n",h_min);
  fprintf(stderr,"%lf\n",h_max);
  fprintf(stderr,"%lf\n",h_stride);
  fprintf(stderr,"%d\n",m_min);
  fprintf(stderr,"%d\n",m_max);
  fprintf(stderr,"%d\n",m_stride);
  fprintf(stderr,"%d\n",l);
  fprintf(stderr,"%d\n",d);
  //return 0;
  
  n_max = 1<<ngpt(m_max);
  
  srand48(time(NULL));
  
  /** Allocate data structures. */
  b = (complex*) malloc(l*sizeof(complex));
  nu = (double*) malloc(2*l*sizeof(double));
  f_hat = (complex**) malloc((2*m_max+1)*sizeof(complex*));
  for (n = -m_max; n <= m_max; n++)
  {
    f_hat[n+m_max] = (complex*) malloc((n_max+1)*sizeof(complex));
  }  
  a = (complex*) malloc((m_max+1)*sizeof(complex));
  xi = (double*) malloc(2*d*sizeof(double));
  f = (complex*) malloc(d*sizeof(complex));
  fd = (complex*) malloc(d*sizeof(complex));
  
  nfsft_precompute(m_max,2000,0U);
  
  /* Create kernels. */
  fprintf(stderr,"kernels:\n");
  for (k = 0; k < l; k++)
  {
    b[k] = 1.0;//drand48();
    nu[2*k] = 0.0;//drand48()-0.5;
    nu[2*k+1] = 0.25;//0.5*drand48();
    fprintf(stderr,"k = %d, b[%d] = %5f, nu[2*%d] = %5f, nu[2*%d+1] = %5f\n",k,k,creal(b[k]),k,nu[2*k],k,nu[2*k+1]);
  }

  /* Create nodes. */
  fprintf(stderr,"nodes:\n");
  for (j = 0; j < d; j++)
  {
    xi[2*j] = ((double)j)/d-0.5;//drand48()-0.5;
    xi[2*j+1] = 0.25;//0.5*drand48();
    fprintf(stderr,"j = %d, xi[2*%d] = %5f, xi[2*%d+1] = %5f\n",j,j,xi[2*j],j,xi[2*j+1]);
  }
  
  for (h = h_min; h <= h_max; h = h + h_stride)
  {  
    fprintf(stderr,"h = %lf\n",h);

    /* Kernel coeffcients up to M */
    fprintf(stderr,"Kernels coefs:\n");
    for (k = 0; k <= m_max; k++)
    {
      a[k] = (2*k+1)*pow(h,k)*((1-h)*(1-h)/(1+h));
    	fprintf(stderr,"h = %f, k = %d, m_max = %d, a[%d] = %.16f\n",h,k,m_max,k,a[k]);
    }

    /* Compute real values of Poisson kernel */
    for (j = 0; j < d; j++)
    {
      fd[j] = 0.0;
      for (k = 0; k < l; k++)
      {
        fd[j] = fd[j] + b[k]*poisson(2*PI*nu[2*k],2*PI*nu[2*k+1],2*PI*xi[2*j],2*PI*xi[2*j+1],h);
        //fprintf(stderr,"%5lf %5lf %5lf %5lf %5lf %lf\n",creal(b[k]),nu[2*k],nu[2*k+1],xi[2*j],xi[2*j+1],creal(fd[j]));
      }
    }

    for (m = m_min; m <= m_max; m++)
    {      
      /* Target nodes */
      for (k = 0; k < d; k++)
      {
        xi[2*k] = drand48()-0.5;
        xi[2*k+1] = 0.5*drand48();        
      }  
      
      /* Adjoint transform */
      plan_adjoint = nfsft_init(l,m,nu,f_hat,b,0U);
      nfsft_adjoint(plan_adjoint);
      nfsft_finalize(plan_adjoint);
      
      /* Multiplication with diagonal matrix. */
      for (k = 0; k <= m; k++)
      {
        for (n = -k; n <= k; n++)
        {
          f_hat[n+m][k] *= a[k];
        }
      }
      
      /* Forward transform */
      plan = nfsft_init(d,m,xi,f_hat,f,0U);
      nfsft_trafo(plan);
      nfsft_finalize(plan);
      
      printf("%d\n",d);
      for (k = 0; k < d; k++)
      {
        printf("%+5f %+5f %16.15f %16.15f\n",xi[2*k],xi[2*k+1],creal(f[k]), creal(fd[k]));
      }  
    }
  }
  
  free(fd);
  free(f);
  free(xi);
  free(a);
  for (n = -m_max; n <= m_max; n++)
  {
    free(f_hat[n+m_max]);
  }   
  free(f_hat);
  free(b);
	
  return EXIT_SUCCESS;
}
