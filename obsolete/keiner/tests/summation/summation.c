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

#define M_MIN 32
#define M_MAX 256
#define M_STRIDE 32

#define H_MIN 0.7
#define H_MAX 0.9
#define H_STRIDE 0.2

#define D_MIN 1000
#define D_MAX 10000
#define D_STRIDE 1000

#define DD_MAX 40000

#define L_MIN 1000
#define L_MAX 10000
#define L_STRIDE 1000

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
  int l_min;
  int l_max;
  int l_stride;
  int d_min;
  int d_max;
  int d_stride;

  /** Next greater power of two with respect to M_MAX */
  int N_MAX;
  /* Bandwidth */
  int M;
  /** Next greater power of two with respect to M */
  int N;
	/** Error */
	double eps;
  double h;
	double err;
	double t_d, t_f;
	
 	complex *b;
  complex **f_hat;
 	complex *a;
 	double *xi;
 	double *nu;
 	complex *f, *f2;
 	nfsft_plan plan, plan_adjoint;
 	int j,k,n,l,d;
	FILE *file;
	 	
  if (argc == 1)
  {
    h_min = H_MIN;
    h_max = H_MAX; 
    h_stride = H_STRIDE;
    m_min = M_MIN; 
    m_max = M_MAX; 
    m_stride = M_STRIDE; 
    l_min = L_MIN; 
    l_max = L_MAX; 
    l_stride = L_STRIDE; 
    d_min = D_MIN; 
    d_max = D_MAX; 
    d_stride = D_STRIDE; 
  }
  else if (argc == 13)
  {
    sscanf(argv[1],"%lf",&h_min);
    sscanf(argv[2],"%lf",&h_max);
    sscanf(argv[3],"%lf",&h_stride);
    sscanf(argv[4],"%d",&m_min);
    sscanf(argv[5],"%d",&m_max);
    sscanf(argv[6],"%d",&m_stride);
    sscanf(argv[7],"%d",&l_min);
    sscanf(argv[8],"%d",&l_max);
    sscanf(argv[9],"%d",&l_stride);
    sscanf(argv[10],"%d",&d_min);
    sscanf(argv[11],"%d",&d_max);
    sscanf(argv[12],"%d",&d_stride);
  }
  else
  {
    fprintf(stderr,"Convolution - Convolution test for NFSFT\n");
    fprintf(stderr,"Usage: convolution H_MIN H_MAX H_STRIDE M_MIN M_MAX M_STRIDE L_MIN L_MAX L_STRIDE D_MIN D_MAX D_STRIDE\n");
    return -1;
  }  
		
  printf("d = %d\n",d);
	N_MAX = 1<<ngpt(m_max);	
		
	nfsft_precompute(m_max,2000);
	
  /** Allocate data structures. */
  b = (complex*) malloc(l_max*sizeof(complex));
  nu = (double*) malloc(2*l_max*sizeof(double));
  f_hat = (complex**) malloc((2*m_max+1)*sizeof(complex*));
  for (n = -m_max; n <= m_max; n++)
  {
    f_hat[n+m_max] = (complex*) malloc((N_MAX+1)*sizeof(complex));
  }  
  a = (complex*) malloc((m_max+1)*sizeof(complex));
  xi = (double*) malloc(2*d_max*sizeof(double));
  f = (complex*) malloc(d_max*sizeof(complex));
  f2 = (complex*) malloc(d_max*sizeof(complex));
	  
  /* Target nodes */
  j = 0;
  for (j = 0; j < d_max; j++)
  {
    xi[2*j] = drand48();
    xi[2*j+1] = 0.5*drand48();
  }

  /* Init kernels */
	/* Create kernels. */
	for (k = 0; k < l_max; k++)
	{
		b[k] = drand48();
		nu[2*k] = drand48()-0.5;
		nu[2*k+1] = 0.5*drand48();
	}

  file = fopen("summation.tex","w");
  fprintf(file,"\\begin{tabular}{l|l|l|l|l|l|l}\n");
	fprintf(file,"$h$ & $M$ & $K$ & $N$ & $t_{\\text{fast}}$ & $t_{\\text{slow}}$ & $\\text{err}_{\\infty}$\\\\\\hline\n");

  for (l = l_min; l <= l_max; l = l + l_stride)
	{
	  for (d = d_min; d <= d_max; d = d + d_stride)
		{
  		int temp = l;
		  l = d/2;
			for (h = h_min; h <= h_max; h = h + h_stride)
			{
				/* Kernel coeffcients up to m_max */
				for (k = 0; k <= m_max; k++)
				{
					a[k] = ABEL_POISSON_KERNEL(k,h);
				}
				
				if (d <= DD_MAX)
				{
					t_d = mysecond();
					for (j = 0; j < d; j++)
					{
						f2[j] = 0.0;
						for (k = 0; k < l; k++)
						{
							f2[j] += b[k]*poisson(2*PI*nu[2*k],2*PI*nu[2*k+1],2*PI*xi[2*j],2*PI*xi[2*j+1],h);
						}
					} 
					t_d = mysecond() - t_d;
				}
				
				for (M = m_min; M <= m_max; M = M + m_stride)
				{
					fprintf(stderr,"L = %d, D = %d, h = %lf, M = %d\n",l,d,h,M);
					fflush(stderr);
					/* Adjoint transform */
					plan_adjoint = nfsft_init(l,M,nu,f_hat,b,0U);
					plan = nfsft_init(d,M,xi,f_hat,f,0U);
					t_f = mysecond();
					nfsft_adjoint(plan_adjoint);
				
					/* Multiplication with diagonal matrix. */
					for (k = 0; k <= M; k++)
					{
						for (n = -k; n <= k; n++)
						{
							f_hat[n+M][k] *= a[k];
						}
					}
				
					/* Forward transform */
					nfsft_trafo(plan);
					t_f = mysecond() - t_f;
					nfsft_finalize(plan_adjoint);
					nfsft_finalize(plan);
					
					if (d <= DD_MAX)
					{
					  fprintf(file,"%.2f & %3d & %6d & %6d & %5.2f & %5.2f & %.4E\\\\\n",h,M,l,d,t_d,t_f,error_complex_inf(f, f2, d));
					}
					else
					{
					  fprintf(file,"%.2f & %3d & %6d & %6d & -- & %5.2f & --\\\\\n",h,M,l,d,t_f);
					}
				}
			}
			l = temp;
		}
		fprintf(file,"\\hline");
	}
		
	fprintf(file,"\\end{tabular}\n");
  fclose(file);
	
  free(f);
	free(f2);
  free(xi);
	free(nu);
  free(a);
  for (n = -m_max; n <= m_max; n++)
  {
    free(f_hat[n+m_max]);
  }   
  free(f_hat);
  free(b);
	
  return EXIT_SUCCESS;
}
