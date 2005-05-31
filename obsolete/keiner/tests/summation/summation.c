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
#include "time.h"

#define M_MIN 64
#define M_STRIDE 4
#define M_MAX 64

#define EPS_START 1E-1
#define EPS_STRIDE 5E-1
#define EPS_END 1E-6 

#define M_ERR_QUADRATIC(eps) (int)ceil(0.5*(sqrt(1+4/eps)-1))
#define M_ERR_ABEL_POISSON(eps) ()

inline double ip(double phi1,double theta1,double phi2,double theta2)
{
  //printf("theta = %f -> cos(theta) = %f\n",theta1,cos(theta1));
  return cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
}

inline double poisson(double phi1,double theta1,double phi2,double theta2,double h)
{
  double t = ip(phi1,theta1,phi2,theta2);
	return pow((1-h)/sqrt(1-2*h*t+h*h),3);
}

inline double gauss_weierstrass(double phi1,double theta1,double phi2,double theta2,double h)
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
  /** Next greater power of two with respect to M_MAX */
  const int N_MAX = 1<<ngpt(M_MAX);
  /* Bandwidth */
  int M;
	/* Use only one kernel */
 	const int L = 1;
	/* Number of nodes in colatitudinal direction */
  const int D_PHI = 120;
	/* Use only one theta */
  const int D_THETA = 1;
	/* Total number of nodes */
 	const int D = D_PHI * D_THETA;
  /** Next greater power of two with respect to M */
  int N;
	/** Error */
	double eps;
  double h;
	double err;
	
 	complex *b;
  complex **f_hat;
 	complex *a;
 	double *xi;
 	double *nu;
 	complex *f, *f2;
 	nfsft_plan plan, plan_adjoint;
 	int k,n,d,d_theta,d_phi;
	 
	/** Test quadratic kernel */
	/*for (eps = EPS_START; eps >= EPS_END; eps = eps * EPS_STRIDE)
	{
	  printf("eps = %.4E => M >= %d\n",eps,M_ERR_QUADRATIC(eps));
	}*/
	
  /** Allocate data structures. */
  b = (complex*) malloc(L*sizeof(complex));
  nu = (double*) malloc(2*L*sizeof(double));
  f_hat = (complex**) malloc((2*M_MAX+1)*sizeof(complex*));
  for (n = -M_MAX; n <= M_MAX; n++)
  {
    f_hat[n+M_MAX] = (complex*) malloc((N_MAX+1)*sizeof(complex));
  }  
  a = (complex*) malloc((M_MAX+1)*sizeof(complex));
  xi = (double*) malloc(2*D*sizeof(double));
  f = (complex*) malloc(D*sizeof(complex));
  f2 = (complex*) malloc(D*sizeof(complex));
	  
	/* One kernel only */
  b[0] = 1.0;
  nu[0] = 0.0;
  nu[1] = 0.25;;
  
	h = 0.90;
	
  /* Kernel coeffcients up to M_MAX */
  for (k = 0; k <= M_MAX; k++)
  {
    //a[k] = QUADRATIC_KERNEL(k);
    //a[k] = ABEL_POISSON_KERNEL(k,h);
    a[k] = GAUSS_WEIERSTRASS(k,1.00);
		//printf("%f\n",GAUSS_WEIERSTRASS(k,10.0));
  }

  /* Target nodes */
  d = 0;
  for (d_phi = 0; d_phi < D_PHI; d_phi++)
  {
    xi[2*d] = ((double)d_phi)/(D_PHI-1)-0.5;
    xi[2*d+1] = 0.25;
    //printf("(%f,%f)\n",xi[2*d],xi[2*d+1]);
    d++;
  }
	
  for (d = 0; d < D; d++)
  {
    f2[d] = b[0]*poisson(2*PI*nu[0],2*PI*nu[1],2*PI*xi[2*d],2*PI*xi[2*d+1],h);
  } 
	
	for (M = M_MIN; M <= M_MAX; M = M + M_STRIDE)
	{
    /* Adjoint transform */
	  plan_adjoint = nfsft_init(L,M,nu,f_hat,b,0U);
	  nfsft_adjoint(plan_adjoint);
	  nfsft_finalize(plan_adjoint);
  
   	/* Multiplication with diagonal matrix. */
 	  for (k = 0; k <= M; k++)
   	{
 	    for (n = -k; n <= k; n++)
	   	{
		    f_hat[n+M][k] *= a[k];
		  }
    }
	
	  /* Forward transform */
	  plan = nfsft_init(D,M,xi,f_hat,f,0U);
	  nfsft_trafo(plan);
    nfsft_finalize(plan);

    /*printf("%d\n",D);
    printf("%d\n",D_PHI);
    printf("%d\n",D_THETA);*/
		
      //printf("%.16E %.16E %.16E %.16E",xi[2*d],xi[2*d+1],creal(f[d]),b[0]*poisson(2*PI*nu[0],2*PI*nu[1],2*PI*xi[2*d],2*PI*xi[2*d+1],h));
			//printf(" %E\n",fabs(creal(f[d])-b[0]*poisson(2*PI*nu[0],2*PI*nu[1],2*PI*xi[2*d],2*PI*xi[2*d+1],h)));
		  //err = 0.0;
		
		for (d = 0; d < D; d++)
		{
		  printf("%.4E\n",creal(f[d]));
		}
		 //printf("%d: err = %.4E\n",M,error_complex_inf(f,f2,D)); 
	}
	
	//printf("%f\n",ip(0.0,PI/2.0,0.0,PI/2.0));
		
  free(f);
	free(f2);
  free(xi);
	free(nu);
  free(a);
  for (n = -M_MAX; n <= M_MAX; n++)
  {
    free(f_hat[n+M_MAX]);
  }   
  free(f_hat);
  free(b);
	
  return EXIT_SUCCESS;
}
