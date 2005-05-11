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
  const int M = 256;
 	const int L = 5;
  const int D_PHI = 40;
  const int D_THETA = D_PHI/2+1;
 	const int D = D_PHI * D_THETA;
  /** Next greater power of two with respect to L */
  const int N = 1<<((int)ceil(log((double)M)/log(2.0)));
 	complex *b;
  complex **f_hat;
 	complex *a;
 	double *xi;
 	double *nu;
 	complex *f;
 	nfsft_plan plan, plan_adjoint;
 	int k,n,d,d_theta,d_phi;
  double h = 0.94;
	
  srand48(time(NULL));
  
  /** Allocate data structures. */
  b = (complex*) malloc(L*sizeof(complex));
  nu = (double*) malloc(2*L*sizeof(double));
  f_hat = (complex**) malloc((2*M+1)*sizeof(complex*));
  for (n = -M; n <= M; n++)
  {
    f_hat[n+M] = (complex*) malloc((N+1)*sizeof(complex));
  }  
  a = (complex*) malloc((M+1)*sizeof(double));
  xi = (double*) malloc(2*D*sizeof(double));
  f = (complex*) malloc(D*sizeof(complex));
  
  for (k = 0; k < L; k++)
  {
    b[k] = 1.0;
    nu[2*k] = drand48()-0.5;
    nu[2*k+1] = 0.5*drand48();
  }
  
  /* Kernel coeffcients up to M */
  for (k = 0; k <= M; k++)
  {
    a[k] = pow(h,k)*(2*k+1);
  }
  
  /* Target nodes */
  d = 0;
  for (d_theta = 0; d_theta < D_THETA; d_theta++)
  {
    for (d_phi = 0; d_phi < D_PHI; d_phi++)
    {
      xi[2*d] = ((double)d_phi)/(D_PHI-1)-0.5;
      xi[2*d+1] = (0.5*d_theta)/(D_THETA-1);
      //printf("(%f,%f)\n",xi[2*d],xi[2*d+1]);
      d++;
    }
  }  
  
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
  
  printf("%d\n",D);
  printf("%d\n",D_PHI);
  printf("%d\n",D_THETA);
  for (d = 0; d < D; d++)
  {
    printf("%.16E %.16E %.16E\n",xi[2*d],xi[2*d+1],f[d]/300.0);
  }  
  
  free(f);
  free(xi);
  free(a);
  for (n = -M; n <= M; n++)
  {
    free(f_hat[n+M]);
  }   
  free(f_hat);
  free(b);
	
  return EXIT_SUCCESS;
}
