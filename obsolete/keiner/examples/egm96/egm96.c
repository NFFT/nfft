/* 
   egm96 - Fast spherical Fourier transform example for the NFSFT

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

#define M 360
#define N 512
#define D_PHI 360
#define D_THETA 180

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
 	complex *f;
  complex **f_hat;
  double *x;
  double *phi;
  double *theta;
 	nfsft_plan plan;
 	int n,k,d,j,l;
  FILE *file;
  double a,b,c,e;
  const int D = D_PHI * D_THETA;
  
  file = fopen("egm96.dat","r");
  if (file != NULL)
  {
    /** Allocate data structures. */
    f = (complex*) calloc(D,sizeof(complex));
    phi = (double*) calloc(D_PHI,sizeof(double));
    theta = (double*) calloc(D_THETA,sizeof(double));
    x = (double*) calloc(2*D,sizeof(double));
    f_hat = (complex**) malloc((2*M+1)*sizeof(complex*));
    for (n = -M; n <= M; n++)
    {
      f_hat[n+M] = (complex*) calloc((N+1),sizeof(complex));
    }  

    while (fscanf(file,"%d %d %le %le %le %le\n",&k,&n,&a,&b,&c,&e) != EOF)
    {
		  if (k <= M && k > 5)
			{
  		  a *= sqrt(2.0*k+1.0); 
	  	  b *= sqrt(2.0*k+1.0); 
        f_hat[n+M][k] = 0.5*(a + I*b);
        f_hat[-n+M][k] = 0.5*(a - I*b);
        //fprintf(stderr,"f_hat[%d][%d] = %E; a = %E\n",n,k,f_hat[n+M][k],creal(f_hat[n+M][k]));
			}
    }
    //f_hat[M][0] = 1.0;
    fclose(file);

    /*for (n = -M; n <= M; n++)
    {
      for (k = 3; k <= M ; k++)
      {
        f_hat[n+M][k] = 0.0;
        fprintf(stderr,"f[%d][%d] = %E\n",n,k,creal(f_hat[n+M][k]));
      }
    }*/  
    //f_hat[4+M][4] = 1E-7; 
    
    for (j = 0; j < D_PHI; j++)
    {
      phi[j] = ((double)j)/D_PHI - 0.5;
      fprintf(stderr,"phi[%d] = %f\n",j,phi[j]);
    }

    for (l = 0; l < D_THETA; l++)
    {
      theta[l] = 0.5*((double)l)/(D_THETA-1);
      fprintf(stderr,"theta[%d] = %f\n",l,theta[l]);
    }
    
    d = 0;
    for (j = 0; j < D_PHI; j++)
    {
      for (l = 0; l < D_THETA; l++)
      {
        x[2*d] = phi[j];
        x[2*d+1] = theta[l];
        d++;
      }
    }
    
    nfsft_precompute(M,2000);

    /* Forward transform */
    plan = nfsft_init(D,M,x,f_hat,f,0U);
    nfsft_trafo(plan);
    nfsft_finalize(plan);
    
    file = fopen("egm96f.dat","w");
    if (file != NULL)
    {
      for (d = 0; d < D; d++)
      {
        //fprintf(file,"%.16f %.16f %.16f %.16f\n",x[2*d],x[2*d+1],creal(f[d]),cimag(f[d]));
        fprintf(file,"%.16f\n",cabs(f[d]));
      }
      fclose(file);
    }

    for (n = -M; n <= M; n++)
    {
      free(f_hat[n+M]);
    }   
    free(f_hat);
    free(x);
    free(theta);
    free(phi);
    free(f);
  }  
  	
  return EXIT_SUCCESS;
}
