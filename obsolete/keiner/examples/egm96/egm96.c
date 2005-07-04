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

#define D_PHI 360
#define D_THETA 180
#define K_MIN 3
#define K_MAX 360

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
  complex **f_hat_orig;
  double *x;
  double *phi;
  double *theta;
  nfsft_plan plan;
  int n,m,k,d,j,l,i, k_min, k_max, d_phi, d_theta, N;
  FILE *file;
  char filename[100];
  double a,b,c,e;
  int D;
  
  if (argc == 1)
  {
    k_min = K_MIN;
    k_max = K_MAX; 
    d_phi = D_PHI;
    d_theta = D_THETA; 
  }
  else if (argc == 5)
  {
    sscanf(argv[1],"%d",&k_min);
    sscanf(argv[2],"%d",&k_max);
    sscanf(argv[3],"%d",&d_phi);
    sscanf(argv[4],"%d",&d_theta);
  }
  else
  {
    fprintf(stderr,"EGM96 - EGM96 example for NFSFT\n");
    fprintf(stderr,"Usage: EGM96 K_MIN K_MAX D_PHI D_THETA\n");
    return -1;
  }    
  
  D = d_phi * d_theta;
  N = 1<<ngpt(k_max);
  
  /** Allocate data structures. */
  f = (complex*) calloc(D,sizeof(complex));
  phi = (double*) calloc(d_phi,sizeof(double));
  theta = (double*) calloc(d_theta,sizeof(double));
  x = (double*) calloc(2*D,sizeof(double));
  f_hat = (complex**) malloc((2*k_max+1)*sizeof(complex*));
  f_hat_orig = (complex**) malloc((2*k_max+1)*sizeof(complex*));
  for (n = -k_max; n <= k_max; n++)
  {
    f_hat[n+k_max] = (complex*) calloc((N+1),sizeof(complex));
    f_hat_orig[n+k_max] = (complex*) calloc((N+1),sizeof(complex));
  }  

  for (j = 0; j < d_phi; j++)
  {
    phi[j] = ((double)j)/d_phi - 0.5;
  }

  for (l = 0; l < d_theta; l++)
  {
    theta[l] = 0.5*((double)l)/(d_theta-1);
  }
  
  d = 0;
  for (j = 0; j < d_phi; j++)
  {
    for (l = 0; l < d_theta; l++)
    {
      x[2*d] = phi[j];
      x[2*d+1] = theta[l];
      d++;
    }
  }

  file = fopen("egm96.dat","r");
  while (fscanf(file,"%d %d %le %le %le %le\n",&j,&m,&a,&b,&c,&e) != EOF)
  {
    if (j <= k_max && j >2)
    {
      a *= sqrt(2.0*j+1.0); 
      b *= sqrt(2.0*j+1.0); 
      f_hat_orig[m+k_max][j] = 0.5*(a + I*b);
      f_hat_orig[-m+k_max][j] = 0.5*(a - I*b);
    }
  }
  fclose(file);
  
  nfsft_precompute(k_max,2000,0U);

  plan = nfsft_init(k_max,D,f_hat,x,f,0U);
  fprintf(stderr,"d_phi = %d, d_theta = %d, D = %d\n",d_phi,d_theta,D);
  i = 0;
  for (k = k_min; k <= k_max; k++)
  {  
    for (n = k; n <= k; n++)
    {  
      fprintf(stderr,"k = %d, n = %d\n",k,n);
      for (m = -k_max; m <= k_max; m++)
      {  
        if (abs(m) <= n)
        {
          memcpy(f_hat[m+k_max],f_hat_orig[m+k_max],(k+1)*sizeof(complex));
          memset(&(f_hat[m+k_max][k+1]),0U,(N-k)*sizeof(complex));
        }
        else
        {
          memset(f_hat[m+k_max],0U,(N+1)*sizeof(complex));
        }
      }

      /* Forward transform */
      nfsft_trafo(plan);
    
      sprintf(filename,"egm96%05d.dat",k-3);
      file = fopen(filename,"w");
      for (d = 0; d < D; d++)
      {
        fprintf(file,"%.16f\n",cabs(f[d]));
      }
      fclose(file);
      i++;
    }
  }  
  nfsft_finalize(plan);
    
  for (n = -k_max; n <= k_max; n++)
  {
    free(f_hat[n+k_max]);
    free(f_hat_orig[n+k_max]);
  }   
  free(f_hat);
  free(f_hat_orig);
  free(x);
  free(theta);
  free(phi);
  free(f);
  
  return EXIT_SUCCESS;
}
