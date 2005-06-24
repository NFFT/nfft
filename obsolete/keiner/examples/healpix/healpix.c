/* 
   Healpix - Healpix example for the NFSFT

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
#include "api.h"

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
  complex *f_orig;
  double *x;
  nfsft_plan plan;
  int n,m,k,d,j,l,i, k_min, k_max, d_phi, d_theta, N, n_d;
  FILE *file;
  char filename[100];
  double a,b,c,e;
  int D;
  infsft_plan iplan;
  
  if (argc == 1)
  {
    k_min = K_MIN;
    k_max = K_MAX; 
  }
  else if (argc == 3)
  {
    sscanf(argv[1],"%d",&k_min);
    sscanf(argv[2],"%d",&k_max);
  }
  else
  {
    fprintf(stderr,"Healpix - Healpix example for NFSFT\n");
    fprintf(stderr,"Usage: EGM96 K_MIN K_MAX\n");
    return -1;
  }    
  
  /* Read Healpix nodes. */  
  sprintf(filename,"healpix.dat");  
  file = fopen(filename,"r");
  
  if (file != NULL)
  {
    fscanf(file,"%d\n",&n_d);
    
    D = 12*n_d *n_d;
    N = 1<<ngpt(k_max);

    /** Allocate data structures. */
    f = (complex*) calloc(D,sizeof(complex));
    f_orig = (complex*) calloc(D,sizeof(complex));
    x = (double*) calloc(2*D,sizeof(double));
    f_hat = (complex**) malloc((2*k_max+1)*sizeof(complex*));
    f_hat_orig = (complex**) malloc((2*k_max+1)*sizeof(complex*));
    for (n = -k_max; n <= k_max; n++)
    {
      f_hat[n+k_max] = (complex*) calloc((N+1),sizeof(complex));
      f_hat_orig[n+k_max] = (complex*) calloc((N+1),sizeof(complex));
    }  
    
    /* Read Healpix nodes */
    for (d = 0; d < D; d++)
    {
  				fscanf(file,"%lf %lf\n",&x[2*d+1],&x[2*d]);
      //fprintf(stderr,"%e %e\n",x[2*d+1],x[2*d]);
    }
    
    fclose(file);

    /* Read EGM96 data. */
    file = fopen("../egm96/egm96.dat","r");
    /*while (fscanf(file,"%d %d %le %le %le %le\n",&j,&m,&a,&b,&c,&e) != EOF)
    {
      if (j <= k_max && j >2)
      {
        a *= sqrt(2.0*j+1.0); 
        b *= sqrt
          (2.0*j+1.0); 
        f_hat_orig[m+k_max][j] = 0.5*(a + I*b);
        f_hat_orig[-m+k_max][j] = 0.5*(a - I*b);
      }
    }*/
    fclose(file);
    for (k = 0; k <= k_max; k++)
    {
      for (n = -k; n <= k ; n++)
      {
        f_hat_orig[n+k_max][k] = drand48() + I * drand48();
      }  
    }  
    
    copyc_hat(f_hat,f_hat_orig,k_max);

    nfsft_precompute(k_max,2000,0U);
    
    plan = nfsft_init(D,k_max,x,f_hat,f,0U);

    nfsft_trafo(plan);
    
    copyc(f_orig,f,D);

    /* Init guess */
    for (n = -k_max; n <= k_max; n++)
    {
      memset(f_hat[n+k_max],0U,(N+1)*sizeof(complex));
    }  
    
    iplan = infsft_make_plan();
    infsft_init_guru(iplan, plan, CGNR_E);
    
    /* inverse trafo */  
    infsft_before_loop(iplan);
    fprintf(stderr,"%e,\n",sqrt(iplan->dot_r_iter));
    for(l=0;l<50;l++)
    { 
      fprintf(stderr,"%e,\n",sqrt(iplan->dot_r_iter));
      infsft_loop_one_step(iplan);
    }
    
    infsft_finalize(iplan);  
    nfsft_finalize(plan);  
    
    nfsft_forget();
    
    for (n = -k_max; n <= k_max; n++)
    {
      free(f_hat[n+k_max]);
      free(f_hat_orig[n+k_max]);
    }   
    free(f_hat);
    free(f_hat_orig);
    free(x);
    free(f);
    free(f_orig);
  }      
  
  return EXIT_SUCCESS;
}
