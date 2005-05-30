/* 
   inverse - Inverse test for the NFSFT

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
#include <time.h>
#include <fftw3.h>
#include "nfsft.h"
#include "util.h"
#include "api.h"
#include "infsft.h"

#define M_DEFAULT 16

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
  /** The current bandwidth */
  int M;
  /** Next greater power of two with respect to M */
  int N;
  /** Maximum number of nodes */
  int D;
  /** Loop counter for Legendre index k. */
  int k;
  /** Loop counter for Legendre index n. */
  int n;
  /** Loop counter for nodes. */
  int d;  
  int l;
  
  /** Arrays for complex Fourier coefficients. */
  complex **f_hat;
  /** Copy the original Fourier coefficients. */
  complex **f_hat_orig;
  
  /** Array of angles defining the nodes. */
  double *angles;
  /** Array for angles \f$\theta\f$ of a grid. */
  double *theta;
  /** Array for angles \f$\phi\f$ of a grid. */
  double* phi;
  /** Array for Gauss-Legendre weights. */
  //double *w;  
  /** Array for function values. */
  complex *f;
  double err_infty, err_1, err_2;
  
  /** Plan for fast spherical fourier transform. */
  nfsft_plan plan;
  infsft_plan iplan;
  
  /** Used to measure computation time. */
  double ctime;
  /** 
    * Used to store the filename of a file containing Gauss-Legendre nodes and 
    * weights.
    */ 
  char filename[100];
  unsigned short seed[3]={1,2,3};
  /** File handle for reading quadrature nodes and weights. */
  FILE *file;
  /** FFTW plan for Clenshaw-Curtis quadrature */
  fftw_plan fplan;

  
  if (argc == 1)
  {
    M = M_DEFAULT;
  }
  else if (argc == 2)
  {
    sscanf(argv[1],"%d",&M);
  }
  else
  {
    fprintf(stderr,"Inverse - Inverse test for NFSFT\n");
    fprintf(stderr,"Usage: inverse M\n");
    return -1;
  }  
  
  N = 1<<ngpt(M);
  D = (2*M+1)*(2*M+2);
  
  /* Initialize random number generator. */
  //srand48(time(NULL));
    
  /* Allocate memory. */
  f_hat = (complex**) malloc((2*M+1)*sizeof(complex*));
  f_hat_orig = (complex**) malloc((2*M+1)*sizeof(complex*));
  for (n = -M; n <= M; n++)
  {
    f_hat[n+M] = (complex*) malloc((N+1)*sizeof(complex));
    f_hat_orig[n+M] = (complex*) malloc((N+1)*sizeof(complex));
  }  
  
  angles = (double*) malloc(2*D*sizeof(double));
  theta = (double*) malloc((2*M+1)*sizeof(double));
  phi = (double*) malloc((2*M+2)*sizeof(double));
  f = (complex*) malloc(2*D*sizeof(complex));
    
  printf("Precomputing wisdom up to M = %d ...",M);
  fflush(stdout);
  ctime = mysecond();
  nfsft_precompute(M,1000);
  printf(" needed %7.2f secs.\n",mysecond()-ctime);
          
  /* Compute random Fourier coefficients. */
  for (n = -M; n <= M; n++)
  {
    for (k = abs(n); k <= M; k++)
    {
      f_hat[n+M][k] = drand48() + I*drand48();
    }
    /* Save a copy. */
    memcpy(f_hat_orig[n+M],f_hat[n+M],(M+1)*sizeof(complex));
  }
        
  /* Compute Clenshaw-Curtis nodes and weights. */  
  for (k = 0; k < 2*M+2; k++)
  {
    phi[k] = k/((double)2*M+2);
  }
  for (n = 0; n < 2*M+1; n++)
  {
    theta[n] = n/((double)4*M);
  }
  /* Create grid nodes. */
  d = 0;
  for (n = 0; n < 2*M+1; n++)
  {
    for (k = 0; k < 2*M+2; k++)
    {
      angles[2*d] = phi[k];
      angles[2*d+1] = theta[n];
      d++;
    }  
  }
    
  plan = nfsft_init(D, M, angles, f_hat, f, 0U);
  nfsft_trafo(plan);

  for (n = -M; n <= M; n++)
  {
    for (k = abs(n); k <= M; k++)
    {
      f_hat[n+M][k] = 0.0;
    }
  }

  iplan = infsft_make_plan();
  infsft_init(iplan, plan);
  
  /* inverse trafo */  
  infsft_before_loop(iplan);
  for(l=0;l<30;l++)
  { 
    fprintf(stderr,"%e,\n",sqrt(iplan->dot_r_iter));
    infsft_loop_one_step(iplan);
  }
    
  infsft_finalize(iplan);  
  nfsft_finalize(plan);  
  
  nfsft_forget();
  for (n = -M; n <= M; n++)
  {
    free(f_hat[n+M]);
    free(f_hat_orig[n+M]);
  }  
  free(f_hat);
  free(f_hat_orig);
  
  free(angles);
  free(theta);
  free(phi);
  free(f);
  return EXIT_SUCCESS;
}
