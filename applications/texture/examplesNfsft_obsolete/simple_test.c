/* 
   accuracy - Accuracy test for the NFSFT

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
#  include <string.h>
#else
#  error Need ANSI-C headers
#endif

/* Auxilliary headers */
#include <complex.h>
#include <time.h>
#include <fftw3.h>
#include "nfsft_old.h"
#include "../../../../lowlevel/nfsft/obsolete/util.h"

#define M 256
#define D 100
#define THRESHOLD 1000.0

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
  /** Arrays for complex Fourier coefficients. */
  complex **f_hat;
  /** Copy the original Fourier coefficients. */
  complex **f_hat_direct;
  /** Array of angles defining the nodes. */
  double *angles;
  /** Array for function values. */
  complex *f;
  /** Array for function values. */
  complex *f_direct;
  double err_infty, err_1, err_2;
  int N = 1<<ngpt(M);
  int k,n,d;
  
  /** Plan for fast spherical fourier transform. */
  nfsft_plan_old plan;
  nfsft_plan_old plan_direct;
  
  unsigned short seed[3]={1,2,3};

  /* Initialize random number generator. */
  seed48(seed);
    
  /* Allocate memory. */
  f_hat = (complex**) malloc((2*M+1)*sizeof(complex*));
  f_hat_direct = (complex**) malloc((2*M+1)*sizeof(complex*));
  for (n = -M; n <= M; n++)
  {
    f_hat[n+M] = (complex*) fftw_malloc((N+1)*sizeof(complex));
    f_hat_direct[n+M] = (complex*) fftw_malloc((N+1)*sizeof(complex));
  }  
  
  angles = (double*) malloc(2*D*sizeof(double));
  f = (complex*) malloc(D*sizeof(complex));
  f_direct = (complex*) malloc(D*sizeof(complex));

  nfsft_precompute_old(M, THRESHOLD, 0U);

  /* Compute random Fourier coefficients. */
  for (n = -M; n <= M; n++)
  {
    for (k = abs(n); k <= M; k++)
    {
      f_hat[n+M][k] = (drand48()-0.5) + I * (drand48()-0.5);
    }
    /* Save a copy. */
    memcpy(f_hat_direct[n+M],f_hat[n+M],(M+1)*sizeof(complex));
  }
  
  /* Random angles */
  for (d = 0; d < D; d++)
  {
    angles[2*d] = drand48()-0.5;
    angles[2*d+1] = 0.5*drand48();
  }  

  plan = nfsft_init_guru_old(M, D, f_hat, angles, f, NFSFT_NORMALIZED_OLD,6);
  plan_direct = nfsft_init_guru_old(M, D, f_hat_direct, angles, f_direct, NFSFT_NORMALIZED_OLD,6);

  //ndsft_trafo(plan);
  nfsft_trafo_old(plan);
  ndsft_trafo_old(plan_direct);
  fprintf(stdout,"reached!\n");
  fflush(stdout);
  
  for (d = 0; d < D; d++)
  {
    fprintf(stdout,"%.4E %.4E\n",cabs(f[d]),cabs(f_direct[d]));
  }

  err_infty = nfft_error_l_infty_complex(f,f_direct,D);
  err_2 = nfft_error_l_2_complex(f,f_direct,D);
            
  fprintf(stdout,"M = %d, D = %d: %.4E %.4E\n",M,D,err_infty, err_2);

  nfsft_finalize_old(plan);      
  nfsft_finalize_old(plan_direct);      
      
  nfsft_forget_old();
  for (n = -M; n <= M; n++)
  {
    free(f_hat[n+M]);
    free(f_hat_direct[n+M]);
  }  
  free(f_hat);
  free(f_hat_direct);
  
  free(angles);
  free(f);
  free(f_direct);
  return EXIT_SUCCESS;
}
