/* 
   flft - FLFT test for the NFSFT

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
#include "api.h"
#include "util.h"

#define THRESHOLD 1000

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
  /** The bandwidth */
  int M = 4;
  int N = 1<<ngpt(M);
  int n = 0;
  int k;
  complex *leg;
  
  leg = (complex*) fftw_malloc((N+1)*sizeof(complex));

		/* Initialize */
 	nfsft_precompute(N,1E6,0U);

  for (k = 0; k <= M ; k++)
  {   
    //leg[k] = (k<abs(n))?(0.0):(drand48() + I*drand48());
    leg[k] = (k == M-1)?(1.0):(0.0);
    fprintf(stderr,"leg[%d] = %lf + %lfi\n",k,creal(leg[k]),cimag(leg[k]));
  }   
  fprintf(stderr,"\n");
 
  nfsft_flft_trafo(leg,n,M);
    
  for (k = 0; k <= M ; k++)
  {   
    fprintf(stderr,"leg[%d] = %lf + %lfi\n",k,creal(leg[k]),cimag(leg[k]));
  } 
  fprintf(stderr,"\n");

  nfsft_flft_trafo_adjoint(leg,n,M);

  for (k = 0; k <= M ; k++)
  {   
    fprintf(stderr,"leg[%d] = %lf + %lfi\n",k,creal(leg[k]),cimag(leg[k]));
  } 
  fprintf(stderr,"\n");

  nfsft_forget();
  
  fftw_free(leg);
    
  return EXIT_SUCCESS;
}
