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
  int M;
  int N;
  int n;
  int k;
  complex *coefs;
  complex *chebc;
  double ctime;
  int re,im;
  struct nfsft_wisdom *wisdom;
  
  fscanf(stdin,"%d\n",&M);
  fscanf(stdin,"%d\n",&n);
  
  N = 1<<ngpt(M);
  
  /* Initialize nodes. */
  coefs = (complex*) fftw_malloc((N+1)*sizeof(complex));
  for (k = abs(n); k <= M ; k++)
  {   
    fscanf(stdin,"%lf\n%lf\n",&re,&im);
    coefs[k] = re + I*im;
    //fprintf(stderr,"%lf + %lfi\n",creal(coefs[k]),cimag(coefs[k]));
  } 
  
  
		/* Initialize */
 	nfsft_precompute_guru(N,THRESHOLD,n);
    
  wisdom = nfsft_get_wisdom();
  
  flft(M,ngpt(M),n,coefs,wisdom);
    
  for (k = abs(n); k <= M ; k++)
  {   
    fprintf(stdout,"%.30lf\n%.30lf\n",creal(coefs[k]),cimag(coefs[k]));
  } 

  /* Forget wisdom. */
  nfsft_forget_guru(n);
  
  fftw_free(coefs);
    
  return EXIT_SUCCESS;
}
