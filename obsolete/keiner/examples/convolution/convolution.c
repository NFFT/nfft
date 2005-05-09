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
  const int M = 128;
	const int L = 64;
	const int D = 64;
	complex **a_tilde;
	complex *b;
	double *nu;
	complex *a;
	double *xi;
	complex *f;
	nfsft_plan plan, plan_adjoint;
	int k,n;
	
	
	/* Adjoint transform */
	plan_adjoint = nfsft_init(L,M,nu,a_tilde,b,0U);
	nfsft_adjoint(plan_adjoint);
	
	/* Multiplication with diagonal matrix. */
	for (k = 0; k <= M; k++)
	{
	  for (n = -k; n <= k; n++)
		{
		  a_tilde[n+M][k] *= a[K];
		}
	}
	
	/* Forward transform */
	plan_adjoint = nfsft_init(D,M,xi,a_tilde,f,0U);
	nfsft_trafo(plan);
	
  return EXIT_SUCCESS;
}
