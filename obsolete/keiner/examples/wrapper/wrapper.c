/* 
   nfsft - Spherical Fourier transform for nonuniform sampled data

   Copyright (C) 2004 Jens Keiner

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
#  include <float.h>
#  include <math.h>
#else
#  error Need ANSI-C headers
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#else
#  error Need sys/types.h
#endif

/* Auxilliary headers */
#include <complex.h>
#include <fftw3.h>
#include "nfsft.h"
#include "util.h"

#define FORWARD 0
#define ADJOINT 1
#define FAST 0
#define SLOW 1 

/** Arrays for complex Fourier-coefficients. */
complex *f_hat;
complex *f_hat2;
complex **F_HAT;
complex **F_HAT2;

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
  /** Next greater power of 2 relative to M */
  int N;
  /** The number of nodes */
  int D; 
	/**Array of angles defining the nodes. */
  double *angles;
  /** Plan for fast spherical fourier transform. */
  nfsft_plan plan;
  complex *result;
  double ctime;
  int i,n,k;
  int mode, type;
    
	/* Read transform type. */
  fscanf(stdin,"%d\n",&type);
	/* Read transform mode. */	
  fscanf(stdin,"%d\n",&mode);

  /* Initialize nodes. */
  fscanf(stdin,"%d\n",&D);
  angles = (double*) calloc(2*D,sizeof(double));
  for (i = 0; i < D; i++)
  {
    fscanf(stdin,"%lf\n%lf\n",&angles[2*i],&angles[2*i+1]);    
  } 
#ifdef DEBUG
  myvpr(angles,2*D,"angles");
#endif

  /* Read bandwidth M. */  
  fscanf(stdin,"%d\n",&M);
#ifdef DEBUG
  printf("M = %d\n",M);
#endif
  /* Calculate N as next greater power of 2 of the bandwidth M. */
  N = M==0?0:1<<(int)ceil(log((double)M)/log(2.0));
#ifdef DEBUG
  printf("N = %d\n",N);
#endif
  
  /* Initialize data structures for Fourier coefficients. */
  F_HAT = (complex**) calloc(2*M+1,sizeof(complex*));
  for (n = -M; n <= M; n++)
  {
    f_hat = (complex*) calloc(N+1,sizeof(complex));
    F_HAT[n+M] = f_hat;
  }

  /* Initialize array for function values. */
  result = (complex*) calloc(D,sizeof(complex));
  
  if (type == FORWARD)
  {      
    /* Read Fourier coefficients from standard input. */    
    double c;
    for (k = 0; k <= M; k++)
    {
      for (n = -k; n <= k; n++)
      {
        fscanf(stdin,"%lf\n",&c);
        F_HAT[n+M][k] = c;
      }  
    }  

    /* Create plan for fast spherical Fourier transform.*/
    plan = nfsft_init(D, M, angles, F_HAT, result, 0U);
			
	  /* Switch by mode. */
		if (mode == FAST)
		{		
			/* Initialize */
			ctime = mysecond();
			nfsft_compute_wisdom(M,1000);
			ctime = mysecond() - ctime;
	    #ifdef DEBUG
			printf("Time for initialization: %f sec.\n",ctime);
	    #endif
			/* Precompute. */
			ctime = mysecond();
			nfsft_compute_wisdom(M,1000);
			ctime = mysecond() - ctime;
	    #ifdef DEBUG
			printf("Time for precomputation: %f sec.\n",ctime);
	    #endif
			/* Execute the plan. */
			ctime = mysecond();
			nfsft_trafo(plan);
			ctime = mysecond() - ctime;
	    #ifdef DEBUG
			printf("Time for fast algorithm: %f sec.\n",ctime);
	    #endif
	  }
	  else if (mode == SLOW)
		{
			/* Execute the plan. */
			ctime = mysecond();
			ndsft_trafo(plan);
			ctime = mysecond() - ctime;
      #ifdef DEBUG
			printf("Time for slow algorithm: %f sec.\n",ctime);
	    #endif
		}
		else
		{
		  fprintf(stderr,"Wrong transform mode!\n");
		}

   	for (k = 0; k < D; k++)
		{
			printf("%17.16f\n%17.16f\n",creal(result[k]),cimag(result[k]));
		}
  }
  else
  {
    /* Read function values from standard input. */    
    double c;
    for (k = 0; k < D; k++)
    {
      fscanf(stdin,"%lf\n",&c);
      result[k] = c;
    }  
    
    /* Create plan for fast spherical Fourier transform.*/
    plan = nfsft_init(D, M, angles, F_HAT, result, 0U);
		
		/* Switch by transform mode. */
		if (mode == FAST)
		{	
      /* Initialize */
      ctime = mysecond();
      nfsft_compute_wisdom(M,1000);
      ctime = mysecond() - ctime;
      #ifdef DEBUG
      printf("Time for initialization: %f sec.\n",ctime);
      #endif
      /* Precompute. */
      ctime = mysecond();
      nfsft_compute_wisdom(M,1000);
      ctime = mysecond() - ctime;
      #ifdef DEBUG
      printf("Time for precomputation: %f sec.\n",ctime);
      #endif
      /* Execute the plan. */
      ctime = mysecond();
      nfsft_adjoint(plan);
      ctime = mysecond() - ctime;
      #ifdef DEBUG
      printf("Time for fast algorithm: %f sec.\n",ctime);
      #endif
		}
		else if (mode == SLOW)
		{
      /* Execute the plan. */
      ctime = mysecond();
      ndsft_adjoint(plan);
      ctime = mysecond() - ctime;
      #ifdef DEBUG
      printf("Time for fast algorithm: %f sec.\n",ctime);
      #endif
		}
		else
		{
		  fprintf(stderr,"Wrong transform mode!\n");
		}

    for (k = 0; k <= M; k++)
    {  
      for (n = -k; n <= k; n++)
      {
        printf("%17.16f\n%17.16f\n",creal(F_HAT[n+M][k]),cimag(F_HAT[n+M][k]));
      }
    }
  }
  
  /* Destroy the plan. */
  nfsft_finalize(plan);

  
  /* Forget wisdom. */
  nfsft_forget_wisdom();
  
  fftw_free(result);
    
  /* Free data vectors. */
  for (n=-M;n<=M;n++)
  {
    fftw_free(F_HAT[n+M]);
  }
  
  fftw_free(F_HAT);
  
  /* Free nodes. */
  free(angles);
  
  return EXIT_SUCCESS;
}
