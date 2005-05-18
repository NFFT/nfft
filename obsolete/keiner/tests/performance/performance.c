/* 
   performance - Performance test for NFSFT

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

#define M_MIN 512
#define M_STRIDE 1
#define M_MAX 512

#define D_MIN 256
#define D_STRIDE 256
#define D_MAX 1024

#define THRESHOLD 1000

#define MACRO_TEST(FILENAME,TRAFO_NAME) \
printf("Testing TRAFO_NAME:\n"); \
fflush(stdout); \
sprintf(filename,FILENAME); \
file = fopen(filename,"w"); \
if (file != NULL)  \
{ \
  fprintf(file,"%d\n",((D_MAX-D_MIN)/D_STRIDE+1)*((M_MAX-M_MIN)/M_STRIDE+1)); \
  for (D = D_MIN; D <= D_MAX; D = D + D_STRIDE) \
  {  \
    for (M = M_MIN; M <= M_MAX; M = M + M_STRIDE) \
    { \
      printf("D = %d, M = %d\n",D,M); \
      plan = nfsft_init(D, M, angles, f_hat, f, 0U); \
      ctime = mysecond(); \
      TRAFO_NAME(plan); \
      ctime = mysecond() - ctime; \
      fprintf(file,"%10d %4d %10.4f\n",D,M,ctime); \
      nfsft_finalize(plan); \
    } \
  } \
  fclose(file); \
} \
printf("done\n");


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
  /** Next greater power of two with respect to M_MAX */
  const int N_MAX = 1<<((int)ceil(log((double)M_MAX)/log(2.0)));
  /* Bandwidth */
  int M;
  /** Next greater power of two with respect to M */
  //int N;

  /* Number of nodes */
  int D;
	
  complex **f_hat;
 	complex *f;
  /** Array of angles defining the nodes. */
  double *angles;
 	nfsft_plan plan;
  int n,k,d;
  /** Used to measure computation time. */
  double ctime;
  /** 
   * Used to store the filename of a file containing Gauss-Legendre nodes and 
   * weights.
   */ 
  char filename[20];
  /** File handle for reading quadrature nodes and weights. */
  FILE *file;
  
	 
  /* Initialize random number generator. */
  srand48(time(NULL));
  
  /* Allocate data structures. */
  f_hat = (complex**) malloc((2*M_MAX+1)*sizeof(complex*));
  for (n = -M_MAX; n <= M_MAX; n++)
  {
    f_hat[n+M_MAX] = (complex*) malloc((N_MAX+1)*sizeof(complex));
  }  
  f = (complex*) malloc(D_MAX*sizeof(complex));
  angles = (double*) malloc(2*D_MAX*sizeof(double));
  
  /* Compute random Fourier coefficients. */
  for (n = -M_MAX; n <= M_MAX; n++)
  {
    for (k = abs(n); k <= N_MAX; k++)
    {
      f_hat[n+M_MAX][k] = drand48() + I*drand48();
    }
  }  
  
  nfsft_precompute(M_MAX,THRESHOLD);  
  
  MACRO_TEST("ndsft.dat",         ndsft_trafo)
  MACRO_TEST("ndsft_adjoint.dat", ndsft_adjoint)
  MACRO_TEST("nfsft.dat",         nfsft_trafo)
  MACRO_TEST("nfsft_adjoint.dat", nfsft_adjoint)
		
  free(angles);
  free(f);
  for (n = -M_MAX; n <= M_MAX; n++)
  {
    free(f_hat[n+M_MAX]);
  }   
	
  return EXIT_SUCCESS;
}
