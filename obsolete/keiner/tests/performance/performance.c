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

#ifdef HAVE_COMPLEX_H
#  include <complex.h>
#else
#  error Need complex.h
#endif

#ifdef HAVE_COMPLEX_H
#  include <time.h>
#else
#  error Need time.h
#endif

/* Package headers */
#include "nfsft.h"
#include "util.h"

/* Default values */
#define M_MIN 1
#define M_STRIDE 1
#define M_MAX 32

#define D_MIN 256
#define D_STRIDE 256
#define D_MAX 1024

#define THRESHOLD 1000

/* Macro for the tests */
#define MACRO_TEST(FILENAME,TRAFO_FUNCTION,TRAFO_NAME) \
printf("Testing %s:\n",TRAFO_NAME); \
fflush(stdout); \
sprintf(filename,FILENAME); \
file = fopen(filename,"w"); \
if (file != NULL)  \
{ \
  fprintf(file,"%d\n",((d_max-d_min)/d_stride+1)*((m_max-m_min)/m_stride+1)); \
  fclose(file); \
  for (D = d_min; D <= d_max; D = D + d_stride) \
  {  \
    for (M = m_min; M <= m_max; M = M + m_stride) \
    { \
      plan = nfsft_init(D, M, angles, f_hat, f, 0U); \
      ctime = mysecond(); \
      TRAFO_FUNCTION(plan); \
      ctime = mysecond() - ctime; \
      printf("D = %10d, M = %4d\n",D,M); \
      fopen(filename,"a"); \
      fprintf(file,"%10d %4d %10.4f\n",D,M,ctime); \
      fclose(file); \
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
  int m_min; 
  int m_max; 
  int m_stride; 
  int d_min; 
  int d_max; 
  int d_stride; 
  int threshold;
  /** Next greater power of two with respect to m_max */
  int N_MAX;
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
  

  if (argc == 1)
  {
    m_min = M_MIN;
    m_max = M_MAX; 
    m_stride = M_STRIDE; 
    d_min = D_MIN; 
    d_max = D_MAX; 
    d_stride = D_STRIDE; 
    threshold = THRESHOLD;    
  }
  else if (argc == 8)
  {
    sscanf(argv[1],"%d",&m_min);
    sscanf(argv[2],"%d",&m_max);
    sscanf(argv[3],"%d",&m_stride);
    sscanf(argv[4],"%d",&d_min);
    sscanf(argv[5],"%d",&d_max);
    sscanf(argv[6],"%d",&d_stride);
    sscanf(argv[7],"%d",&threshold);
  }
  else
  {
    fprintf(stderr,"Performance - Performance test for NFSFT\n");
    fprintf(stderr,"Usage: performance M_MIN M_MAX M_STRIDE D_MIN D_MAX D_STRIDE THRESHOLD\n");
    return -1;
  }  
      
  /** Next greater power of two with respect to m_max */
  N_MAX = 1<<((int)ceil(log((double)m_max)/log(2.0)));
  
  /* Initialize random number generator. */
  srand48(time(NULL));
  
  /* Allocate data structures. */
  f_hat = (complex**) malloc((2*m_max+1)*sizeof(complex*));
  for (n = -m_max; n <= m_max; n++)
  {
    f_hat[n+m_max] = (complex*) malloc((N_MAX+1)*sizeof(complex));
  }  
  f = (complex*) malloc(d_max*sizeof(complex));
  angles = (double*) malloc(2*d_max*sizeof(double));
  
  /* Compute random Fourier coefficients. */
  for (n = -m_max; n <= m_max; n++)
  {
    for (k = abs(n); k <= N_MAX; k++)
    {
      f_hat[n+m_max][k] = drand48() + I*drand48();
    }
  }  
  
  nfsft_precompute(m_max,threshold);  
  
  MACRO_TEST("ndsft.dat",         ndsft_trafo,   "NDSFT")
  MACRO_TEST("ndsft_adjoint.dat", ndsft_adjoint, "adjoint NDSFT")
  MACRO_TEST("nfsft.dat",         nfsft_trafo,   "NFSFT")
  MACRO_TEST("nfsft_adjoint.dat", nfsft_adjoint, "adjoint NFSFT")
		
  free(angles);
  free(f);
  for (n = -m_max; n <= m_max; n++)
  {
    free(f_hat[n+m_max]);
  }   
	
  return EXIT_SUCCESS;
}
