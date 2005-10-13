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

#include "nfft3.h"

#include <stdlib.h>
#include <stdio.h>

#define N 256
#define M 100
#define THRESHOLD 1000.0

/**
 * Test function for direct NDSFT
 */
void test_ndsft_trafo(void)
{
  /** The plan */
  nfsft_plan plan;
  /** The file containg the testcase data */
  FILE file;
  /** The bandwidth */
  int N;
  /** The number of nodes */
  int M;
  /** The original samples */
  complex* f_orig;
  int k,n,m;
  double re,im;

  /* Read input data from file. */ 
  file = fopen("test001.dat","r");
  if (file != NULL)
  {
    /* Read in bandwidth. */
    fscanf(file,"%d",&N);
    printf(stdout,"N = %d",N);
    /* Read in number of nodes. */
    fscanf(file,"%d",&M);
    printf(stdout,", M = %d",M);
    /* Precompute. */
    nfsft_precompute(N,THRESHOLD,0U);
    /* Initialise plan. */
    nfsft_plan = nfsft_init(&plan,N,M);
    /* Fill in spherical Fourier coefficients */
    fprintf("/n");
    for (k = 0; k <= N; k++)
    {
      for (n = -k; n <= k; n++)
      {
        fscanf(file,"%fl",$re);
        fscanf(file,"%fl",$im);
        plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k] = re + I*im;
        fprintf(stdout,"f[%d] = %f + I*%f\n",
          (2*plan.NPT+1)*(n+plan.N)+plan.NPT+k,re,im);
      }
    }
    fclose(file);
  }
} 
 
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
  test_ndsft_trafo();
  return EXIT_SUCCESS;
}
