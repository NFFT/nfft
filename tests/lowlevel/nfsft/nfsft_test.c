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

#define THRESHOLD 1000.0

/**
 * Test function for direct NDSFT
 */
void test_ndsft_trafo(void)
{
  /** The plan */
  nfsft_plan plan;
  /** The file containg the testcase data */
  FILE *file;
  /** The bandwidth */
  int N;
  /** The number of nodes */
  int M;
  /** The original samples */
  complex* f_orig;
  int k,n,m;
  double d1,d2;

  /* Read input data from file. */ 
  file = fopen("test001.dat","r");
  fprintf(stdout,"file = %p\n",file);
  if (file != NULL)
  {
    /* Read in bandwidth. */
    fscanf(file,"%d\n",&N);
    fprintf(stdout,"N = %d",N);
    /* Read in number of nodes. */
    fscanf(file,"%d\n",&M);
    fprintf(stdout,", M = %d",M);
    /* Precompute. */
    nfsft_precompute(N,THRESHOLD,0U);
    /* Initialise plan. */
    nfsft_init(&plan,N,M);
    /* Read in spherical Fourier coefficients. */
    fprintf(stdout,"\n");
    for (k = 0; k <= N; k++)
    {
      for (n = -k; n <= k; n++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        /*fprintf(stdout,"Scanned: %fl + I*%fl, Index: %d\n",re,im,
          2*plan.NPT+1);*/
        plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k] = d1 + I*d2;
        fprintf(stdout,"f[%d] = %f + I*%f\n",
          (2*plan.NPT+1)*(n+plan.N)+plan.NPT+k,
          creal(plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k]),
          cimag(plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k]));
      }
    }
    /* Read in nodes. */
    for (m = 0; m < plan.M_total; m++)
    {
      fscanf(file,"%lf",&d1);
      fscanf(file,"%lf",&d2);
      plan.x[2*m] = d1;
      plan.x[2*m+1] = d2;
    }
    /* Execute the plan. */
    ndsft_trafo(&plan);
    /* Display result. */
    for (m = 0; m < plan.M_total; m++)
    {
      fprintf(stdout,"f[%d] = %lf + I*%lf\n",m,creal(plan.f[m]),cimag(plan.f[m]));
    }    
    /* Destroy the plan. */
    nfsft_finalize(&plan);
    /* CLose the file. */
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
