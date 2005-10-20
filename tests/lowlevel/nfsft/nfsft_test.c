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

/* Include NFFT3 header. */
#include "nfft3.h"

/* Include ANSI-C headers. */
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

/* Include CUnit header. */
#include <CUnit/CUnit.h>

/** Maximum filename length */
#define FILENAME_LENGTH_MAX 50

/** Name of the file containing the test data filenames for NDSFT. */
const char TESTFILES_NDSFT[] = "ndsft.txt\0";

/** Name of the file containing the test data filenames for adjoint NDSFT. */
const char TESTFILES_ADJOINT_NDSFT[] = "adjoint_ndsft.txt\0";

/* Default threshold. Not important here since NDSFT-algorithm is used. */
#define THRESHOLD 1000.0

void test_ndsft_trafo(void)
{
  /** The plan */
  nfsft_plan plan;
  /** The bandwidth */
  int N;
  /** The number of nodes */
  int M;
  /** The original samples */
  complex* f_orig;
  /** The degree \f$k\f$ */
  int k;
  /** The order \f$n\f$. */
  int n;
  /** The node index \f$m\f$*/
  int m;
  /** Auxilliary variables used to read in complex numbers. */
  double d1,d2;
  /** The file containing the names of the testdata files. */
  FILE *testfiles;
  /** The file containg the testcase data */
  FILE *file;
  /** Name of file containing test data. */
  char filename[FILENAME_LENGTH_MAX+1];
  //char filename[50] = "../../../../../../test.dat";
  //char filename[50] = "data/test_ndsft_0002_00010.dat";

  /* Tell what we're doing. */
  fprintf(stdout,"ndsft_trafo: Testing ndsft_trafo ...\n");

  /* Try to open file containing the names of the test data files. */
  testfiles = fopen(TESTFILES_NDSFT,"r");

  fprintf(stdout,"%d\n",testfiles);
  fflush(stdout);

  /* Test if successful. */
  if (testfiles == NULL)
  {
    CU_FAIL("Couldn't open %s to read test data filenames!\n");      
    return;
  }

  while (fscanf(testfiles,"%s",filename) == 1)
  { 
    fprintf(stdout,"filename = %s",filename);
    /* Open input file. */ 
    file = fopen(filename,"r");
    /* Check if file was opened successfully. */
    if (file != NULL)
    {
      /* Read in bandwidth. */
      fscanf(file,"%d",&N);
      fprintf(stdout,", N = %d",N);
      /* Read in number of nodes. */
      fscanf(file,"%d",&M);
      fprintf(stdout,", M = %d ...",M);
      /* Precompute. */
      nfsft_precompute(N,THRESHOLD,0U);
      /* Initialise plan. */
      nfsft_init_advanced(&plan,N,M,NFSFT_MALLOC_X | NFSFT_MALLOC_F | 
        NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED);
      /* Read in spherical Fourier coefficients. */
      for (k = 0; k <= N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          fscanf(file,"%lf",&d1);
          fscanf(file,"%lf",&d2);
          plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k] = d1 + I*d2;
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
      /* Read in reference samples. */
      f_orig = (complex*) malloc(M*sizeof(complex));
      for (m = 0; m < M; m++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        f_orig[m] = d1 + I*d2;
        //fprintf(stdout,"f_orig[%d] = %lf + I*%lf\n",m,creal(f_orig[m]),cimag(f_orig[m]));
      }
      
      /* CLose the file. */
      fclose(file);
      /* Execute the plan. */
      ndsft_trafo(&plan);
      /* Check result */
      for (m = 0; m < M; m++)
      {
        /*fprintf(stdout,"f[%d] = %lf + I*%lf, f_orig[%d] = %lf + I*%lf\n",
          m,creal(plan.f[m]),cimag(plan.f[m]),m,creal(f_orig[m]),cimag(f_orig[m]));*/
        if (cabs(plan.f[m]-f_orig[m]) > 0.0001)
        {
          fprintf(stdout," failed\n  f[%d] = %lf + I*%lf, f_orig[%d] = %lf + I*%lf\n",
            m,creal(plan.f[m]),cimag(plan.f[m]),m,creal(f_orig[m]),cimag(f_orig[m]));
          CU_FAIL("Wrong result");  
        }
      }    
      /* Destroy the plan. */
      nfsft_finalize(&plan);
      /* Forget precomputed data. */
      nfsft_forget();
      /* Free memory. */
      free(f_orig);
      /* Test passed. */
      fprintf(stdout," ok\n");
      CU_PASS("ok");
    }
    else
    {
      fprintf(stdout," failed: Couldn't open file %s.\n",filename);
      CU_FAIL("Couldn't open file!\n");
    }
  }  
}

void test_ndsft_adjoint(void)
{
  /** The plan */
  nfsft_plan plan;
  /** The bandwidth */
  int N;
  /** The number of nodes */
  int M;
  /** The original samples */
  complex* f_hat_orig;
  /** The degree \f$k\f$ */
  int k;
  /** The order \f$n\f$. */
  int n;
  /** The node index \f$m\f$*/
  int m;
  /** Auxilliary variables used to read in complex numbers. */
  double d1,d2;
  /** The file containing the names of the testdata files. */
  FILE *testfiles;
  /** The file containg the testcase data */
  FILE *file;
  /** Name of file containing test data. */
  char filename[FILENAME_LENGTH_MAX+1];
  //char filename[50] = "../../../../../../test.dat";
  //char filename[50] = "data/test_ndsft_0002_00010.dat";

  /* Tell what we're doing. */
  fprintf(stdout,"Testing ndsft_adjoint ...\n");

  /* Try to open file containing the names of the test data files. */
  testfiles = fopen(TESTFILES_ADJOINT_NDSFT,"r");

  fprintf(stdout,"%d\n",testfiles);
  fflush(stdout);

  /* Test if successful. */
  if (testfiles == NULL)
  {
    CU_FAIL("Couldn't open %s to read test data filenames!\n");      
    return;
  }

  while (fscanf(testfiles,"%s",filename) == 1)
  { 
    fprintf(stdout,"filename = %s",filename);
    /* Open input file. */ 
    file = fopen(filename,"r");
    /* Check if file was opened successfully. */
    if (file != NULL)
    {
      /* Read in bandwidth. */
      fscanf(file,"%d",&N);
      fprintf(stdout,", N = %d",N);
      /* Read in number of nodes. */
      fscanf(file,"%d",&M);
      fprintf(stdout,", M = %d ...",M);
      /* Precompute. */
      nfsft_precompute(N,THRESHOLD,0U);
      /* Initialise plan. */
      nfsft_init_advanced(&plan,N,M,NFSFT_MALLOC_X | NFSFT_MALLOC_F | 
        NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED);
      /* Read in function samples. */
      for (m = 0; m < M; m++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        plan.f[m] = d1 + I*d2;
      }
      /* Read in nodes. */
      for (m = 0; m < plan.M_total; m++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        plan.x[2*m] = d1;
        plan.x[2*m+1] = d2;
      }
      /* Read in reference Fourier coefficients. */
      f_hat_orig = (complex*) malloc((N+1)*(N+1)*sizeof(complex));
      for (k = 0; k <= N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          fscanf(file,"%lf",&d1);
          fscanf(file,"%lf",&d2);
          f_hat_orig[k*k+n+k] = d1 + I*d2;
          //fprintf(stdout,"f_orig[%d] = %lf + I*%lf\n",m,creal(f_orig[m]),cimag(f_orig[m]));
        }
      }
      
      /* CLose the file. */
      fclose(file);
      /* Execute the plan. */
      ndsft_adjoint(&plan);
      /* Check result */
      //fprintf(stdout,"\n");
      for (k = 0; k <= N; k++)
      {
        for (n = -k; n <= k; n++)
        {
            /*fprintf(stdout,"f_hat[%d] = %lf + I*%lf, f_hat_orig[%d] = %lf + I*%lf\n",
              (2*plan.NPT+1)*(n+plan.N)+plan.NPT+k,
              creal(plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k]),
              cimag(plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k]),
              k*k+n+k,creal(f_hat_orig[k*k+n+k]),cimag(f_hat_orig[k*k+n+k]));*/
          if (cabs(plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k]-
            f_hat_orig[k*k+n+k]) > 0.0001)
          {
            fprintf(stdout," failed\n  f_hat[%d] = %lf + I*%lf, f_hat_orig[%d] = %lf + I*%lf\n",
              (2*plan.NPT+1)*(n+plan.N)+plan.NPT+k,
              creal(plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k]),
              cimag(plan.f_hat[(2*plan.NPT+1)*(n+plan.N)+plan.NPT+k]),
              k*k+n+k,creal(f_hat_orig[k*k+n+k]),cimag(f_hat_orig[k*k+n+k]));
            CU_FAIL("Wrong result");  
          }
        }
      }    
      /* Destroy the plan. */
      nfsft_finalize(&plan);
      /* Forget precomputed data. */
      nfsft_forget();
      /* Free memory. */
      free(f_hat_orig);
      /* Test passed. */
      fprintf(stdout," ok\n");
      CU_PASS("ok");
    }
    else
    {
      fprintf(stdout," failed: Couldn't open file %s.\n",filename);
      CU_FAIL("Couldn't open file!\n");
    }
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
  /* Initialise registry. */
  CU_initialize_registry();
  /* Create test suite. */
  CU_pSuite suite = CU_add_suite("NDSFT", NULL, NULL);
  /* Add test for ndsft_trafo. */
  CU_add_test(suite,"NDSFT",test_ndsft_trafo);
  /* Add test for ndsft_adjoint. */
  CU_add_test(suite,"adjoint NDSFT",test_ndsft_adjoint);
  /* Run the tests. */
  CU_automated_run_tests();
  /* Cleanup registry. */
  CU_cleanup_registry();   
  return EXIT_SUCCESS;
}
