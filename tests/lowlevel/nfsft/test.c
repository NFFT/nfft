/*
 * Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* $Id$ */

/* Include ANSI-C headers. */
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include <complex.h>

/* Include NFFT3 header. */
#include "nfft3.h"
#include "util.h"

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
  /** */
  //int t;
  /** The bandwidth */
  int N;
  /** The number of nodes */
  int M;
  /** The original samples */
  double _Complex* f_orig;
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

  /*fprintf(stdout,"%d\n",testfiles);
  fflush(stdout);*/

  /* Test if successful. */
  if (testfiles == NULL)
  {
    fprintf(stderr,"Couldn't open %s to read test data filenames!\n",
      TESTFILES_NDSFT);
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
      /* Read in exponent. */
      //fscanf(file,"%d",&t);
      //fprintf(stdout,", t = %d",t);
      /* Read in bandwidth. */
      fscanf(file,"%d",&N);
      fprintf(stdout,", N = %4d",N);
      /* Read in number of nodes. */
      fscanf(file,"%d",&M);
      fprintf(stdout,", M = %5d ...",M);
      /* Precompute. */
      nfsft_precompute(N, THRESHOLD, 0U, 0U);
      /* Initialise plan. */
      nfsft_init_guru(&plan,N,M, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
        NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED /*| NFSFT_USE_NDFT*/,
        ((N>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT |
        FFT_OUT_OF_PLACE, 6);

      /* Read in spherical Fourier coefficients. */
      for (k = 0; k <= N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          fscanf(file,"%lf",&d1);
          fscanf(file,"%lf",&d2);
          plan.f_hat[NFSFT_INDEX(k,n,&plan)] = d1 +  I*d2;
        }
      }

      /*fprintf(stdout,"\n");
      for (k = 0; k < 2*N+1; k++)
      {
        for (n = 0; n < 2*N+1; n++)
        {
          fprintf(stdout,"%d ",(cabs(plan.f_hat[k*(2*N+1)+n])<1e-10)?0:1);
        }
        fprintf(stdout,"\n");
      }*/

      /* Read in nodes. */
      for (m = 0; m < plan.M_total; m++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        plan.x[2*m+1] = d1;
        plan.x[2*m] = d2;
      }

      /* Do precomputation for nodes. */
      nfsft_precompute_x(&plan);

      /* Read in reference samples. */
      f_orig = (double _Complex*)nfft_malloc(M*sizeof(double _Complex));
      for (m = 0; m < M; m++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        f_orig[m] = d1 + _Complex_I*d2;
        //fprintf(stdout,"f_orig[%d] = %lf + I*%lf\n",m,creal(f_orig[m]),cimag(f_orig[m]));
      }

      /* Close the file. */
      fclose(file);
      file = NULL;

      /* Execute the plan. */
      nfsft_trafo(&plan);

      /* Check result */
      fprintf(stdout," e_infty = %le,",nfft_error_l_infty_complex(f_orig,plan.f,M));
      fprintf(stdout," e_2 = %le",nfft_error_l_2_complex(f_orig,plan.f,M));

      //fprintf(stdout,"\n");
      /*for (m = 0; m < M; m++)
      {
        fprintf(stdout,"f[%d] = %lf + I*%lf, f_orig[%d] = %lf + I*%lf\n",
          m,creal(plan.f[m]),cimag(plan.f[m]),m,creal(f_orig[m]),cimag(f_orig[m]));
        if (cabs(plan.f[m]-f_orig[m]) > 0.0001)
        {
          fprintf(stdout," failed\n  f[%d] = %lf + I*%lf, f_orig[%d] = %lf + I*%lf\n",
            m,creal(plan.f[m]),cimag(plan.f[m]),m,creal(f_orig[m]),cimag(f_orig[m]));
          CU_FAIL("Wrong result");
        }
      }*/

      /* Destroy the plan. */
      nfsft_finalize(&plan);

      /* Forget precomputed data. */
      nfsft_forget();

      /* Free memory. */
      nfft_free(f_orig);
      f_orig = NULL;

      /* Test passed. */
      fprintf(stdout,"\n");
      fprintf(stderr,"ok");
    }
    else
    {
      fprintf(stdout," failed: Couldn't open file %s.\n",filename);
    }
  }
  fclose(testfiles);
  testfiles = NULL;
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
  double _Complex* f_hat_orig;
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

  /*fprintf(stdout,"%d\n",testfiles);
  fflush(stdout);*/

  /* Test if successful. */
  if (testfiles == NULL)
  {
    fprintf(stderr,"Couldn't open %s to read test data filenames!\n",
      TESTFILES_ADJOINT_NDSFT);
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
      fprintf(stdout,", N = %4d",N);

      /* Read in number of nodes. */
      fscanf(file,"%d",&M);
      fprintf(stdout,", M = %5d ...",M);

      /* Precompute. */
      nfsft_precompute(N,THRESHOLD,0U,0U);

      /* Initialise plan. */
      nfsft_init_advanced(&plan,N,M, NFSFT_MALLOC_X |
        NFSFT_MALLOC_F | NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED |
        NFSFT_ZERO_F_HAT | NFSFT_PRESERVE_F_HAT);

      /* Read in function samples. */
      for (m = 0; m < M; m++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        plan.f[m] = d1 + _Complex_I*d2;
      }

      /* Read in nodes. */
      for (m = 0; m < plan.M_total; m++)
      {
        fscanf(file,"%lf",&d1);
        fscanf(file,"%lf",&d2);
        plan.x[2*m+1] = d1;
        plan.x[2*m] = d2;
      }

      /* Do precomputation for nodes. */
      nfsft_precompute_x(&plan);

      /* Read in reference Fourier coefficients. */
      f_hat_orig = (double _Complex*) nfft_malloc(plan.N_total*sizeof(double _Complex));
      for (k = 0; k <= N; k++)
      {
        for (n = -k; n <= k; n++)
        {
          fscanf(file,"%lf",&d1);
          fscanf(file,"%lf",&d2);
          f_hat_orig[NFSFT_INDEX(k,n,&plan)] = d1 + _Complex_I*d2;
          //fprintf(stdout,"f_orig[%d] = %lf + I*%lf\n",m,creal(f_orig[m]),cimag(f_orig[m]));
        }
      }

      /* CLose the file. */
      fclose(file);

      /* Execute the plan. */
      nfsft_adjoint(&plan);

      /* Check result */
      fprintf(stdout," e_infty = %le,",nfft_error_l_infty_complex(f_hat_orig,
        plan.f_hat, plan.N_total));
      fprintf(stdout," e_2 = %le",nfft_error_l_2_complex(f_hat_orig,plan.f_hat,
        plan.N_total));
      //fprintf(stdout,"\n");
      /*for (n = -N; n <= N; n++)
      {
        for (k = -N; k <= N; k++)
        {
            fprintf(stdout,"f_hat[%3d,%3d] = %+lf + I*%+lf,\tf_hat_orig[%3d,%3d] = %+lf + I*%+lf,%5d,\teps=%le\n",
              k,n,
              creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]),
              cimag(plan.f_hat[NFSFT_INDEX(k,n,&plan)]),
              k,n,
              creal(f_hat_orig[NFSFT_INDEX(k,n,&plan)]),
              cimag(f_hat_orig[NFSFT_INDEX(k,n,&plan)]),
              NFSFT_INDEX(k,n,&plan),
              cabs(plan.f_hat[NFSFT_INDEX(k,n,&plan)]-f_hat_orig[NFSFT_INDEX(k,n,&plan)]));
        }
      }*/
      //nfft_vpr_complex(f_hat_orig, plan.N_total, "f_hat_orig");
      //nfft_vpr_complex(plan.f_hat, plan.N_total, "f_hat");

      /*fprintf(stdout,"\n");
      for (n = 0; n < plan.N_total; n++)
      {
        fprintf(stdout,"f_hat[%3d] = %+lf + I*%+lf,\tf_hat_orig[%3d] = %+lf + I*%+lf,\teps=%le\n",
          n,
          creal(plan.f_hat[n]),
          cimag(plan.f_hat[n]),
          n,
          creal(f_hat_orig[n]),
          cimag(f_hat_orig[n]),
          cabs(plan.f_hat[n]-f_hat_orig[n]));
      }*/
      /* Destroy the plan. */
      nfsft_finalize(&plan);
      /* Forget precomputed data. */
      nfsft_forget();
      /* Free memory. */
      nfft_free(f_hat_orig);
      /* Test passed. */
      fprintf(stdout,"\n");
      fprintf(stderr,"ok");
    }
    else
    {
      fprintf(stdout," failed: Couldn't open file %s.\n",filename);
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
  /* Add test for ndsft_trafo. */
  test_ndsft_trafo();
  /* Add test for ndsft_adjoint. */
  test_ndsft_adjoint();
  /* Exit the program. */
  return EXIT_SUCCESS;
}
