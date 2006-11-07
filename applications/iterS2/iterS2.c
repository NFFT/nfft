/* $Id$
 *
 * iterS2 - Iterative reconstruction on the sphere S2
 *
 * Copyright (C) 2006 Jens Keiner
 */

/** 
 * \defgroup applications_iterS2_matlab iterS2_matlab
 * \ingroup applications_iterS2
 * \{
 */

/* Include standard C headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Include NFFT3 library header. */
#include "nfft3.h"

/* Include NFFT 3 utilities headers. */
#include "util.h"

/** Enumeration for parameter values */
enum boolean {NO = 0, YES = 1};

/**
 * The main program.
 *
 * \param argc The number of arguments
 * \param argv An array containing the arguments as C-strings
 *
 * \return Exit code
 *
 * \author Jens Keiner
 */
int main (int argc, char **argv)
{
  int T;
  int N;
  int M;
  
  int t;                       /* Index variable for testcases                */
  double time;                 /* Time for fast algorithm in seconds          */
  double err_f;                /* Error E_infty for fast algorithm            */

  double *w;                   /* The weights from the normal equation        */
  complex double *fo;          /* The original function values                */
  complex double *f;           /* The function values                         */
  complex double *f_hat;       /* The spherical Fourier coefficients          */
  double *x;                   /* The nodes                                   */
  nfsft_plan plan;             /* NFSFT plan                                  */
  int j;                       /*                                             */
  int k;                       /*                                             */
  int m;                       /*                                             */
  int n;                       /*                                             */
  int use_nfsft;               /*                                             */
  int use_nfft;                /*                                             */
  int use_fpt;                 /*                                             */
  int cutoff;                  /**< The current NFFT cut-off parameter        */
  double threshold;            /**< The current NFSFT threshold parameter     */
  double re;
  double im;

  /* Read the number of testcases. */
  fscanf(stdin,"testcases=%d\n",&T);
  fprintf(stdout,"%d\n",T);

  /* Process each testcase. */
  for (t = 0; t < T; t++)
  {
    /* Read the bandwidth. */
    fscanf(stdin,"bandwidth=%d\n",&N);
    fprintf(stdout,"%d\n",N);
  
    /* Read the number of nodes. */
    fscanf(stdin,"nodes=%d\n",&M);
    fprintf(stdout,"%d\n",M);

    /* Check if the fast transform shall be used. */
    fscanf(stdin,"nfsft=%d\n",&use_nfsft);
    fprintf(stdout,"%d\n",use_nfsft);
    if (use_nfsft != NO)
    {
      /* Check if the NFFT shall be used. */
      fscanf(stdin,"nfft=%d\n",&use_nfft);
      fprintf(stdout,"%d\n",use_nfsft);
      if (use_nfft != NO)
      {
        /* Read the cut-off parameter. */
        fscanf(stdin,"cutoff=%d\n",&cutoff);
        fprintf(stdout,"%d\n",cutoff);
      }
      else
      {
        /* TODO remove this */
        /* Initialize unused variable with dummy value. */
        cutoff = 1;
      }
      /* Check if the fast polynomial transform shall be used. */
      fscanf(stdin,"fpt=%d\n",&use_fpt);
      fprintf(stdout,"%d\n",use_fpt);
      if (use_fpt != NO)
      {
        /* Read the NFSFT threshold parameter. */
        fscanf(stdin,"threshold=%lf\n",&threshold);
        fprintf(stdout,"%lf\n",threshold);
      }
      else
      {
        /* TODO remove this */
        /* Initialize unused variable with dummy value. */
        threshold = 1000.0;
      }
    }
    else
    {
      /* TODO remove this */
      /* Set dummy values. */
      use_nfft = NO;
      use_fpt = NO;
      cutoff = 3;
      threshold = 1000.0;
    }

    /* Read the bandwidth. */
    fscanf(stdin,"bandwidth=%d\n",&N);
    fprintf(stdout,"%d\n",N);

    /* Read the number of nodes. */
    fscanf(stdin,"nodes=%d\n",&M);
    fprintf(stdout,"%d\n",M);

    /* Allocate memory for data structures. */
    w = (double*) malloc(M*sizeof(double));
    fo = (complex double*) malloc(M*sizeof(complex double));
    f = (complex double*) malloc(M*sizeof(complex double));
    x = (double*) malloc(2*M*sizeof(double));
    f_hat = (complex double*) malloc(NFSFT_F_HAT_SIZE(N)*sizeof(complex double));

    /* Read the nodes and function values. */
    for (j = 0; j < M; j++)
    {
      fscanf(stdin,"%le %le %le %le\n",&x[2*j+1],&x[2*j],&re,&im);
      f[j] = re + I * im;
      fprintf(stdout,"%le %le %le %le\n",x[2*j+1],x[2*j],creal(f[j]),cimag(f[j]));
    }

    /* Do precomputation. */
    /*nfsft_precompute(m_max,threshold,
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U+++NFSFT_NO_DIRECT_ALGORITHM)), 0U);*/


    /* Init transform plans. */
    /*nfsft_init_guru(&plan_adjoint, m[im],ld[ild][0],
      ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
      ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)),
      PRE_PHI_HUT | PRE_PSI | FFTW_INIT |
      FFT_OUT_OF_PLACE, cutoff);
    nfsft_init_guru(&plan,m[im],ld[ild][1],
      ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
      ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)),
      PRE_PHI_HUT | PRE_PSI | FFTW_INIT |
      FFT_OUT_OF_PLACE,
       cutoff);
    plan_adjoint.f_hat = f_hat;
    plan_adjoint.x = eta;
    plan_adjoint.f = b;
    plan.f_hat = f_hat;
    plan.x = xi;
    plan.f = f_m;
    nfsft_precompute_x(&plan_adjoint);
    nfsft_precompute_x(&plan);*/

    /* Delete precomputed data. */
    /*nfsft_forget();*/

    /* Free data arrays. */
    free(w);
    free(fo);
    free(f);
    free(x);
    free(f_hat);

  } /* Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
/* \} */
