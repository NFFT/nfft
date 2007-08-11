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

#include "legendre.h"

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
  int M2;

  int t;                       /* Index variable for testcases                */
  double time;                 /* Time for fast algorithm in seconds          */
  double err_f;                /* Error E_infty for fast algorithm            */
  nfsft_plan plan;             /* NFSFT plan                                  */
  nfsft_plan plan2;            /* NFSFT plan                                  */
  infsft_plan iplan;           /* NFSFT plan                                  */
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
  double a,b;
  double *scratch;
  double xs;
  double *ys;
  double *temp;
  double _Complex *temp2;
  int qlength;
  double *qweights;
  fftw_plan fplan;
  fpt_set set;
  int npt;
  int npt_exp;
  double *alpha, *beta, *gamma;

  /* Read the number of testcases. */
  fscanf(stdin,"testcases=%d\n",&T);
  fprintf(stderr,"%d\n",T);

  /* Process each testcase. */
  for (t = 0; t < T; t++)
  {
    /* Check if the fast transform shall be used. */
    fscanf(stdin,"nfsft=%d\n",&use_nfsft);
    fprintf(stderr,"%d\n",use_nfsft);
    if (use_nfsft != NO)
    {
      /* Check if the NFFT shall be used. */
      fscanf(stdin,"nfft=%d\n",&use_nfft);
      fprintf(stderr,"%d\n",use_nfsft);
      if (use_nfft != NO)
      {
        /* Read the cut-off parameter. */
        fscanf(stdin,"cutoff=%d\n",&cutoff);
        fprintf(stderr,"%d\n",cutoff);
      }
      else
      {
        /* TODO remove this */
        /* Initialize unused variable with dummy value. */
        cutoff = 1;
      }
      /* Check if the fast polynomial transform shall be used. */
      fscanf(stdin,"fpt=%d\n",&use_fpt);
      fprintf(stderr,"%d\n",use_fpt);
      if (use_fpt != NO)
      {
        /* Read the NFSFT threshold parameter. */
        fscanf(stdin,"threshold=%lf\n",&threshold);
        fprintf(stderr,"%lf\n",threshold);
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
    fprintf(stderr,"%d\n",N);

    /* Do precomputation. */
    nfsft_precompute(N,threshold,
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U/*NFSFT_NO_DIRECT_ALGORITHM*/)), 0U);

    /* Read the number of nodes. */
    fscanf(stdin,"nodes=%d\n",&M);
    fprintf(stderr,"%d\n",M);

    /* */
    if ((N+1)*(N+1) > M)
    {
      nfft_next_power_of_2_exp(N, &npt, &npt_exp);
      fprintf(stderr, "npt = %d, npt_exp = %d\n", npt, npt_exp);
      fprintf(stderr,"Optimal interpolation!\n");
      scratch = (double*) nfft_malloc(4*sizeof(double));
      ys = (double*) nfft_malloc((N+1)*sizeof(double));
      temp = (double*) nfft_malloc((2*N+1)*sizeof(double));
      temp2 = (double _Complex*) nfft_malloc((N+1)*sizeof(double _Complex));

      a = 0.0;
      for (j = 0; j <= N; j++)
      {
        xs = 2.0 + (2.0*j)/(N+1);
        ys[j] = (2.0-((j == 0)?(1.0):(0.0)))*4.0*nfft_bspline(4,xs,scratch);
        //fprintf(stdout,"%3d: g(%le) = %le\n",j,xs,ys[j]);
        a += ys[j];
      }
      //fprintf(stdout,"a = %le\n",a);
      for (j = 0; j <= N; j++)
      {
        ys[j] *= 1.0/a;
      }

      qlength = 2*N+1;
      qweights = (double*) nfft_malloc(qlength*sizeof(double));

      fplan = fftw_plan_r2r_1d(N+1, qweights, qweights, FFTW_REDFT00, 0U);
      for (j = 0; j < N+1; j++)
      {
        qweights[j] = -2.0/(4*j*j-1);
      }
      fftw_execute(fplan);
      qweights[0] *= 0.5;

      for (j = 0; j < N+1; j++)
      {
        qweights[j] *= 1.0/(2.0*N+1.0);
        qweights[2*N+1-1-j] = qweights[j];
      }

      fplan = fftw_plan_r2r_1d(2*N+1, temp, temp, FFTW_REDFT00, 0U);
      for (j = 0; j <= N; j++)
      {
        temp[j] = ((j==0 || j == 2*N)?(1.0):(0.5))*ys[j];
      }
      for (j = N+1; j < 2*N+1; j++)
      {
        temp[j] = 0.0;
      }
      fftw_execute(fplan);

      for (j = 0; j < 2*N+1; j++)
      {
        temp[j] *= qweights[j];
      }

      fftw_execute(fplan);

      for (j = 0; j < 2*N+1; j++)
      {
        temp[j] *= ((j==0 || j == 2*N)?(1.0):(0.5));
        if (j <= N)
        {
          temp2[j] = temp[j];
        }
      }

      set = fpt_init(0, npt_exp, 0U);

      alpha = (double*) nfft_malloc((N+2)*sizeof(double));
      beta = (double*) nfft_malloc((N+2)*sizeof(double));
      gamma = (double*) nfft_malloc((N+2)*sizeof(double));

      alpha_al_row(alpha, N, 0);
      beta_al_row(beta, N, 0);
      gamma_al_row(gamma, N, 0);

      fpt_precompute(set, 0, alpha, beta, gamma, 0, 1000.0);

      fpt_transposed(set,0, temp2, temp2, N, 0U);

      for (j = 0; j <= N; j++)
      {
        //temp2[j] *= (2*j+1.0)/2.0;
        //fprintf(stdout,"temp2[%d] = %le\n",j,cabs(temp2[j]));
      }

      fpt_finalize(set);

      nfft_free(alpha);
      nfft_free(beta);
      nfft_free(gamma);

      fftw_destroy_plan(fplan);

      nfft_free(scratch);
      nfft_free(qweights);
      nfft_free(ys);
      nfft_free(temp);
      //nfft_free(temp2);
      //return EXIT_SUCCESS;
    }

    /* Init transform plans. */
    nfsft_init_guru(&plan, N, M,
      ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
      ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)) | NFSFT_MALLOC_F | NFSFT_MALLOC_X |
      NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_ZERO_F_HAT,
      PRE_PHI_HUT | PRE_PSI | FFTW_INIT |
      FFT_OUT_OF_PLACE,
      cutoff);

    if ((N+1)*(N+1) > M)
    {
      infsft_init_advanced (&iplan, &plan, CGNE | PRECOMPUTE_DAMP);
    }
    else
    {
      infsft_init_advanced (&iplan, &plan, CGNR | PRECOMPUTE_WEIGHT | PRECOMPUTE_DAMP);
    }

    /* Read the nodes and function values. */
    for (j = 0; j < M; j++)
    {
      fscanf(stdin,"%le %le %le %le\n",&plan.x[2*j+1],&plan.x[2*j],&re,&im);
      plan.x[2*j+1] = plan.x[2*j+1]/(2.0*PI);
      plan.x[2*j] = plan.x[2*j]/(2.0*PI);
      if (plan.x[2*j] >= 0.5)
      {
        plan.x[2*j] = plan.x[2*j] - 1;
      }
      iplan.y[j] = re + _Complex_I * im;
      fprintf(stderr,"%le %le %le %le\n",plan.x[2*j+1],plan.x[2*j],
        creal(iplan.y[j]),cimag(iplan.y[j]));
    }

    /* Read the number of nodes. */
    fscanf(stdin,"nodes_eval=%d\n",&M2);
    fprintf(stderr,"%d\n",M2);

    /* Init transform plans. */
    nfsft_init_guru(&plan2, N, M2,
      ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
      ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)) | NFSFT_MALLOC_F | NFSFT_MALLOC_X |
      NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_ZERO_F_HAT,
      PRE_PHI_HUT | PRE_PSI | FFTW_INIT |
      FFT_OUT_OF_PLACE,
      cutoff);

    /* Read the nodes and function values. */
    for (j = 0; j < M2; j++)
    {
      fscanf(stdin,"%le %le\n",&plan2.x[2*j+1],&plan2.x[2*j]);
      plan2.x[2*j+1] = plan2.x[2*j+1]/(2.0*PI);
      plan2.x[2*j] = plan2.x[2*j]/(2.0*PI);
      if (plan2.x[2*j] >= 0.5)
      {
        plan2.x[2*j] = plan2.x[2*j] - 1;
      }
      fprintf(stderr,"%le %le\n",plan2.x[2*j+1],plan2.x[2*j]);
    }

    nfsft_precompute_x(&plan);

    nfsft_precompute_x(&plan2);

    /* Frequency weights. */
    if ((N+1)*(N+1) > M)
    {
      /* Compute Voronoi weights. */
      //nfft_voronoi_weights_S2(iplan.w, plan.x, M);

      /* Print out Voronoi weights. */
      /*a = 0.0;
      for (j = 0; j < plan.M_total; j++)
      {
        fprintf(stderr,"%le\n",iplan.w[j]);
        a += iplan.w[j];
      }
      fprintf(stderr,"sum = %le\n",a);*/

      for (j = 0; j < plan.N_total; j++)
      {
        iplan.w_hat[j] = 0.0;
      }

      for (k = 0; k <= N; k++)
      {
        for (j = -k; j <= k; j++)
        {
          iplan.w_hat[NFSFT_INDEX(k,j,&plan)] = 1.0/(pow(k+1.0,2.0)); /*temp2[j]*/;
        }
      }
    }
    else
    {
      for (j = 0; j < plan.N_total; j++)
      {
        iplan.w_hat[j] = 0.0;
      }

      for (k = 0; k <= N; k++)
      {
        for (j = -k; j <= k; j++)
        {
          iplan.w_hat[NFSFT_INDEX(k,j,&plan)] = 1/(pow(k+1.0,2.5));
        }
      }

      /* Compute Voronoi weights. */
      nfft_voronoi_weights_S2(iplan.w, plan.x, M);

      /* Print out Voronoi weights. */
      a = 0.0;
      for (j = 0; j < plan.M_total; j++)
      {
        fprintf(stderr,"%le\n",iplan.w[j]);
        a += iplan.w[j];
      }
      fprintf(stderr,"sum = %le\n",a);
    }

    fprintf(stderr, "N_total = %d\n", plan.N_total);
    fprintf(stderr, "M_total = %d\n", plan.M_total);

    /* init some guess */
    for (k = 0; k < plan.N_total; k++)
    {
      iplan.f_hat_iter[k] = 0.0;
    }

    /* inverse trafo */
    infsft_before_loop(&iplan);

    /*for (k = 0; k < plan.M_total; k++)
    {
      printf("%le %le\n",creal(iplan.r_iter[k]),cimag(iplan.r_iter[k]));
    }*/

    for (m = 0; m < 29; m++)
    {
      fprintf(stderr,"Residual ||r||=%e,\n",sqrt(iplan.dot_r_iter));
      infsft_loop_one_step(&iplan);
    }

    /*NFFT_SWAP_complex(iplan.f_hat_iter, plan.f_hat);
    nfsft_trafo(&plan);
    NFFT_SWAP_complex(iplan.f_hat_iter, plan.f_hat);

    a = 0.0;
    b = 0.0;
    for (k = 0; k < plan.M_total; k++)
    {
      printf("%le %le %le\n",cabs(iplan.y[k]),cabs(plan.f[k]),
        cabs(iplan.y[k]-plan.f[k]));
      a += cabs(iplan.y[k]-plan.f[k])*cabs(iplan.y[k]-plan.f[k]);
      b += cabs(iplan.y[k])*cabs(iplan.y[k]);
    }

    fprintf(stderr,"relative error in 2-norm: %le\n",a/b);*/

    NFFT_SWAP_complex(iplan.f_hat_iter, plan2.f_hat);
    nfsft_trafo(&plan2);
    NFFT_SWAP_complex(iplan.f_hat_iter, plan2.f_hat);
    for (k = 0; k < plan2.M_total; k++)
    {
      fprintf(stdout,"%le\n",cabs(plan2.f[k]));
    }

    infsft_finalize (&iplan);

    nfsft_finalize(&plan);

    nfsft_finalize(&plan2);

    /* Delete precomputed data. */
    nfsft_forget();

    if ((N+1)*(N+1) > M)
    {
      nfft_free(temp2);
    }

  } /* Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
/* \} */
