/* $Id: fastsumS2.c 741 2006-04-05 07:00:11Z keiner $
 *
 * fastsumS2 - Fast summation of spherical radial functions
 *
 * Copyright (C) 2005 Jens Keiner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/* Include standard C headers. */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/* Include NFFT3 library header. */
#include "nfft3.h"

/* Include NFFT 3 utilities headers. */
#include "util.h"

/** The Fourier-Legendre coefficients of the Abel-Poisson kernel */
#define SYMBOL_ABEL_POISSON(k,h) (pow(h,k))

/** The Fourier-Legendre coefficients of the singularity kernel */
#define SYMBOL_SINGULARITY(k,h) ((2.0/(2*k+1))*pow(h,k))

/* Flags for the different quadrature grid types */

/** Enumeration for parameter values */
enum boolean {NO = 0, YES = 1};

/** Enumeration for quadrature grid types */
enum gridtype {GRID_GAUSS_LEGENDRE = 0, GRID_CLENSHAW_CURTIS = 1,
  GRID_HEALPIX = 2, GRID_EQUIDISTRIBUTION = 3};

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
  int tc_max;                  /**< The number of testcases                   */
  int tc;                      /**< Index variable for testcases              */
  int *m;                      /**< The array containing the cut-off degrees  *
                                    \f$M\f$.                                  */
  int m_max;                   /**< The maximum cut-off degree \f$M\f$ for the*
                                    current testcase                          */
  int im_max;                  /**< The maximum index for \code m             */
  int im;                      /**< Index variable for \code m                */
  int use_nfsft;               /**< */
  int use_nfft;                /**< */
  int use_fpt;                 /**< */
  int gridtype;                /**< */
  int cutoff;                  /**< The current NFFT cut-off parameter        */
  double threshold;            /**< The current NFSFT threshold parameter     */
  int repetitions;             /**< */
  int n_max;                   /**< Next greater power of two with respect to *
                                    m_max                                     */
  double t_avg;                /**< The average computation time needed       */
  double err;                  /**< Error \f$E_\infty\f$                      */
  double t;                    /**<                                           */
  double *w;                   /**< The quadrature weights                    */
  double *x;                   /**< The quadrature nodes                      *
                                    \f$\left(w_k\right)_{k=0}^{M}\f$          */
  complex *f_hat;              /**< The spherical Fourier coefficients        */
  complex *f;                  /**< The function values                       */
  complex *f_bak;              /**< A copy of the exact function values       */
  nfsft_plan plan;             /**< NFSFT plan                                */
  nfsft_plan plan_adjoint;     /**< adjoint NFSFT plan                        */
  int i;                       /**< */
  int j;                       /**< */
  int k;                       /**< */
  int n;                       /**< */
  int d;                       /**< */
  int grid_max_theta;          /**< Number of different colatitudinal angles. */
  int grid_max_phi;            /**< Number of different colatitudinal angles. */
  int grid_max_total;          /**< The total maximum number of grid nodes.   */
  int grid_theta;              /**< The current number of different           *
                                    colatitudinal angles.                     */
  int grid_phi;                /**< The current number of different           *
                                    colatitudinal angles.                     */
  int grid_total;              /**< The total maximum number of grid nodes.   */
  double *theta;               /**< */
  double *phi;                 /**< */

  /* Read the number of testcases. */
  fscanf(stdin,"testcases=%d\n",&tc_max);
  fprintf(stdout,"%d\n",tc_max);

  /* Process each testcase. */
  for (tc = 0; tc < tc_max; tc++)
  {
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
        threshold = 1000000000000.0;
      }
    }
    else
    {
      /* TODO remove this */
      /* Set dummy values. */
      cutoff = 3;
      threshold = 1000000000000.0;
    }

    /* Read the quadrature grid type. */
    fscanf(stdin,"gridtype=%d\n",&gridtype);

    /* Read the number of repetitions. */
    fscanf(stdin,"repetitions=%d\n",&repetitions);

    /* Initialize bandwidth bound. */
    m_max = 0;

    /* Read the number of cut-off degrees. */
    fscanf(stdin,"bandwidths=%d\n",&im_max);
    fprintf(stdout,"%d\n",im_max);

    /* Allocate memory for the cut-off degrees. */
    m = (int*) malloc(im_max*sizeof(int));

    /* Read the cut-off degrees. */
    for (im = 0; im < im_max; im++)
    {
      /* Read cut-off degree. */
      fscanf(stdin,"%d\n",&m[im]);
      fprintf(stdout,"%d\n",m[im]);
      m_max = MAX(m_max,m[im]);
    }

    /* Determine the maximum number of nodes in co-latitudinal and longitudinal
     * direction. */
    switch (gridtype)
    {
      case GRID_GAUSS_LEGENDRE:
        grid_max_theta = m_max+1;
        grid_max_phi = 2*m_max+2;
        break;
      case GRID_CLENSHAW_CURTIS:
        grid_max_theta = 2*m_max+1;
        grid_max_phi = 2*m_max+2;
        break;
      default:
        grid_max_theta = m_max+1;
        grid_max_phi = 2*m_max+2;
    }
    grid_max_total = grid_max_theta * grid_max_phi;

    /* Allocate memory for data structures. */
    theta = (double*) malloc(grid_max_theta*sizeof(double));
    phi = (double*) malloc(grid_max_phi*sizeof(double));
    x = (double*) malloc(2*grid_max_total*sizeof(double));
    w = (double*) malloc(grid_max_theta*sizeof(double));
    f = (complex*) malloc(grid_max_total*sizeof(complex));
    f_bak = (complex*) malloc(grid_max_total*sizeof(complex));
    f_hat = (complex*) malloc(NFSFT_F_HAT_SIZE(m_max)*sizeof(complex));

    /* Generate the quadrature nodes. */
    /*for (l = 0; l < l_max; l++)
    {
      b[l] = drand48() - 0.5;
      eta[2*l] = drand48() - 0.5;
      eta[2*l+1] = 0.5*drand48();
    }*/

    /* Do precomputation. */
    nfsft_precompute(m_max,threshold, 0U |
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U)));

    /* Initialize error and cumulative time variable. */
    err = -1.0;
    t_avg = -1.0;

    /* Process all cut-off bandwidths. */
    for (im = 0; im < im_max; im++)
    {
      switch (gridtype)
      {
        case GRID_GAUSS_LEGENDRE:
          /* Calculate grid dimensions. */
          grid_theta = m[im] + 1;
          grid_phi = 2*m[im] + 2;
          grid_total = grid_theta*grid_phi;

          fprintf(stdout,"\n");

          /* Read quadrature weights. */
          for (k = 0; k < grid_theta; k++)
          {
            fscanf(stdin,"%le\n",&w[k]);
            fprintf(stdout,"%le\n",w[k]);
          }

          fprintf(stdout,"\n");

          /* Read grid angles theta. */
          for (k = 0; k < grid_theta; k++)
          {
            fscanf(stdin,"%le\n",&theta[k]);
            fprintf(stdout,"%le\n",theta[k]);
          }

          fprintf(stdout,"\n");

          /* Read grid angles phi. */
          for (n = 0; n < grid_phi; n++)
          {
            phi[n] = n/((double)grid_phi)-0.5;
            fprintf(stdout,"%le\n",phi[n]);
          }

          /* Generate grid nodes. */
          d = 0;
          for (k = 0; k < grid_theta; k++)
          {
            for (n = 0; n < grid_phi; n++)
            {
              x[2*d] = phi[n];
              x[2*d+1] = theta[k];
              d++;
            }
          }
          break;
        case GRID_CLENSHAW_CURTIS:
          /* Calculate grid dimensions. */
          grid_theta = 2*m[im] + 1;
          grid_phi = 2*m[im] + 2;
          grid_total = grid_theta*grid_phi;
          /* Generate quadrature nodes. */
          /*x = ...;/
          /* Generate quadrature weights. */
          /*w = ...;*/
          break;
      }

      /* Generate random function samples. */
      for (d = 0; d < grid_total; d++)
      {
        f_bak[d] = drand48() - 0.5 + I*(drand48() - 0.5);
      }

      /* Init transform plans. */
      nfsft_init_guru(&plan_adjoint,m[im],grid_total,
        ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
        ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)), cutoff);
      nfsft_init_guru(&plan,m[im],grid_total,
        ((use_nfft!=0)?(0U):(NFSFT_USE_NDFT)) |
        ((use_fpt!=0)?(0U):(NFSFT_USE_DPT)), cutoff);
      plan_adjoint.f_hat = f_hat;
      plan_adjoint.x = x;
      plan_adjoint.f = f;
      plan.f_hat = f_hat;
      plan.x = x;
      plan.f = f;
      /*nfsft_precompute_x(&plan_adjoint);
      nfsft_precompute_x(&plan);*/

      /* Initialize cumulative time variable. */
      t_avg = 0.0;
      err = 0.0;

      /* Cycle through all runs. */
      for (i = 0; i < repetitions; i++)
      {
        /* Copy exact funtion values to working array. */
//        memcpy(f,f_bak,grid_total*sizeof(complex));

        /* Initialize time measurement. */
//        t = second();

        /* Multiplication with the quadrature weights. */
/*        d = 0;
        for (k = 0; k < grid_theta; k++)
        {
          for (n = 0; n < grid_phi; n++)
          {
            f[d] *= w[k];
          }
          d++;
        }*/

        /* Check if the fast NFSFT algorithm shall be tested. */
//        if (use_nfsft != NO)
//        {
          /* Execute the adjoint NFSFT transformation. */
//          nfsft_adjoint(&plan_adjoint);
          /* Execute the NFSFT transformation. */
/*          nfsft_trafo(&plan);
        }
        else
        {*/
          /* Execute the adjoint direct NDSFT transformation. */
//          ndsft_adjoint(&plan_adjoint);
          /* Execute the direct NDSFT transformation. */
//          ndsft_trafo(&plan);
//        }

/*        t_avg += second() - t;

        err += error_l_2_complex(f, f_bak, grid_total);*/
      }

      /* Calculate average time needed. */
      t_avg = t_avg/((double)repetitions);

      /* Calculate the average error. */
      err = err/((double)repetitions);

      /* Print out the error measurements. */
      fprintf(stdout,"%e\n%e\n", t_avg, err);

      /* Finalize the NFSFT plans */
      nfsft_finalize(&plan_adjoint);
      nfsft_finalize(&plan);
    } /* for (im = 0; im < im_max; im++) - Process all cut-off
       * bandwidths.*/

    /* Delete precomputed data. */
    nfsft_forget();

    /* Free data arrays. */
    free(theta);
    free(phi);
    free(x);
    free(w);
    free(f);
    free(f_bak);
    free(f_hat);

    /* Free memory for cut-off bandwidths. */
    free(m);
  } /* for (tc = 0; tc < tc_max; tc++) - Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
