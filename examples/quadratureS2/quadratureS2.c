/* $Id$
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
#include <string.h>
#include <time.h>

/* Include FFTW header. */
#include <fftw3.h>

/* Include NFFT3 library header. */
#include "nfft3.h"

/* Include NFFT 3 utilities headers. */
#include "util.h"

/** Enumeration for parameter values */
enum boolean {NO = 0, YES = 1};

/** Enumeration for quadrature grid types */
enum gridtype {GRID_GAUSS_LEGENDRE = 0, GRID_CLENSHAW_CURTIS = 1,
  GRID_HEALPIX = 2, GRID_EQUIDISTRIBUTION_7_1_11 = 3};

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
  int tc;                      /**< The index variable for testcases          */
  int tc_max;                  /**< The number of testcases                   */
  int *m;                      /**< The array containing the cut-off degrees  *
                                    \f$M\f$.                                  */
  int m_max;                   /**< The maximum cut-off degree \f$M\f$ for the*
                                    current testcase                          */
  int n_max;                   /**< Next greater power of two with respect to *
                                    m_max                                     */
  int im;                      /**< Index variable for cut-off degrees        */
  int im_max;                  /**< The maximum number of cut-off degrees     */
  int use_nfsft;               /**< Whether to use the NFSFT algorithm or not */
  int use_nfft;                /**< Whether to use the NFFT algorithm or not  */
  int use_fpt;                 /**< Whether to use the FPT algorithm or not   */
  int cutoff;                  /**< The current NFFT cut-off parameter        */
  double threshold;            /**< The current NFSFT threshold parameter     */
  int gridtype;                /**< The type of quadrature grid to be used    */
  int repetitions;             /**< The number of repetitions to be performed */
  double t_avg;                /**< The average computation time needed       */
  double err_infty;            /**< Error \f$E_\infty\f$                      */
  double err_2;                /**< Error \f$E_2\f$                           */
  double t;                    /**< A variable for time measurements          */
  double *w;                   /**< The quadrature weights                    */
  double *x;                   /**< The quadrature nodes                      *
                                    \f$\left(w_k\right)_{k=0}^{M}\f$          */
  complex *f_hat;              /**< The spherical Fourier coefficients        */
  complex *f_hat_ref;          /**< The spherical Fourier coefficients        */
  complex *f;                  /**< The function values                       */
  complex *f_ref;              /**< A copy of the exact function values       */
  nfsft_plan plan;             /**< The NFSFT plan                            */
  nfsft_plan plan_gen;         /**< The NFSFT plan                            */
  int i;                       /**< A loop variable                           */
  int j;                       /**< A loop variable                           */
  int k;                       /**< A loop variable                           */
  int n;                       /**< A loop variable                           */
  int d;                       /**< A loop variable                           */
  int grid_theta;              /**< The current number of different           *
                                    colatitudinal angles.                     */
  int grid_phi;                /**< The current number of different           *
                                    longitudinal angles.                      */
  int grid_total;              /**< The total number of grid nodes.           */
  int grid_max_theta;          /**< Number of different colatitudinal angles. */
  int grid_max_phi;            /**< Number of different longitudinal angles.  */
  int grid_max_total;          /**< The total maximum number of grid nodes.   */
  double *theta;               /**< */
  double *phi;                 /**< */
  fftw_plan fplan;             /**< */
  int nside;                   /**< */
  int d2;
  int gamma;
  double thetai;
  double gammai;
  FILE *file;

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

    /* Read the quadrature grid type. */
    fscanf(stdin,"gridtype=%d\n",&gridtype);
    fprintf(stdout,"%d\n",gridtype);

    /* Read the number of repetitions. */
    fscanf(stdin,"repetitions=%d\n",&repetitions);
    fprintf(stdout,"%d\n",repetitions);

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
    f_ref = (complex*) malloc(grid_max_total*sizeof(complex));
    f_hat = (complex*) malloc(NFSFT_F_HAT_SIZE(m_max)*sizeof(complex));
    f_hat_ref = (complex*) malloc(NFSFT_F_HAT_SIZE(m_max)*sizeof(complex));

    /* Do precomputation. */
    nfsft_precompute(m_max,threshold, 0U |
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U)));

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

          /* Read quadrature weights. */
          for (k = 0; k < grid_theta; k++)
          {
            fscanf(stdin,"%le\n",&w[k]);
            w[k] *= 2.0*PI/(2*m[im]+2);
          }

          /* Read angles theta. */
          for (k = 0; k < grid_theta; k++)
          {
            fscanf(stdin,"%le\n",&theta[k]);
          }

          /* Generate the grid angles phi. */
          for (n = 0; n < grid_phi; n++)
          {
            phi[n] = n/((double)grid_phi);
            phi[n] -= ((phi[n]>=0.5)?(1.0):(0.0));
          }

          /* Generate the grid. */
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
          d = 0;
          for (k = 0; k < grid_theta; k++)
          {
            for (n = 0; n < grid_phi; n++)
            {
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
          for (n = 0; n < 2*m[im]+2; n++)
          {
            phi[n] = n/((double)2*m[im]+2)-0.5;
          }

          for (k = 0; k < 2*m[im]+1; k++)
          {
            theta[k] = k/((double)4*m[im]);
          }

          /* Generate quadrature weights. */
          fplan = fftw_plan_r2r_1d(m[im]+1, w, w, FFTW_REDFT00, 0U);
          for (k = 0; k < m[im]+1; k++)
          {
            w[k] = -2.0/(4*k*k-1);
          }
          fftw_execute(fplan);
          w[0] *= 0.5;

          for (k = 0; k < m[im]+1; k++)
          {
            w[k] *= 1/((double)2*m[im]);
            w[2*m[im]-k] = w[k];
            w[k] *= 2.0*PI/(2*m[im]+2);
            w[2*m[im]-k] *= 2.0*PI/(2*m[im]+2);
          }
          fftw_destroy_plan(fplan);

          d = 0;
          for (k = 0; k < grid_theta; k++)
          {
            for (n = 0; n < grid_phi; n++)
            {
              x[2*d] = phi[n];
              x[2*d+1] = theta[k];
              /*fprintf(stdout,"x[%d] = %le, x[%d] = %le\n",2*d,x[2*d],2*d+1,
                x[2*d+1]);*/
              d++;
            }
          }
          break;

        case GRID_HEALPIX:
          /* Calculate grid dimensions. */
          nside = next_power_of_2(ceil((2.0*m[im])));
          grid_theta = 1;
          grid_phi = 12*nside*nside;
          grid_total = grid_theta*grid_phi;

          d = 0;
          for (k = 1; k <= nside-1; k++)
          {
            for (n = 0; n <= 4*k-1; n++)
            {
              x[2*d+1] = 1 - (k*k)/((double)(3.0*nside*nside));
              x[2*d] =  ((n+0.5)/(4*k));
              x[2*d] -= (x[2*d]>=0.5)?(1.0):(0.0);
              d++;
            }
          }

          d2 = d-1;

          for (k = nside; k <= 3*nside; k++)
          {
            for (n = 0; n <= 4*nside-1; n++)
            {
              x[2*d+1] = 2.0/(3*nside)*(2*nside-k);
              x[2*d] = (n+((k%2==0)?(0.5):(0.0)))/(4*nside);
              x[2*d] -= (x[2*d]>=0.5)?(1.0):(0.0);
              d++;
            }
          }

          for (k = 1; k <= nside-1; k++)
          {
            for (n = 0; n <= 4*k-1; n++)
            {
              x[2*d+1] = -x[2*d2+1];
              x[2*d] =  x[2*d2];
              d++;
              d2--;
            }
          }

          for (d = 0; d < grid_total; d++)
          {
            x[2*d+1] = acos(x[2*d+1])/(2.0*PI);
          }

          w[0] = (4.0*PI)/(grid_total);
          break;

        case GRID_EQUIDISTRIBUTION_7_1_11:
          /* Calculate grid dimensions. */
          gamma = 2*m[im];
          grid_theta = 1;
          grid_phi = 2+4*((int)floor((gamma+1)/2.0))*((int)floor(gamma/2.0));
          grid_total = grid_theta*grid_phi;
          x[0] = 0.0;
          x[1] = 0.0;
          d = 1;
          for (k = 1; k <= gamma-1; k++)
          {
            thetai = (k*PI)/gamma;
            gammai = ((k<=(gamma/2.0))?(4*k):(4*(gamma-k)));
            for(n = 1; n <= gammai; n++)
            {
              x[2*d+1] = thetai/(2.0*PI);
              x[2*d] = ((n-0.5)*((2.0*PI)/gammai))/(2.0*PI);
              x[2*d] -= (x[2*d]>=0.5)?(1.0):(0.0);
              d++;
            }
          }
          x[2*d+1] = 0.5;
          x[2*d] = 0.0;
          w[0] = (4.0*PI)/(grid_total);
          break;

        default:
          break;
      }

      /* Generate random function samples. */
      fprintf(stderr,"quadratureS2: flags = %d\n",NFSFT_NORMALIZED | ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
                      ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)));
      nfsft_init_guru(&plan_gen,m[im],grid_total,
        NFSFT_NORMALIZED | ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
                      ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)), cutoff);
      plan_gen.f_hat = f_hat_ref;
      plan_gen.x = x;
      plan_gen.f = f_ref;
      nfsft_precompute_x(&plan_gen);
      for (k = 0; k < plan_gen.N_total; k++)
      {
        f_hat_ref[k] = 0.0;
      }
      for (k = 0; k < m[im]; k++)
      {
        for (n = -k; n <= k; n++)
        {
          f_hat_ref[NFSFT_INDEX(k,n,&plan_gen)] = drand48()-0.5 + I*(drand48()-0.5);
        }
      }
        if (use_nfsft != NO)
        {
          /* Execute the NFSFT transformation. */
          nfsft_trafo(&plan_gen);
        }
        else
        {
          /* Execute the direct NDSFT transformation. */
          ndsft_trafo(&plan_gen);
        }
      nfsft_finalize(&plan_gen);

      /* Init transform plans. */
      nfsft_init_guru(&plan,m[im],grid_total, NFSFT_NORMALIZED | ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
                      ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)), cutoff);
      plan.f_hat = f_hat;
      plan.x = x;
      plan.f = f;
      nfsft_precompute_x(&plan);

      /* Initialize cumulative time variable. */
      t_avg = 0.0;
      err_infty = 0.0;
      err_2 = 0.0;

      /* Cycle through all runs. */
      for (i = 0; i < repetitions; i++)
      {
        /* Copy exact funtion values to working array. */
        memcpy(f,f_ref,grid_total*sizeof(complex));

        /* Initialize time measurement. */
        t = second();

        /* Multiplication with the quadrature weights. */
        /*fprintf(stderr,"\n");*/
        d = 0;
        for (k = 0; k < grid_theta; k++)
        {
          for (n = 0; n < grid_phi; n++)
          {
            /*fprintf(stderr,"f_bak[%d] = %le + I*%le,\t f[%d] = %le + I*%le,  \t w[%d] = %le\n",
            d,creal(f_bak[d]),cimag(f_bak[d]),d,creal(f[d]),cimag(f[d]),k,
            w[k]);*/
            f[d] *= w[k];
            d++;
          }
        }

        /*fprintf(stderr,"\n");
        d = 0;
        for (d = 0; d < grid_total; d++)
        {
          fprintf(stderr,"f[%d] = %le + I*%le, theta[%d] = %le, phi[%d] = %le\n",
                  d,creal(f[d]),cimag(f[d]),d,x[2*d+1],d,x[2*d]);
        }*/

        /* Check if the fast NFSFT algorithm shall be tested. */
        if (use_nfsft != NO)
        {
          /* Execute the adjoint NFSFT transformation. */
          nfsft_adjoint(&plan);
        }
        else
        {
          /* Execute the adjoint direct NDSFT transformation. */
          ndsft_adjoint(&plan);
        }

        /* Multiplication with the Fourier-Legendre coefficients. */
        /*for (k = 0; k <= m[im]; k++)
        {
          for (n = -k; n <= k; n++)
          {
            fprintf(stderr,"f_hat[%d,%d] = %le\t + I*%le\n",k,n,
                    creal(f_hat[NFSFT_INDEX(k,n,&plan_adjoint)]),
                    cimag(f_hat[NFSFT_INDEX(k,n,&plan_adjoint)]));
          }
        }*/

        if (use_nfsft != NO)
        {
          /* Execute the NFSFT transformation. */
          nfsft_trafo(&plan);
        }
        else
        {
          /* Execute the direct NDSFT transformation. */
          ndsft_trafo(&plan);
        }

        t_avg += second() - t;

        err_infty += error_l_infty_complex(f, f_ref, grid_total);
        err_2 += error_l_2_complex(f, f_ref, grid_total);

        /*for (d = 0; d < grid_total; d++)
        {
          fprintf(stderr,"f_bak[%d] = %le + I*%le,\t f[%d] = %le + I*%le\n",
                  d,creal(f_bak[d]),cimag(f_bak[d]),d,creal(f[d]),cimag(f[d]));
        }*/
      }

      /* Calculate average time needed. */
      t_avg = t_avg/((double)repetitions);

      /* Calculate the average error. */
      err_infty = err_infty/((double)repetitions);

      /* Calculate the average error. */
      err_2 = err_2/((double)repetitions);

      /* Print out the error measurements. */
      fprintf(stdout,"%4d: %e\t%e\t%e\n", m[im], t_avg, err_infty, err_2);

      /* Finalize the NFSFT plans */
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
    free(f_ref);
    free(f_hat);

    /* Free memory for cut-off bandwidths. */
    free(m);
  } /* for (tc = 0; tc < tc_max; tc++) - Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
