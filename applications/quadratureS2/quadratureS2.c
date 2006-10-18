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
/* extern double drand48 (void) __THROW; */
#include <stdio.h>
#include <string.h>
#include <time.h>

/* Include FFTW header. */
#include <fftw3.h>

/* Include NFFT 3 library header. */
#include "nfft3.h"

/* Include NFFT 3 utilities headers. */
#include "util.h"

/** Enumeration for parameter values */
enum boolean {NO = 0, YES = 1};

/** Enumeration for test modes */
enum testtype {ERROR = 0, TIMING = 1};

/** Enumeration for quadrature grid types */
enum gridtype {GRID_GAUSS_LEGENDRE = 0, GRID_CLENSHAW_CURTIS = 1,
  GRID_HEALPIX = 2, GRID_EQUIDISTRIBUTION = 3, GRID_EQUIDISTRIBUTION_UNIFORM = 4};

/** Enumeration for test functions */
enum functiontype {FUNCTION_RANDOM_BANDLIMITED = 0, FUNCTION_F1 = 1,
  FUNCTION_F2 = 2, FUNCTION_F3 = 3, FUNCTION_F4 = 4, FUNCTION_F5 = 5,
  FUNCTION_F6 = 6};

/** TODO Add comment here */
enum modes {USE_GRID = 0, RANDOM = 1};

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

  int *NQ;                     /**< The array containing the cut-off degrees  *
                                    \f$N\f$                                   */
  int NQ_max;                  /**< The maximum cut-off degree \f$N\f$ for the*
                                    current testcase                          */
  int *SQ;                     /**< The array containing the grid size
                                    parameters                                */
  int SQ_max;                  /**< The maximum grid size parameter           */
  int *RQ;                     /**< The array containing the grid size
                                    parameters                                */
  int iNQ;                     /**< Index variable for cut-off degrees        */
  int iNQ_max;                 /**< The maximum number of cut-off degrees     */
  int testfunction;            /**< The testfunction                          */
  int N;                       /**< The test function's bandwidth             */

  int use_nfsft;               /**< Whether to use the NFSFT algorithm or not */
  int use_nfft;                /**< Whether to use the NFFT algorithm or not  */
  int use_fpt;                 /**< Whether to use the FPT algorithm or not   */
  int cutoff;                  /**< The current NFFT cut-off parameter        */
  double threshold;            /**< The current NFSFT threshold parameter     */

  int gridtype;                /**< The type of quadrature grid to be used    */
  int repetitions;             /**< The number of repetitions to be performed */
  int mode;                    /**< The number of repetitions to be performed */

  double *w;                   /**< The quadrature weights                    */
  double *x_grid;              /**< The quadrature nodes                      */
  double *x_compare;           /**< The quadrature nodes                      */
  double complex *f_grid;             /**< The reference function values             */
  double complex *f_compare;          /**< The function values                       */
  double complex *f;                  /**< The function values                       */
  double complex *f_hat_gen;         /**< The reference spherical Fourier           *
                                    coefficients                              */
  double complex *f_hat;              /**< The spherical Fourier coefficients        */

  nfsft_plan plan_adjoint;     /**< The NFSFT plan                            */
  nfsft_plan plan;             /**< The NFSFT plan                            */
  nfsft_plan plan_gen;         /**< The NFSFT plan                            */

  double t;                    /**< The computation time needed for a single  *
                                    run                                       */
  double t_avg;                /**< The average computation time needed       */
  double err_infty_avg;        /**< The average error \f$E_\infty\f$          */
  double err_2_avg;            /**< The average error \f$E_2\f$               */

  int i;                       /**< A loop variable                           */
  int j;                       /**< A loop variable                           */
  int k;                       /**< A loop variable                           */
  int n;                       /**< A loop variable                           */
  int d;                       /**< A loop variable                           */

  int m_theta;                 /**< The current number of different           *
                                    colatitudinal angles (for grids)          */
  int m_phi;                   /**< The current number of different           *
                                    longitudinal angles (for grids).          */
  int m_total;                 /**< The total number nodes.                   */
  int m_theta_max;             /**< The maximum number of different           *
                                    colatitudinal angles (for grids).         */
  int m_phi_max;               /**< The maximum number of different           *
                                    longitudinal angles (for grids).          */
  int m_total_max;          /**< The total maximum number of grid nodes.   */
  double *theta;               /**< An array for saving the angles theta of a *
                                    grid                                      */
  double *phi;                 /**< An array for saving the angles phi of a   *
                                    grid                                      */
  fftw_plan fplan;             /**< An FFTW plan for computing Clenshaw-Curtis
                                    quadrature weights                        */
  //int nside;                   /**< The size parameter for the HEALPix grid   */
  int d2;
  int gamma;
  double thetai;
  double gammai;
  FILE *file;
  int M;
  double theta_s;
  double x1,x2,x3,temp;
  int m_compare;
  nfsft_plan *plan_adjoint_ptr;
  nfsft_plan *plan_ptr;
  double *w_temp;
  int testmode;

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

    /* Read the testmode type. */
    fscanf(stdin,"testmode=%d\n",&testmode);
    fprintf(stdout,"%d\n",testmode);

    if (testmode == ERROR)
    {
      /* Read the quadrature grid type. */
      fscanf(stdin,"gridtype=%d\n",&gridtype);
      fprintf(stdout,"%d\n",gridtype);

      /* Read the test function. */
      fscanf(stdin,"testfunction=%d\n",&testfunction);
      fprintf(stdout,"%d\n",testfunction);

      /* Check if random bandlimited function has been chosen. */
      if (testfunction == FUNCTION_RANDOM_BANDLIMITED)
      {
        /* Read the bandwidht. */
        fscanf(stdin,"bandlimit=%d\n",&N);
        fprintf(stdout,"%d\n",N);
      }
      else
      {
        N = 1;
      }

      /* Read the number of repetitions. */
      fscanf(stdin,"repetitions=%d\n",&repetitions);
      fprintf(stdout,"%d\n",repetitions);

      fscanf(stdin,"mode=%d\n",&mode);
      fprintf(stdout,"%d\n",mode);

      if (mode == RANDOM)
      {
        /* Read the bandwidht. */
        fscanf(stdin,"points=%d\n",&m_compare);
        fprintf(stdout,"%d\n",m_compare);
        x_compare = (double*) malloc(2*m_compare*sizeof(double));
        d = 0;
        while (d < m_compare)
        {
          x1 = 2.0*drand48() - 1.0;
          x2 = 2.0*drand48() - 1.0;
          x3 = 2.0*drand48() - 1.0;
          temp = sqrt(x1*x1+x2*x2+x3*x3);
          if (temp <= 1)
          {
            x_compare[2*d+1] = acos(x3);
            if (x_compare[2*d+1] == 0 || x_compare[2*d+1] == PI)
            {
              x_compare[2*d] = 0.0;
            }
            else
            {
              x_compare[2*d] = atan2(x2/sin(x_compare[2*d+1]),x1/sin(x_compare[2*d+1]));
            }
            x_compare[2*d] *= 1.0/(2.0*PI);
            x_compare[2*d+1] *= 1.0/(2.0*PI);
            d++;
          }
        }
        f_compare = (double complex*) malloc(m_compare*sizeof(double complex));
        f = (double complex*) malloc(m_compare*sizeof(double complex));
      }
    }

    /* Initialize maximum cut-off degree and grid size parameter. */
    NQ_max = 0;
    SQ_max = 0;

    /* Read the number of cut-off degrees. */
    fscanf(stdin,"bandwidths=%d\n",&iNQ_max);
    fprintf(stdout,"%d\n",iNQ_max);

    /* Allocate memory for the cut-off degrees and grid size parameters. */
    NQ = (int*) malloc(iNQ_max*sizeof(int));
    SQ = (int*) malloc(iNQ_max*sizeof(int));
    if (testmode == TIMING)
    {
      RQ = (int*) malloc(iNQ_max*sizeof(int));
    }

    /* Read the cut-off degrees and grid size parameters. */
    for (iNQ = 0; iNQ < iNQ_max; iNQ++)
    {
      if (testmode == TIMING)
      {
        /* Read cut-off degree and grid size parameter. */
        fscanf(stdin,"%d %d %d\n",&NQ[iNQ],&SQ[iNQ],&RQ[iNQ]);
        fprintf(stdout,"%d %d %d\n",NQ[iNQ],SQ[iNQ],RQ[iNQ]);
        NQ_max = NFFT_MAX(NQ_max,NQ[iNQ]);
        SQ_max = NFFT_MAX(SQ_max,SQ[iNQ]);
      }
      else
      {
        /* Read cut-off degree and grid size parameter. */
        fscanf(stdin,"%d %d\n",&NQ[iNQ],&SQ[iNQ]);
        fprintf(stdout,"%d %d\n",NQ[iNQ],SQ[iNQ]);
        NQ_max = NFFT_MAX(NQ_max,NQ[iNQ]);
        SQ_max = NFFT_MAX(SQ_max,SQ[iNQ]);
      }
    }

    /* Do precomputation. */
    //fprintf(stderr,"NFSFT Precomputation\n");
    //fflush(stderr);
    nfsft_precompute(NQ_max, threshold,
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U)), 0U);

    if (testmode == TIMING)
    {
      /* Allocate data structures. */
      f_hat = (double complex*) malloc(NFSFT_F_HAT_SIZE(NQ_max)*sizeof(double complex));
      f = (double complex*) malloc(SQ_max*sizeof(double complex));
      x_grid = (double*) malloc(2*SQ_max*sizeof(double));
      for (d = 0; d < SQ_max; d++)
      {
        f[d] = drand48()-0.5 + I*(drand48()-0.5);
        x_grid[2*d] = drand48() - 0.5;
        x_grid[2*d+1] = drand48() * 0.5;
      }
    }

    //fprintf(stderr,"Entering loop\n");
    //fflush(stderr);
    /* Process all cut-off bandwidths. */
    for (iNQ = 0; iNQ < iNQ_max; iNQ++)
    {
      if (testmode == TIMING)
      {
        nfsft_init_guru(&plan,NQ[iNQ],SQ[iNQ], NFSFT_NORMALIZED |
          ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
          ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
          PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFTW_MEASURE | FFT_OUT_OF_PLACE,
          cutoff);

        plan.f_hat = f_hat;
        plan.x = x_grid;
        plan.f = f;

        nfsft_precompute_x(&plan);

        t_avg = 0.0;

        for (i = 0; i < RQ[iNQ]; i++)
        {
          t = nfft_second();

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

          t_avg += nfft_second() - t;
        }

        t_avg = t_avg/((double)RQ[iNQ]);

        nfsft_finalize(&plan);

        fprintf(stdout,"%+le\n", t_avg);
        fprintf(stderr,"%d: %4d %4d %+le\n", tc, NQ[iNQ], SQ[iNQ], t_avg);
      }
      else
      {
        /* Determine the maximum number of nodes. */
        switch (gridtype)
        {
          case GRID_GAUSS_LEGENDRE:
            /* Calculate grid dimensions. */
            m_theta = SQ[iNQ] + 1;
            m_phi = 2*SQ[iNQ] + 2;
            m_total = m_theta*m_phi;
            break;
          case GRID_CLENSHAW_CURTIS:
            /* Calculate grid dimensions. */
            m_theta = 2*SQ[iNQ] + 1;
            m_phi = 2*SQ[iNQ] + 2;
            m_total = m_theta*m_phi;
            break;
          case GRID_HEALPIX:
            m_theta = 1;
            m_phi = 12*SQ[iNQ]*SQ[iNQ];
            m_total = m_theta * m_phi;
            //fprintf("HEALPix: SQ = %d, m_theta = %d, m_phi= %d, m");
            break;
          case GRID_EQUIDISTRIBUTION:
          case GRID_EQUIDISTRIBUTION_UNIFORM:
            m_theta = 2;
            //fprintf(stderr,"ed: m_theta = %d\n",m_theta);
            for (k = 1; k < SQ[iNQ]; k++)
            {
              m_theta += (int)floor((2*PI)/acos((cos(PI/(double)SQ[iNQ])-
                cos(k*PI/(double)SQ[iNQ])*cos(k*PI/(double)SQ[iNQ]))/
                (sin(k*PI/(double)SQ[iNQ])*sin(k*PI/(double)SQ[iNQ]))));
              //fprintf(stderr,"ed: m_theta = %d\n",m_theta);
            }
            //fprintf(stderr,"ed: m_theta final = %d\n",m_theta);
            m_phi = 1;
            m_total = m_theta * m_phi;
            break;
        }

        /* Allocate memory for data structures. */
        w = (double*) malloc(m_theta*sizeof(double));
        x_grid = (double*) malloc(2*m_total*sizeof(double));

        //fprintf(stderr,"NQ = %d\n",NQ[iNQ]);
        //fflush(stderr);
        switch (gridtype)
        {
          case GRID_GAUSS_LEGENDRE:
            //fprintf(stderr,"Generating grid for NQ = %d, SQ = %d\n",NQ[iNQ],SQ[iNQ]);
            //fflush(stderr);

            /* Read quadrature weights. */
            for (k = 0; k < m_theta; k++)
            {
              fscanf(stdin,"%le\n",&w[k]);
              w[k] *= (2.0*PI)/((double)m_phi);
            }

            //fprintf(stderr,"Allocating theta and phi\n");
            //fflush(stderr);
            /* Allocate memory to store the grid's angles. */
            theta = (double*) malloc(m_theta*sizeof(double));
            phi = (double*) malloc(m_phi*sizeof(double));

            //if (theta == NULL || phi == NULL)
            //{
              //fprintf(stderr,"Couldn't allocate theta and phi\n");
              //fflush(stderr);
            //}


            /* Read angles theta. */
            for (k = 0; k < m_theta; k++)
            {
              fscanf(stdin,"%le\n",&theta[k]);
            }

            /* Generate the grid angles phi. */
            for (n = 0; n < m_phi; n++)
            {
              phi[n] = n/((double)m_phi);
              phi[n] -= ((phi[n]>=0.5)?(1.0):(0.0));
            }

            //fprintf(stderr,"Generating grid nodes\n");
            //fflush(stderr);

            /* Generate the grid's nodes. */
            d = 0;
            for (k = 0; k < m_theta; k++)
            {
              for (n = 0; n < m_phi; n++)
              {
                x_grid[2*d] = phi[n];
                x_grid[2*d+1] = theta[k];
                d++;
              }
            }

            //fprintf(stderr,"Freeing theta and phi\n");
            //fflush(stderr);
            /* Free the arrays for the grid's angles. */
            free(theta);
            free(phi);

            break;

          case GRID_CLENSHAW_CURTIS:

            /* Allocate memory to store the grid's angles. */
            theta = (double*) malloc(m_theta*sizeof(double));
            phi = (double*) malloc(m_phi*sizeof(double));

            /* Generate the grid angles theta. */
            for (k = 0; k < m_theta; k++)
            {
              theta[k] = k/((double)2*(m_theta-1));
            }

            /* Generate the grid angles phi. */
            for (n = 0; n < m_phi; n++)
            {
              phi[n] = n/((double)m_phi);
              phi[n] -= ((phi[n]>=0.5)?(1.0):(0.0));
            }

            /* Generate quadrature weights. */
            fplan = fftw_plan_r2r_1d(SQ[iNQ]+1, w, w, FFTW_REDFT00, 0U);
            for (k = 0; k < SQ[iNQ]+1; k++)
            {
              w[k] = -2.0/(4*k*k-1);
            }
            fftw_execute(fplan);
            w[0] *= 0.5;

            for (k = 0; k < SQ[iNQ]+1; k++)
            {
              w[k] *= (2.0*PI)/((double)(m_theta-1)*m_phi);
              w[m_theta-1-k] = w[k];
            }
            fftw_destroy_plan(fplan);

            /* Generate the grid's nodes. */
            d = 0;
            for (k = 0; k < m_theta; k++)
            {
              for (n = 0; n < m_phi; n++)
              {
                x_grid[2*d] = phi[n];
                x_grid[2*d+1] = theta[k];
                d++;
              }
            }

            /* Free the arrays for the grid's angles. */
            free(theta);
            free(phi);

            break;

          case GRID_HEALPIX:

            d = 0;
            for (k = 1; k <= SQ[iNQ]-1; k++)
            {
              for (n = 0; n <= 4*k-1; n++)
              {
                x_grid[2*d+1] = 1 - (k*k)/((double)(3.0*SQ[iNQ]*SQ[iNQ]));
                x_grid[2*d] =  ((n+0.5)/(4*k));
                x_grid[2*d] -= (x_grid[2*d]>=0.5)?(1.0):(0.0);
                d++;
              }
            }

            d2 = d-1;

            for (k = SQ[iNQ]; k <= 3*SQ[iNQ]; k++)
            {
              for (n = 0; n <= 4*SQ[iNQ]-1; n++)
              {
                x_grid[2*d+1] = 2.0/(3*SQ[iNQ])*(2*SQ[iNQ]-k);
                x_grid[2*d] = (n+((k%2==0)?(0.5):(0.0)))/(4*SQ[iNQ]);
                x_grid[2*d] -= (x_grid[2*d]>=0.5)?(1.0):(0.0);
                d++;
              }
            }

            for (k = 1; k <= SQ[iNQ]-1; k++)
            {
              for (n = 0; n <= 4*k-1; n++)
              {
                x_grid[2*d+1] = -x_grid[2*d2+1];
                x_grid[2*d] =  x_grid[2*d2];
                d++;
                d2--;
              }
            }

            for (d = 0; d < m_total; d++)
            {
              x_grid[2*d+1] = acos(x_grid[2*d+1])/(2.0*PI);
            }

            w[0] = (4.0*PI)/(m_total);
            break;

          case GRID_EQUIDISTRIBUTION:
          case GRID_EQUIDISTRIBUTION_UNIFORM:
            /* TODO Compute the weights. */

            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              w_temp = (double*) malloc((SQ[iNQ]+1)*sizeof(double));
              fplan = fftw_plan_r2r_1d(SQ[iNQ]/2+1, w_temp, w_temp, FFTW_REDFT00, 0U);
              for (k = 0; k < SQ[iNQ]/2+1; k++)
              {
                w_temp[k] = -2.0/(4*k*k-1);
              }
              fftw_execute(fplan);
              w_temp[0] *= 0.5;

              for (k = 0; k < SQ[iNQ]/2+1; k++)
              {
                w_temp[k] *= (2.0*PI)/((double)(SQ[iNQ]));
                w_temp[SQ[iNQ]-k] = w_temp[k];
              }
              fftw_destroy_plan(fplan);
            }

            d = 0;
            x_grid[2*d] = -0.5;
            x_grid[2*d+1] = 0.0;
            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              w[d] = w_temp[0];
            }
            else
            {
              w[d] = (4.0*PI)/(m_total);
            }
            d = 1;
            x_grid[2*d] = -0.5;
            x_grid[2*d+1] = 0.5;
            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              w[d] = w_temp[SQ[iNQ]];
            }
            else
            {
              w[d] = (4.0*PI)/(m_total);
            }
            d = 2;

            for (k = 1; k < SQ[iNQ]; k++)
            {
              theta_s = (double)k*PI/(double)SQ[iNQ];
              M = (int)floor((2.0*PI)/acos((cos(PI/(double)SQ[iNQ])-
                cos(theta_s)*cos(theta_s))/(sin(theta_s)*sin(theta_s))));

              for (n = 0; n < M; n++)
              {
                x_grid[2*d] = (n + 0.5)/M;
                x_grid[2*d] -= (x_grid[2*d]>=0.5)?(1.0):(0.0);
                x_grid[2*d+1] = theta_s/(2.0*PI);
                if (gridtype == GRID_EQUIDISTRIBUTION)
                {
                  w[d] = w_temp[k]/((double)(M));
                }
                else
                {
                  w[d] = (4.0*PI)/(m_total);
                }
                d++;
              }
            }

            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              free(w_temp);
            }
            break;

          default:
            break;
        }

        /* Allocate memory for grid values. */
        f_grid = (double complex*) malloc(m_total*sizeof(double complex));

        if (mode == RANDOM)
        {
        }
        else
        {
          m_compare = m_total;
          f_compare = (double complex*) malloc(m_compare*sizeof(double complex));
          x_compare = x_grid;
          f = f_grid;
        }

        //fprintf(stderr,"Generating test function\n");
        //fflush(stderr);
        switch (testfunction)
        {
          case FUNCTION_RANDOM_BANDLIMITED:
            f_hat_gen = (double complex*) malloc(NFSFT_F_HAT_SIZE(N)*sizeof(double complex));
            //fprintf(stderr,"Generating random test function\n");
            //fflush(stderr);
            /* Generate random function samples by sampling a bandlimited
             * function. */
            nfsft_init_guru(&plan_gen,N,m_total, NFSFT_NORMALIZED |
              ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
              ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
              ((N>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT |
              FFT_OUT_OF_PLACE, cutoff);

            plan_gen.f_hat = f_hat_gen;
            plan_gen.x = x_grid;
            plan_gen.f = f_grid;

            nfsft_precompute_x(&plan_gen);

            for (k = 0; k < plan_gen.N_total; k++)
            {
              f_hat_gen[k] = 0.0;
            }

            for (k = 0; k <= N; k++)
            {
              for (n = -k; n <= k; n++)
              {
                f_hat_gen[NFSFT_INDEX(k,n,&plan_gen)] =
                drand48()-0.5 + I*(drand48()-0.5);
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

            if (mode == RANDOM)
            {
              nfsft_init_guru(&plan_gen,N,m_compare, NFSFT_NORMALIZED |
                ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
                ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
                ((N>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT |
                FFT_OUT_OF_PLACE, cutoff);

              plan_gen.f_hat = f_hat_gen;
              plan_gen.x = x_compare;
              plan_gen.f = f_compare;

              nfsft_precompute_x(&plan_gen);

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
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }

            free(f_hat_gen);

            break;

          case FUNCTION_F1:
            for (d = 0; d < m_total; d++)
            {
              x1 = sin(x_grid[2*d+1]*2.0*PI)*cos(x_grid[2*d]*2.0*PI);
              x2 = sin(x_grid[2*d+1]*2.0*PI)*sin(x_grid[2*d]*2.0*PI);
              x3 = cos(x_grid[2*d+1]*2.0*PI);
              f_grid[d] = x1*x2*x3;
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = sin(x_compare[2*d+1]*2.0*PI)*cos(x_compare[2*d]*2.0*PI);
                x2 = sin(x_compare[2*d+1]*2.0*PI)*sin(x_compare[2*d]*2.0*PI);
                x3 = cos(x_compare[2*d+1]*2.0*PI);
                f_compare[d] = x1*x2*x3;
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }
            break;
          case FUNCTION_F2:
            for (d = 0; d < m_total; d++)
            {
              x1 = sin(x_grid[2*d+1]*2.0*PI)*cos(x_grid[2*d]*2.0*PI);
              x2 = sin(x_grid[2*d+1]*2.0*PI)*sin(x_grid[2*d]*2.0*PI);
              x3 = cos(x_grid[2*d+1]*2.0*PI);
              f_grid[d] = 0.1*exp(x1+x2+x3);
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = sin(x_compare[2*d+1]*2.0*PI)*cos(x_compare[2*d]*2.0*PI);
                x2 = sin(x_compare[2*d+1]*2.0*PI)*sin(x_compare[2*d]*2.0*PI);
                x3 = cos(x_compare[2*d+1]*2.0*PI);
                f_compare[d] = 0.1*exp(x1+x2+x3);
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }
            break;
          case FUNCTION_F3:
            for (d = 0; d < m_total; d++)
            {
              x1 = sin(x_grid[2*d+1]*2.0*PI)*cos(x_grid[2*d]*2.0*PI);
              x2 = sin(x_grid[2*d+1]*2.0*PI)*sin(x_grid[2*d]*2.0*PI);
              x3 = cos(x_grid[2*d+1]*2.0*PI);
              temp = sqrt(x1*x1)+sqrt(x2*x2)+sqrt(x3*x3);
              f_grid[d] = 0.1*temp;
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = sin(x_compare[2*d+1]*2.0*PI)*cos(x_compare[2*d]*2.0*PI);
                x2 = sin(x_compare[2*d+1]*2.0*PI)*sin(x_compare[2*d]*2.0*PI);
                x3 = cos(x_compare[2*d+1]*2.0*PI);
                temp = sqrt(x1*x1)+sqrt(x2*x2)+sqrt(x3*x3);
                f_compare[d] = 0.1*temp;
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }
            break;
          case FUNCTION_F4:
            for (d = 0; d < m_total; d++)
            {
              x1 = sin(x_grid[2*d+1]*2.0*PI)*cos(x_grid[2*d]*2.0*PI);
              x2 = sin(x_grid[2*d+1]*2.0*PI)*sin(x_grid[2*d]*2.0*PI);
              x3 = cos(x_grid[2*d+1]*2.0*PI);
              temp = sqrt(x1*x1)+sqrt(x2*x2)+sqrt(x3*x3);
              f_grid[d] = 1.0/(temp);
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = sin(x_compare[2*d+1]*2.0*PI)*cos(x_compare[2*d]*2.0*PI);
                x2 = sin(x_compare[2*d+1]*2.0*PI)*sin(x_compare[2*d]*2.0*PI);
                x3 = cos(x_compare[2*d+1]*2.0*PI);
                temp = sqrt(x1*x1)+sqrt(x2*x2)+sqrt(x3*x3);
                f_compare[d] = 1.0/(temp);
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }
            break;
          case FUNCTION_F5:
            for (d = 0; d < m_total; d++)
            {
              x1 = sin(x_grid[2*d+1]*2.0*PI)*cos(x_grid[2*d]*2.0*PI);
              x2 = sin(x_grid[2*d+1]*2.0*PI)*sin(x_grid[2*d]*2.0*PI);
              x3 = cos(x_grid[2*d+1]*2.0*PI);
              temp = sqrt(x1*x1)+sqrt(x2*x2)+sqrt(x3*x3);
              f_grid[d] = 0.1*sin(1+temp)*sin(1+temp);
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = sin(x_compare[2*d+1]*2.0*PI)*cos(x_compare[2*d]*2.0*PI);
                x2 = sin(x_compare[2*d+1]*2.0*PI)*sin(x_compare[2*d]*2.0*PI);
                x3 = cos(x_compare[2*d+1]*2.0*PI);
                temp = sqrt(x1*x1)+sqrt(x2*x2)+sqrt(x3*x3);
                f_compare[d] = 0.1*sin(1+temp)*sin(1+temp);
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }
            break;
          case FUNCTION_F6:
            for (d = 0; d < m_total; d++)
            {
              if (x_grid[2*d+1] <= 0.25)
              {
                f_grid[d] = 1.0;
              }
              else
              {
                f_grid[d] = 1.0/(sqrt(1+3*cos(2.0*PI*x_grid[2*d+1])*cos(2.0*PI*x_grid[2*d+1])));
              }
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                if (x_compare[2*d+1] <= 0.25)
                {
                  f_compare[d] = 1.0;
                }
                else
                {
                  f_compare[d] = 1.0/(sqrt(1+3*cos(2.0*PI*x_compare[2*d+1])*cos(2.0*PI*x_compare[2*d+1])));
                }
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }
            break;
          default:
            //fprintf(stderr,"Generating one function\n");
            //fflush(stderr);
            for (d = 0; d < m_total; d++)
            {
              f_grid[d] = 1.0;
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                f_compare[d] = 1.0;
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(double complex));
            }
            break;
        }

        //fprintf(stderr,"Initializing trafo\n");
        //fflush(stderr);
        /* Init transform plan. */
        nfsft_init_guru(&plan_adjoint,NQ[iNQ],m_total, NFSFT_NORMALIZED |
          ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
          ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
          ((NQ[iNQ]>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT |
          FFT_OUT_OF_PLACE, cutoff);

        plan_adjoint_ptr = &plan_adjoint;

        if (mode == RANDOM)
        {
          nfsft_init_guru(&plan,NQ[iNQ],m_compare, NFSFT_NORMALIZED |
            ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
            ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
            ((NQ[iNQ]>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT |
            FFT_OUT_OF_PLACE, cutoff);
          plan_ptr = &plan;
        }
        else
        {
          plan_ptr = &plan_adjoint;
        }

        f_hat = (double complex*) malloc(NFSFT_F_HAT_SIZE(NQ[iNQ])*sizeof(double complex));

        plan_adjoint_ptr->f_hat = f_hat;
        plan_adjoint_ptr->x = x_grid;
        plan_adjoint_ptr->f = f_grid;

        plan_ptr->f_hat = f_hat;
        plan_ptr->x = x_compare;
        plan_ptr->f = f;

        //fprintf(stderr,"Precomputing for x\n");
        //fflush(stderr);
        nfsft_precompute_x(plan_adjoint_ptr);
        if (plan_adjoint_ptr != plan_ptr)
        {
          nfsft_precompute_x(plan_ptr);
        }

        /* Initialize cumulative time variable. */
        t_avg = 0.0;
        err_infty_avg = 0.0;
        err_2_avg = 0.0;

        /* Cycle through all runs. */
        for (i = 0; i < 1/*repetitions*/; i++)
        {
          //fprintf(stderr,"Copying original values\n");
          //fflush(stderr);
          /* Copy exact funtion values to working array. */
          //memcpy(f,f_grid,m_total*sizeof(double complex));

          /* Initialize time measurement. */
          t = nfft_second();

          //fprintf(stderr,"Multiplying with quadrature weights\n");
          //fflush(stderr);
          /* Multiplication with the quadrature weights. */
          /*fprintf(stderr,"\n");*/
          d = 0;
          for (k = 0; k < m_theta; k++)
          {
            for (n = 0; n < m_phi; n++)
            {
              /*fprintf(stderr,"f_ref[%d] = %le + I*%le,\t f[%d] = %le + I*%le,  \t w[%d] = %le\n",
              d,creal(f_ref[d]),cimag(f_ref[d]),d,creal(f[d]),cimag(f[d]),k,
              w[k]);*/
              f_grid[d] *= w[k];
              d++;
            }
          }

          t_avg += nfft_second() - t;

          free(w);

          t = nfft_second();

          /*fprintf(stderr,"\n");
          d = 0;
          for (d = 0; d < grid_total; d++)
          {
            fprintf(stderr,"f[%d] = %le + I*%le, theta[%d] = %le, phi[%d] = %le\n",
                    d,creal(f[d]),cimag(f[d]),d,x[2*d+1],d,x[2*d]);
          }*/

          //fprintf(stderr,"Executing adjoint\n");
          //fflush(stderr);
          /* Check if the fast NFSFT algorithm shall be tested. */
          if (use_nfsft != NO)
          {
            /* Execute the adjoint NFSFT transformation. */
            nfsft_adjoint(plan_adjoint_ptr);
          }
          else
          {
            /* Execute the adjoint direct NDSFT transformation. */
            ndsft_adjoint(plan_adjoint_ptr);
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

          //fprintf(stderr,"Executing trafo\n");
          //fflush(stderr);
          if (use_nfsft != NO)
          {
            /* Execute the NFSFT transformation. */
            nfsft_trafo(plan_ptr);
          }
          else
          {
            /* Execute the direct NDSFT transformation. */
            ndsft_trafo(plan_ptr);
          }

          t_avg += nfft_second() - t;

          //fprintf(stderr,"Finalizing\n");
          //fflush(stderr);
          /* Finalize the NFSFT plans */
          nfsft_finalize(plan_adjoint_ptr);
          if (plan_ptr != plan_adjoint_ptr)
          {
            nfsft_finalize(plan_ptr);
          }

          /* Free data arrays. */
          free(f_hat);
          free(x_grid);

          err_infty_avg += nfft_error_l_infty_complex(f, f_compare, m_compare);
          err_2_avg += nfft_error_l_2_complex(f, f_compare, m_compare);

          free(f_grid);

          if (mode == RANDOM)
          {
          }
          else
          {
            free(f_compare);
          }

          /*for (d = 0; d < m_total; d++)
          {
            fprintf(stderr,"f_ref[%d] = %le + I*%le,\t f[%d] = %le + I*%le\n",
              d,creal(f_ref[d]),cimag(f_ref[d]),d,creal(f[d]),cimag(f[d]));
          }*/
        }

        //fprintf(stderr,"Calculating the error\n");
        //fflush(stderr);
        /* Calculate average time needed. */
        t_avg = t_avg/((double)repetitions);

        /* Calculate the average error. */
        err_infty_avg = err_infty_avg/((double)repetitions);

        /* Calculate the average error. */
        err_2_avg = err_2_avg/((double)repetitions);

        /* Print out the error measurements. */
        fprintf(stdout,"%+le %+le %+le\n", t_avg, err_infty_avg, err_2_avg);
        fprintf(stderr,"%d: %4d %4d %+le %+le %+le\n", tc, NQ[iNQ], SQ[iNQ],
          t_avg, err_infty_avg, err_2_avg);
      }
    } /* for (im = 0; im < im_max; im++) - Process all cut-off
       * bandwidths.*/
    fprintf(stderr,"\n");

    /* Delete precomputed data. */
    nfsft_forget();

    /* Free memory for cut-off bandwidths and grid size parameters. */
    free(NQ);
    free(SQ);
    if (testmode == TIMING)
    {
      free(RQ);
    }

    if (mode == RANDOM)
    {
      free(x_compare);
      free(f_compare);
      free(f);
    }

    if (testmode == TIMING)
    {
      /* Allocate data structures. */
      free(f_hat);
      free(f);
      free(x_grid);
    }

  } /* for (tc = 0; tc < tc_max; tc++) - Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
