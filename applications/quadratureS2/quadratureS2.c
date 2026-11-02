/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

/**
 * \defgroup applications_quadratureS2_test quadratureS2_test
 * \ingroup applications_quadratureS2
 * \{
 */
/* Include standard C headers. */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <complex.h>

/* Include NFFT 3 utilities headers. */

/* Include NFFT 3 library header. */
#include "nfft3.h"
/* NFFT_PRECISION_SINGLE, NFFT_PRECISION_DOUBLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

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
  NFFT_R threshold;            /**< The current NFSFT threshold parameter     */

  int gridtype;                /**< The type of quadrature grid to be used    */
  int repetitions;             /**< The number of repetitions to be performed */
  int mode;                    /**< The number of repetitions to be performed */

  NFFT_R *w;                   /**< The quadrature weights                    */
  NFFT_R *x_grid;              /**< The quadrature nodes                      */
  NFFT_R *x_compare;           /**< The quadrature nodes                      */
  NFFT_C *f_grid;             /**< The reference function values             */
  NFFT_C *f_compare;          /**< The function values                       */
  NFFT_C *f;                  /**< The function values                       */
  NFFT_C *f_hat_gen;         /**< The reference spherical Fourier           *
                                    coefficients                              */
  NFFT_C *f_hat;              /**< The spherical Fourier coefficients        */

  NFSFT(plan) plan_adjoint;     /**< The NFSFT plan                            */
  NFSFT(plan) plan;             /**< The NFSFT plan                            */
  NFSFT(plan) plan_gen;         /**< The NFSFT plan                            */

  NFFT_R t_avg;                /**< The average computation time needed       */
  NFFT_R err_infty_avg;        /**< The average error \f$E_\infty\f$          */
  NFFT_R err_2_avg;            /**< The average error \f$E_2\f$               */

  int i;                       /**< A loop variable                           */
  int k;                       /**< A loop variable                           */
  int n;                       /**< A loop variable                           */
  int d;                       /**< A loop variable                           */

  int m_theta;                 /**< The current number of different           *
                                    colatitudinal angles (for grids)          */
  int m_phi;                   /**< The current number of different           *
                                    longitudinal angles (for grids).          */
  int m_total;                 /**< The total number nodes.                   */
  NFFT_R *theta;               /**< An array for saving the angles theta of a *
                                    grid                                      */
  NFFT_R *phi;                 /**< An array for saving the angles phi of a   *
                                    grid                                      */
  FFTW(plan) fplan;             /**< An FFTW plan for computing Clenshaw-Curtis
                                    quadrature weights                        */
  //int nside;                   /**< The size parameter for the HEALPix grid   */
  int d2;
  int M;
  NFFT_R theta_s;
  NFFT_R x1,x2,x3,temp;
  int m_compare;
  NFSFT(plan) *plan_adjoint_ptr;
  NFSFT(plan) *plan_ptr;
  NFFT_R *w_temp;
  int testmode;
  NFFT_R t0, t1;

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
        fscanf(stdin,"threshold=%" NFFT__FIS__ "\n",&threshold);
        fprintf(stdout,"%" NFFT__FIS__ "\n",threshold);
      }
      else
      {
        /* TODO remove this */
        /* Initialize unused variable with dummy value. */
        threshold = NFFT_K(1000.0);
      }
    }
    else
    {
      /* TODO remove this */
      /* Set dummy values. */
      use_nfft = NO;
      use_fpt = NO;
      cutoff = 3;
      threshold = NFFT_K(1000.0);
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
        x_compare = (NFFT_R*) NFFT(malloc)(2*m_compare*sizeof(NFFT_R));
        d = 0;
        while (d < m_compare)
        {
          x1 = NFFT_K(2.0)*(((NFFT_R)rand())/RAND_MAX) - NFFT_K(1.0);
          x2 = NFFT_K(2.0)*(((NFFT_R)rand())/RAND_MAX) - NFFT_K(1.0);
          x3 = NFFT_K(2.0)*(((NFFT_R)rand())/RAND_MAX) - NFFT_K(1.0);
          temp = NFFT_SQRT(x1*x1+x2*x2+x3*x3);
          if (temp <= 1)
          {
            x_compare[2*d+1] = NFFT_M(acos)(x3);
            if (x_compare[2*d+1] == 0 || x_compare[2*d+1] == NFFT_KPI)
            {
              x_compare[2*d] = NFFT_K(0.0);
            }
            else
            {
              x_compare[2*d] = NFFT_M(atan2)(x2/NFFT_SIN(x_compare[2*d+1]),x1/NFFT_SIN(x_compare[2*d+1]));
            }
            x_compare[2*d] *= NFFT_K(1.0)/(NFFT_K(2.0)*NFFT_KPI);
            x_compare[2*d+1] *= NFFT_K(1.0)/(NFFT_K(2.0)*NFFT_KPI);
            d++;
          }
        }
        f_compare = (NFFT_C*) NFFT(malloc)(m_compare*sizeof(NFFT_C));
        f = (NFFT_C*) NFFT(malloc)(m_compare*sizeof(NFFT_C));
      }
    }

    /* Initialize maximum cut-off degree and grid size parameter. */
    NQ_max = 0;
    SQ_max = 0;

    /* Read the number of cut-off degrees. */
    fscanf(stdin,"bandwidths=%d\n",&iNQ_max);
    fprintf(stdout,"%d\n",iNQ_max);

    /* Allocate memory for the cut-off degrees and grid size parameters. */
    NQ = (int*) NFFT(malloc)(iNQ_max*sizeof(int));
    SQ = (int*) NFFT(malloc)(iNQ_max*sizeof(int));
    if (testmode == TIMING)
    {
      RQ = (int*) NFFT(malloc)(iNQ_max*sizeof(int));
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
    NFSFT(precompute)(NQ_max, threshold,
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U)), 0U);

    if (testmode == TIMING)
    {
      /* Allocate data structures. */
      f_hat = (NFFT_C*) NFFT(malloc)(NFSFT_F_HAT_SIZE(NQ_max)*sizeof(NFFT_C));
      f = (NFFT_C*) NFFT(malloc)(SQ_max*sizeof(NFFT_C));
      x_grid = (NFFT_R*) NFFT(malloc)(2*SQ_max*sizeof(NFFT_R));
      for (d = 0; d < SQ_max; d++)
      {
        f[d] = (((NFFT_R)rand())/RAND_MAX)-NFFT_K(0.5) + _Complex_I*((((NFFT_R)rand())/RAND_MAX)-NFFT_K(0.5));
        x_grid[2*d] = (((NFFT_R)rand())/RAND_MAX) - NFFT_K(0.5);
        x_grid[2*d+1] = (((NFFT_R)rand())/RAND_MAX) * NFFT_K(0.5);
      }
    }

    //fprintf(stderr,"Entering loop\n");
    //fflush(stderr);
    /* Process all cut-off bandwidths. */
    for (iNQ = 0; iNQ < iNQ_max; iNQ++)
    {
      if (testmode == TIMING)
      {
        NFSFT(init_guru)(&plan,NQ[iNQ],SQ[iNQ], NFSFT_NORMALIZED |
          ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
          ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
          PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFTW_MEASURE,
          cutoff);

        plan.f_hat = f_hat;
        plan.x = x_grid;
        plan.f = f;

        NFSFT(precompute_x)(&plan);

        t_avg = NFFT_K(0.0);

        for (i = 0; i < RQ[iNQ]; i++)
        {
          t0 = NFFT(clock_gettime_seconds)();

          if (use_nfsft != NO)
          {
            /* Execute the adjoint NFSFT transformation. */
            NFSFT(adjoint)(&plan);
          }
          else
          {
            /* Execute the adjoint direct NDSFT transformation. */
            NFSFT(adjoint_direct)(&plan);
          }

          t1 = NFFT(clock_gettime_seconds)();
          t_avg += t1 - t0;
        }

        t_avg = t_avg/((NFFT_R)RQ[iNQ]);

        NFSFT(finalize)(&plan);

        fprintf(stdout,"%+" NFFT__FES__ "\n", t_avg);
        fprintf(stderr,"%d: %4d %4d %+" NFFT__FES__ "\n", tc, NQ[iNQ], SQ[iNQ], t_avg);
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
              m_theta += (int)NFFT_FLOOR((2*NFFT_KPI)/NFFT_M(acos)((NFFT_COS(NFFT_KPI/(NFFT_R)SQ[iNQ])-
                NFFT_COS(k*NFFT_KPI/(NFFT_R)SQ[iNQ])*NFFT_COS(k*NFFT_KPI/(NFFT_R)SQ[iNQ]))/
                (NFFT_SIN(k*NFFT_KPI/(NFFT_R)SQ[iNQ])*NFFT_SIN(k*NFFT_KPI/(NFFT_R)SQ[iNQ]))));
              //fprintf(stderr,"ed: m_theta = %d\n",m_theta);
            }
            //fprintf(stderr,"ed: m_theta final = %d\n",m_theta);
            m_phi = 1;
            m_total = m_theta * m_phi;
            break;
        }

        /* Allocate memory for data structures. */
        w = (NFFT_R*) NFFT(malloc)(m_theta*sizeof(NFFT_R));
        x_grid = (NFFT_R*) NFFT(malloc)(2*m_total*sizeof(NFFT_R));

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
              fscanf(stdin,NFFT__FR__ "\n",&w[k]);
              w[k] *= (NFFT_K(2.0)*NFFT_KPI)/((NFFT_R)m_phi);
            }

            //fprintf(stderr,"Allocating theta and phi\n");
            //fflush(stderr);
            /* Allocate memory to store the grid's angles. */
            theta = (NFFT_R*) NFFT(malloc)(m_theta*sizeof(NFFT_R));
            phi = (NFFT_R*) NFFT(malloc)(m_phi*sizeof(NFFT_R));

            //if (theta == NULL || phi == NULL)
            //{
              //fprintf(stderr,"Couldn't allocate theta and phi\n");
              //fflush(stderr);
            //}


            /* Read angles theta. */
            for (k = 0; k < m_theta; k++)
            {
              fscanf(stdin,NFFT__FR__ "\n",&theta[k]);
            }

            /* Generate the grid angles phi. */
            for (n = 0; n < m_phi; n++)
            {
              phi[n] = n/((NFFT_R)m_phi);
              phi[n] -= ((phi[n]>=NFFT_K(0.5))?(NFFT_K(1.0)):(NFFT_K(0.0)));
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
            NFFT(free)(theta);
            NFFT(free)(phi);

            break;

          case GRID_CLENSHAW_CURTIS:

            /* Allocate memory to store the grid's angles. */
            theta = (NFFT_R*) NFFT(malloc)(m_theta*sizeof(NFFT_R));
            phi = (NFFT_R*) NFFT(malloc)(m_phi*sizeof(NFFT_R));

            /* Generate the grid angles theta. */
            for (k = 0; k < m_theta; k++)
            {
              theta[k] = k/((NFFT_R)2*(m_theta-1));
            }

            /* Generate the grid angles phi. */
            for (n = 0; n < m_phi; n++)
            {
              phi[n] = n/((NFFT_R)m_phi);
              phi[n] -= ((phi[n]>=NFFT_K(0.5))?(NFFT_K(1.0)):(NFFT_K(0.0)));
            }

            /* Generate quadrature weights. */
            fplan = FFTW(plan_r2r_1d)(SQ[iNQ]+1, w, w, FFTW_REDFT00, 0U);
            for (k = 0; k < SQ[iNQ]+1; k++)
            {
              w[k] = -NFFT_K(2.0)/(4*k*k-1);
            }
            FFTW(execute)(fplan);
            w[0] *= NFFT_K(0.5);

            for (k = 0; k < SQ[iNQ]+1; k++)
            {
              w[k] *= (NFFT_K(2.0)*NFFT_KPI)/((NFFT_R)(m_theta-1)*m_phi);
              w[m_theta-1-k] = w[k];
            }
            FFTW(destroy_plan)(fplan);

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
            NFFT(free)(theta);
            NFFT(free)(phi);

            break;

          case GRID_HEALPIX:

            d = 0;
            for (k = 1; k <= SQ[iNQ]-1; k++)
            {
              for (n = 0; n <= 4*k-1; n++)
              {
                x_grid[2*d+1] = 1 - (k*k)/((NFFT_R)(NFFT_K(3.0)*SQ[iNQ]*SQ[iNQ]));
                x_grid[2*d] =  ((n+NFFT_K(0.5))/(4*k));
                x_grid[2*d] -= (x_grid[2*d]>=NFFT_K(0.5))?(NFFT_K(1.0)):(NFFT_K(0.0));
                d++;
              }
            }

            d2 = d-1;

            for (k = SQ[iNQ]; k <= 3*SQ[iNQ]; k++)
            {
              for (n = 0; n <= 4*SQ[iNQ]-1; n++)
              {
                x_grid[2*d+1] = NFFT_K(2.0)/(3*SQ[iNQ])*(2*SQ[iNQ]-k);
                x_grid[2*d] = (n+((k%2==0)?(NFFT_K(0.5)):(NFFT_K(0.0))))/(4*SQ[iNQ]);
                x_grid[2*d] -= (x_grid[2*d]>=NFFT_K(0.5))?(NFFT_K(1.0)):(NFFT_K(0.0));
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
              x_grid[2*d+1] = NFFT_M(acos)(x_grid[2*d+1])/(NFFT_K(2.0)*NFFT_KPI);
            }

            w[0] = (NFFT_K(4.0)*NFFT_KPI)/(m_total);
            break;

          case GRID_EQUIDISTRIBUTION:
          case GRID_EQUIDISTRIBUTION_UNIFORM:
            /* TODO Compute the weights. */

            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              w_temp = (NFFT_R*) NFFT(malloc)((SQ[iNQ]+1)*sizeof(NFFT_R));
              fplan = FFTW(plan_r2r_1d)(SQ[iNQ]/2+1, w_temp, w_temp, FFTW_REDFT00, 0U);
              for (k = 0; k < SQ[iNQ]/2+1; k++)
              {
                w_temp[k] = -NFFT_K(2.0)/(4*k*k-1);
              }
              FFTW(execute)(fplan);
              w_temp[0] *= NFFT_K(0.5);

              for (k = 0; k < SQ[iNQ]/2+1; k++)
              {
                w_temp[k] *= (NFFT_K(2.0)*NFFT_KPI)/((NFFT_R)(SQ[iNQ]));
                w_temp[SQ[iNQ]-k] = w_temp[k];
              }
              FFTW(destroy_plan)(fplan);
            }

            d = 0;
            x_grid[2*d] = -NFFT_K(0.5);
            x_grid[2*d+1] = NFFT_K(0.0);
            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              w[d] = w_temp[0];
            }
            else
            {
              w[d] = (NFFT_K(4.0)*NFFT_KPI)/(m_total);
            }
            d = 1;
            x_grid[2*d] = -NFFT_K(0.5);
            x_grid[2*d+1] = NFFT_K(0.5);
            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              w[d] = w_temp[SQ[iNQ]];
            }
            else
            {
              w[d] = (NFFT_K(4.0)*NFFT_KPI)/(m_total);
            }
            d = 2;

            for (k = 1; k < SQ[iNQ]; k++)
            {
              theta_s = (NFFT_R)k*NFFT_KPI/(NFFT_R)SQ[iNQ];
              M = (int)NFFT_FLOOR((NFFT_K(2.0)*NFFT_KPI)/NFFT_M(acos)((NFFT_COS(NFFT_KPI/(NFFT_R)SQ[iNQ])-
                NFFT_COS(theta_s)*NFFT_COS(theta_s))/(NFFT_SIN(theta_s)*NFFT_SIN(theta_s))));

              for (n = 0; n < M; n++)
              {
                x_grid[2*d] = (n + NFFT_K(0.5))/M;
                x_grid[2*d] -= (x_grid[2*d]>=NFFT_K(0.5))?(NFFT_K(1.0)):(NFFT_K(0.0));
                x_grid[2*d+1] = theta_s/(NFFT_K(2.0)*NFFT_KPI);
                if (gridtype == GRID_EQUIDISTRIBUTION)
                {
                  w[d] = w_temp[k]/((NFFT_R)(M));
                }
                else
                {
                  w[d] = (NFFT_K(4.0)*NFFT_KPI)/(m_total);
                }
                d++;
              }
            }

            if (gridtype == GRID_EQUIDISTRIBUTION)
            {
              NFFT(free)(w_temp);
            }
            break;

          default:
            break;
        }

        /* Allocate memory for grid values. */
        f_grid = (NFFT_C*) NFFT(malloc)(m_total*sizeof(NFFT_C));

        if (mode == RANDOM)
        {
        }
        else
        {
          m_compare = m_total;
          f_compare = (NFFT_C*) NFFT(malloc)(m_compare*sizeof(NFFT_C));
          x_compare = x_grid;
          f = f_grid;
        }

        //fprintf(stderr,"Generating test function\n");
        //fflush(stderr);
        switch (testfunction)
        {
          case FUNCTION_RANDOM_BANDLIMITED:
            f_hat_gen = (NFFT_C*) NFFT(malloc)(NFSFT_F_HAT_SIZE(N)*sizeof(NFFT_C));
            //fprintf(stderr,"Generating random test function\n");
            //fflush(stderr);
            /* Generate random function samples by sampling a bandlimited
             * function. */
            NFSFT(init_guru)(&plan_gen,N,m_total, NFSFT_NORMALIZED |
              ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
              ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
              ((N>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT,
              cutoff);

            plan_gen.f_hat = f_hat_gen;
            plan_gen.x = x_grid;
            plan_gen.f = f_grid;

            NFSFT(precompute_x)(&plan_gen);

            for (k = 0; k < plan_gen.N_total; k++)
            {
              f_hat_gen[k] = NFFT_K(0.0);
            }

            for (k = 0; k <= N; k++)
            {
              for (n = -k; n <= k; n++)
              {
                f_hat_gen[NFSFT_INDEX(k,n,&plan_gen)] =
                (((NFFT_R)rand())/RAND_MAX)-NFFT_K(0.5) + _Complex_I*((((NFFT_R)rand())/RAND_MAX)-NFFT_K(0.5));
              }
            }

            if (use_nfsft != NO)
            {
              /* Execute the NFSFT transformation. */
              NFSFT(trafo)(&plan_gen);
            }
            else
            {
              /* Execute the direct NDSFT transformation. */
              NFSFT(trafo_direct)(&plan_gen);
            }

            NFSFT(finalize)(&plan_gen);

            if (mode == RANDOM)
            {
              NFSFT(init_guru)(&plan_gen,N,m_compare, NFSFT_NORMALIZED |
                ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
                ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
                ((N>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT,
                cutoff);

              plan_gen.f_hat = f_hat_gen;
              plan_gen.x = x_compare;
              plan_gen.f = f_compare;

              NFSFT(precompute_x)(&plan_gen);

              if (use_nfsft != NO)
              {
                /* Execute the NFSFT transformation. */
                NFSFT(trafo)(&plan_gen);
              }
              else
              {
                /* Execute the direct NDSFT transformation. */
                NFSFT(trafo_direct)(&plan_gen);
              }

              NFSFT(finalize)(&plan_gen);
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }

            NFFT(free)(f_hat_gen);

            break;

          case FUNCTION_F1:
            for (d = 0; d < m_total; d++)
            {
              x1 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x2 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x3 = NFFT_COS(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
              f_grid[d] = x1*x2*x3;
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x2 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x3 = NFFT_COS(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
                f_compare[d] = x1*x2*x3;
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }
            break;
          case FUNCTION_F2:
            for (d = 0; d < m_total; d++)
            {
              x1 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x2 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x3 = NFFT_COS(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
              f_grid[d] = NFFT_K(0.1)*NFFT_EXP(x1+x2+x3);
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x2 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x3 = NFFT_COS(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
                f_compare[d] = NFFT_K(0.1)*NFFT_EXP(x1+x2+x3);
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }
            break;
          case FUNCTION_F3:
            for (d = 0; d < m_total; d++)
            {
              x1 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x2 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x3 = NFFT_COS(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
              temp = NFFT_SQRT(x1*x1)+NFFT_SQRT(x2*x2)+NFFT_SQRT(x3*x3);
              f_grid[d] = NFFT_K(0.1)*temp;
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x2 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x3 = NFFT_COS(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
                temp = NFFT_SQRT(x1*x1)+NFFT_SQRT(x2*x2)+NFFT_SQRT(x3*x3);
                f_compare[d] = NFFT_K(0.1)*temp;
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }
            break;
          case FUNCTION_F4:
            for (d = 0; d < m_total; d++)
            {
              x1 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x2 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x3 = NFFT_COS(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
              temp = NFFT_SQRT(x1*x1)+NFFT_SQRT(x2*x2)+NFFT_SQRT(x3*x3);
              f_grid[d] = NFFT_K(1.0)/(temp);
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x2 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x3 = NFFT_COS(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
                temp = NFFT_SQRT(x1*x1)+NFFT_SQRT(x2*x2)+NFFT_SQRT(x3*x3);
                f_compare[d] = NFFT_K(1.0)/(temp);
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }
            break;
          case FUNCTION_F5:
            for (d = 0; d < m_total; d++)
            {
              x1 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x2 = NFFT_SIN(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_grid[2*d]*NFFT_K(2.0)*NFFT_KPI);
              x3 = NFFT_COS(x_grid[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
              temp = NFFT_SQRT(x1*x1)+NFFT_SQRT(x2*x2)+NFFT_SQRT(x3*x3);
              f_grid[d] = NFFT_K(0.1)*NFFT_SIN(1+temp)*NFFT_SIN(1+temp);
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                x1 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_COS(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x2 = NFFT_SIN(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI)*NFFT_SIN(x_compare[2*d]*NFFT_K(2.0)*NFFT_KPI);
                x3 = NFFT_COS(x_compare[2*d+1]*NFFT_K(2.0)*NFFT_KPI);
                temp = NFFT_SQRT(x1*x1)+NFFT_SQRT(x2*x2)+NFFT_SQRT(x3*x3);
                f_compare[d] = NFFT_K(0.1)*NFFT_SIN(1+temp)*NFFT_SIN(1+temp);
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }
            break;
          case FUNCTION_F6:
            for (d = 0; d < m_total; d++)
            {
              if (x_grid[2*d+1] <= NFFT_K(0.25))
              {
                f_grid[d] = NFFT_K(1.0);
              }
              else
              {
                f_grid[d] = NFFT_K(1.0)/(NFFT_SQRT(1+3*NFFT_COS(NFFT_K(2.0)*NFFT_KPI*x_grid[2*d+1])*NFFT_COS(NFFT_K(2.0)*NFFT_KPI*x_grid[2*d+1])));
              }
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                if (x_compare[2*d+1] <= NFFT_K(0.25))
                {
                  f_compare[d] = NFFT_K(1.0);
                }
                else
                {
                  f_compare[d] = NFFT_K(1.0)/(NFFT_SQRT(1+3*NFFT_COS(NFFT_K(2.0)*NFFT_KPI*x_compare[2*d+1])*NFFT_COS(NFFT_K(2.0)*NFFT_KPI*x_compare[2*d+1])));
                }
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }
            break;
          default:
            //fprintf(stderr,"Generating one function\n");
            //fflush(stderr);
            for (d = 0; d < m_total; d++)
            {
              f_grid[d] = NFFT_K(1.0);
            }
            if (mode == RANDOM)
            {
              for (d = 0; d < m_compare; d++)
              {
                f_compare[d] = NFFT_K(1.0);
              }
            }
            else
            {
              memcpy(f_compare,f_grid,m_total*sizeof(NFFT_C));
            }
            break;
        }

        //fprintf(stderr,"Initializing trafo\n");
        //fflush(stderr);
        /* Init transform plan. */
        NFSFT(init_guru)(&plan_adjoint,NQ[iNQ],m_total, NFSFT_NORMALIZED |
          ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
          ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
          ((NQ[iNQ]>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT,
          cutoff);

        plan_adjoint_ptr = &plan_adjoint;

        if (mode == RANDOM)
        {
          NFSFT(init_guru)(&plan,NQ[iNQ],m_compare, NFSFT_NORMALIZED |
            ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
            ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)),
            ((NQ[iNQ]>512)?(0U):(PRE_PHI_HUT | PRE_PSI)) | FFTW_INIT,
            cutoff);
          plan_ptr = &plan;
        }
        else
        {
          plan_ptr = &plan_adjoint;
        }

        f_hat = (NFFT_C*) NFFT(malloc)(NFSFT_F_HAT_SIZE(NQ[iNQ])*sizeof(NFFT_C));

        plan_adjoint_ptr->f_hat = f_hat;
        plan_adjoint_ptr->x = x_grid;
        plan_adjoint_ptr->f = f_grid;

        plan_ptr->f_hat = f_hat;
        plan_ptr->x = x_compare;
        plan_ptr->f = f;

        //fprintf(stderr,"Precomputing for x\n");
        //fflush(stderr);
        NFSFT(precompute_x)(plan_adjoint_ptr);
        if (plan_adjoint_ptr != plan_ptr)
        {
          NFSFT(precompute_x)(plan_ptr);
        }

        /* Initialize cumulative time variable. */
        t_avg = NFFT_K(0.0);
        err_infty_avg = NFFT_K(0.0);
        err_2_avg = NFFT_K(0.0);

        /* Cycle through all runs. */
        for (i = 0; i < 1/*repetitions*/; i++)
        {
          //fprintf(stderr,"Copying original values\n");
          //fflush(stderr);
          /* Copy exact funtion values to working array. */
          //memcpy(f,f_grid,m_total*sizeof(NFFT_C));

          /* Initialize time measurement. */
          t0 = NFFT(clock_gettime_seconds)();

          //fprintf(stderr,"Multiplying with quadrature weights\n");
          //fflush(stderr);
          /* Multiplication with the quadrature weights. */
          /*fprintf(stderr,"\n");*/
          d = 0;
          for (k = 0; k < m_theta; k++)
          {
            for (n = 0; n < m_phi; n++)
            {
              /*fprintf(stderr,"f_ref[%d] = %" NFFT__FES__ " + I*%" NFFT__FES__ ",\t f[%d] = %" NFFT__FES__ " + I*%" NFFT__FES__ ",  \t w[%d] = %" NFFT__FES__ "\n",
              d,NFFT_CREAL(f_ref[d]),NFFT_CIMAG(f_ref[d]),d,NFFT_CREAL(f[d]),NFFT_CIMAG(f[d]),k,
              w[k]);*/
              f_grid[d] *= w[k];
              d++;
            }
          }

          t1 = NFFT(clock_gettime_seconds)();
          t_avg += t1 - t0;

          NFFT(free)(w);

          t0 = NFFT(clock_gettime_seconds)();

          /*fprintf(stderr,"\n");
          d = 0;
          for (d = 0; d < grid_total; d++)
          {
            fprintf(stderr,"f[%d] = %" NFFT__FES__ " + I*%" NFFT__FES__ ", theta[%d] = %" NFFT__FES__ ", phi[%d] = %" NFFT__FES__ "\n",
                    d,NFFT_CREAL(f[d]),NFFT_CIMAG(f[d]),d,x[2*d+1],d,x[2*d]);
          }*/

          //fprintf(stderr,"Executing adjoint\n");
          //fflush(stderr);
          /* Check if the fast NFSFT algorithm shall be tested. */
          if (use_nfsft != NO)
          {
            /* Execute the adjoint NFSFT transformation. */
            NFSFT(adjoint)(plan_adjoint_ptr);
          }
          else
          {
            /* Execute the adjoint direct NDSFT transformation. */
            NFSFT(adjoint_direct)(plan_adjoint_ptr);
          }

          /* Multiplication with the Fourier-Legendre coefficients. */
          /*for (k = 0; k <= m[im]; k++)
          {
            for (n = -k; n <= k; n++)
            {
              fprintf(stderr,"f_hat[%d,%d] = %" NFFT__FES__ "\t + I*%" NFFT__FES__ "\n",k,n,
                      NFFT_CREAL(f_hat[NFSFT_INDEX(k,n,&plan_adjoint)]),
                      NFFT_CIMAG(f_hat[NFSFT_INDEX(k,n,&plan_adjoint)]));
            }
          }*/

          //fprintf(stderr,"Executing trafo\n");
          //fflush(stderr);
          if (use_nfsft != NO)
          {
            /* Execute the NFSFT transformation. */
            NFSFT(trafo)(plan_ptr);
          }
          else
          {
            /* Execute the direct NDSFT transformation. */
            NFSFT(trafo_direct)(plan_ptr);
          }

          t1 = NFFT(clock_gettime_seconds)();
          t_avg += t1 - t0;

          //fprintf(stderr,"Finalizing\n");
          //fflush(stderr);
          /* Finalize the NFSFT plans */
          NFSFT(finalize)(plan_adjoint_ptr);
          if (plan_ptr != plan_adjoint_ptr)
          {
            NFSFT(finalize)(plan_ptr);
          }

          /* Free data arrays. */
          NFFT(free)(f_hat);
          NFFT(free)(x_grid);

          err_infty_avg += NFFT(error_l_infty_complex)(f, f_compare, m_compare);
          err_2_avg += NFFT(error_l_2_complex)(f, f_compare, m_compare);

          NFFT(free)(f_grid);

          if (mode == RANDOM)
          {
          }
          else
          {
            NFFT(free)(f_compare);
          }

          /*for (d = 0; d < m_total; d++)
          {
            fprintf(stderr,"f_ref[%d] = %" NFFT__FES__ " + I*%" NFFT__FES__ ",\t f[%d] = %" NFFT__FES__ " + I*%" NFFT__FES__ "\n",
              d,NFFT_CREAL(f_ref[d]),NFFT_CIMAG(f_ref[d]),d,NFFT_CREAL(f[d]),NFFT_CIMAG(f[d]));
          }*/
        }

        //fprintf(stderr,"Calculating the error\n");
        //fflush(stderr);
        /* Calculate average time needed. */
        t_avg = t_avg/((NFFT_R)repetitions);

        /* Calculate the average error. */
        err_infty_avg = err_infty_avg/((NFFT_R)repetitions);

        /* Calculate the average error. */
        err_2_avg = err_2_avg/((NFFT_R)repetitions);

        /* Print out the error measurements. */
        fprintf(stdout,"%+" NFFT__FES__ " %+" NFFT__FES__ " %+" NFFT__FES__ "\n", t_avg, err_infty_avg, err_2_avg);
        fprintf(stderr,"%d: %4d %4d %+" NFFT__FES__ " %+" NFFT__FES__ " %+" NFFT__FES__ "\n", tc, NQ[iNQ], SQ[iNQ],
          t_avg, err_infty_avg, err_2_avg);
      }
    } /* for (im = 0; im < im_max; im++) - Process all cut-off
       * bandwidths.*/
    fprintf(stderr,"\n");

    /* Delete precomputed data. */
    NFSFT(forget)();

    /* Free memory for cut-off bandwidths and grid size parameters. */
    NFFT(free)(NQ);
    NFFT(free)(SQ);
    if (testmode == TIMING)
    {
      NFFT(free)(RQ);
    }

    if (mode == RANDOM)
    {
      NFFT(free)(x_compare);
      NFFT(free)(f_compare);
      NFFT(free)(f);
    }

    if (testmode == TIMING)
    {
      /* Allocate data structures. */
      NFFT(free)(f_hat);
      NFFT(free)(f);
      NFFT(free)(x_grid);
    }

  } /* for (tc = 0; tc < tc_max; tc++) - Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}
/* \} */
