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

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "nfft3.h"

/** Swap two vectors. */
#define CSWAP(x,y) {double _Complex * NFFT_SWAP_temp__; \
  NFFT_SWAP_temp__=(x); (x)=(y); (y)=NFFT_SWAP_temp__;}

void simple_test_nnfft_1d(void)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nnfft_plan my_plan;                    /**< plan for the nfft                */

  int N[1];
  N[0]=10;

  /** init an one dimensional plan */
  nnfft_init(&my_plan, 1, 3, 19, N);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[j]=((double)rand())/((double)RAND_MAX)-0.5;
  }
  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[j]=((double)rand())/((double)RAND_MAX)-0.5;
  }

  /** precompute psi, the entries of the matrix B */
/*  if(my_plan.nnfft_flags & PRE_PSI)
        nnfft_precompute_psi(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
 	nnfft_precompute_full_psi(&my_plan);
  if(my_plan.nnfft_flags & PRE_LIN_PSI)
	 nnfft_precompute_lin_psi(&my_plan);
  /** precompute phi_hut, the entries of the matrix D */
/*  if(my_plan.nnfft_flags & PRE_PHI_HUT)
	  nnfft_precompute_phi_hut(&my_plan);
*/

nnfft_precompute_one_psi(&my_plan);


  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/((double)RAND_MAX) + _Complex_I*((double)rand())/((double)RAND_MAX);

  nfft_vpr_complex(my_plan.f_hat,my_plan.N_total,"given Fourier coefficients, vector f_hat");

  /** direct trafo and show the result */
  nnfft_trafo_direct(&my_plan);
  nfft_vpr_complex(my_plan.f,my_plan.M_total,"nndft, vector f");

  /** approx. trafo and show the result */
  nnfft_trafo(&my_plan);
  nfft_vpr_complex(my_plan.f,my_plan.M_total,"nnfft, vector f");

  /** finalise the one dimensional plan */
  nnfft_finalize(&my_plan);
}

static void simple_test_adjoint_nnfft_1d(void)
{
  int j;                                 /**< index for nodes and freqencies   */
  nnfft_plan my_plan;                    /**< plan for the nfft                */

  int N[1];
  N[0]=12;

  /** init an one dimensional plan */
  nnfft_init(&my_plan, 1, 20, 33, N);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[j]=((double)rand())/((double)RAND_MAX)-0.5;
  }
  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[j]=((double)rand())/((double)RAND_MAX)-0.5;
  }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    nnfft_precompute_psi(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
    nnfft_precompute_full_psi(&my_plan);

  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    nnfft_precompute_lin_psi(&my_plan);

  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(j=0;j<my_plan.M_total;j++)
    my_plan.f[j] = ((double)rand())/((double)RAND_MAX) + _Complex_I*((double)rand())/((double)RAND_MAX);

  nfft_vpr_complex(my_plan.f,my_plan.M_total,"given Samples, vector f");

  /** direct trafo and show the result */
  nnfft_adjoint_direct(&my_plan);
  nfft_vpr_complex(my_plan.f_hat,my_plan.N_total,"adjoint nndft, vector f_hat");

  /** approx. trafo and show the result */
  nnfft_adjoint(&my_plan);
  nfft_vpr_complex(my_plan.f_hat,my_plan.N_total,"adjoint nnfft, vector f_hat");

  /** finalise the one dimensional plan */
  nnfft_finalize(&my_plan);
}

static void simple_test_nnfft_2d(void)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nnfft_plan my_plan;                    /**< plan for the nfft                */

  int N[2];
  N[0]=12;
  N[1]=14;

  /** init an one dimensional plan */
  nnfft_init(&my_plan, 2,12*14,19, N);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[2*j]=((double)rand())/((double)RAND_MAX)-0.5;
    my_plan.x[2*j+1]=((double)rand())/((double)RAND_MAX)-0.5;
  }

  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[2*j]=((double)rand())/((double)RAND_MAX)-0.5;
    my_plan.v[2*j+1]=((double)rand())/((double)RAND_MAX)-0.5;
  }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    nnfft_precompute_psi(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
    nnfft_precompute_full_psi(&my_plan);

  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    nnfft_precompute_lin_psi(&my_plan);

  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((double)rand())/((double)RAND_MAX) + _Complex_I*((double)rand())/((double)RAND_MAX);

  nfft_vpr_complex(my_plan.f_hat,12,
        "given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  nnfft_trafo_direct(&my_plan);
  nfft_vpr_complex(my_plan.f,my_plan.M_total,"ndft, vector f");

  /** approx. trafo and show the result */
  nnfft_trafo(&my_plan);
  nfft_vpr_complex(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** finalise the one dimensional plan */
  nnfft_finalize(&my_plan);
}

static void simple_test_innfft_1d(void)
{
  int j,k,l,N=8;                        /**< index for nodes, freqencies, iter*/
  nnfft_plan my_plan;                   /**< plan for the nnfft               */
  solver_plan_complex my_iplan;         /**< plan for the inverse nnfft       */

  /** initialise an one dimensional plan */
  nnfft_init(&my_plan,1 ,8 ,8 ,&N);

  /** initialise my_iplan */
  solver_init_advanced_complex(&my_iplan,(nfft_mv_plan_complex*)(&my_plan),CGNR);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    my_plan.x[j]=((double)rand())/((double)RAND_MAX)-0.5;

  /** init pseudo random nodes */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.v[k]=((double)rand())/((double)RAND_MAX)-0.5;

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    nnfft_precompute_psi(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
      nnfft_precompute_full_psi(&my_plan);

  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);

  /** init pseudo random samples (real) and show them */
  for(j=0;j<my_plan.M_total;j++)
    my_iplan.y[j] = ((double)rand())/((double)RAND_MAX);

  nfft_vpr_complex(my_iplan.y,my_plan.M_total,"given data, vector given_f");

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k] = 0.0;

  nfft_vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
        "approximate solution, vector f_hat_iter");

  /** solve the system */
  solver_before_loop_complex(&my_iplan);

  for(l=0;l<8;l++)
  {
    printf("iteration l=%d\n",l);
    solver_loop_one_step_complex(&my_iplan);
    nfft_vpr_complex(my_iplan.f_hat_iter,my_plan.N_total,
          "approximate solution, vector f_hat_iter");

    CSWAP(my_iplan.f_hat_iter,my_plan.f_hat);
    nnfft_trafo(&my_plan);
    nfft_vpr_complex(my_plan.f,my_plan.M_total,"fitting the data, vector f");
    CSWAP(my_iplan.f_hat_iter,my_plan.f_hat);
  }

  solver_finalize_complex(&my_iplan);
  nnfft_finalize(&my_plan);
}

static void measure_time_nnfft_1d(void)
{
  int j,k;                              /**< index for nodes and freqencies   */
  nnfft_plan my_plan;                    /**< plan for the nfft                */
  int my_N,N=4;
  double t;
  double t0, t1;

  for(my_N=16; my_N<=16384; my_N*=2)
  {
    nnfft_init(&my_plan,1,my_N,my_N,&N);

    for(j=0;j<my_plan.M_total;j++)
      my_plan.x[j]=((double)rand())/((double)RAND_MAX)-0.5;

    for(j=0;j<my_plan.N_total;j++)
      my_plan.v[j]=((double)rand())/((double)RAND_MAX)-0.5;

    if(my_plan.nnfft_flags & PRE_PSI)
      nnfft_precompute_psi(&my_plan);

    if(my_plan.nnfft_flags & PRE_FULL_PSI)
        nnfft_precompute_full_psi(&my_plan);

    if(my_plan.nnfft_flags & PRE_PHI_HUT)
      nnfft_precompute_phi_hut(&my_plan);

    for(k=0;k<my_plan.N_total;k++)
      my_plan.f_hat[k] = ((double)rand())/((double)RAND_MAX) + _Complex_I*((double)rand())/((double)RAND_MAX);

    t0 = nfft_clock_gettime_seconds();
    nnfft_trafo_direct(&my_plan);
    t1 = nfft_clock_gettime_seconds();
    t = t1 - t0;
    printf("t_nndft=%e,\t",t);

    t0 = nfft_clock_gettime_seconds();
    nnfft_trafo(&my_plan);
    t1 = nfft_clock_gettime_seconds();
    t = t1 - t0;
    printf("t_nnfft=%e\t",t);

    printf("(N=M=%d)\n",my_N);

    nnfft_finalize(&my_plan);
  }
}

int main(void)
{
  system("clear");
  printf("1) computing a one dimensional nndft, nnfft\n\n");
  simple_test_nnfft_1d();

  /*getc(stdin);

  system("clear");
  printf("1a) computing a one dimensional adjoint nndft, nnfft\n\n");
  simple_test_adjoint_nnfft_1d();

  getc(stdin);

  system("clear");
  printf("2) computing a two dimensional nndft, nnfft\n\n");
  simple_test_nnfft_2d();

  getc(stdin);

  system("clear");
  printf("3) computing a one dimensional innfft\n\n");
  simple_test_innfft_1d();

  getc(stdin);

  system("clear");
  printf("4) computing times for one dimensional nnfft\n\n");
  measure_time_nnfft_1d();
*/
  return 1;
}
