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
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

void simple_test_nnfft_1d(void)
{
  int j,k;                              /**< index for nodes and freqencies   */
  NNFFT(plan) my_plan;                    /**< plan for the nfft                */

  int N[1];
  N[0]=10;

  /** init an one dimensional plan */
  NNFFT(init)(&my_plan, 1, 3, 19, N);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
  }
  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
  }

  /** precompute psi, the entries of the matrix B */
/*  if(my_plan.nnfft_flags & PRE_PSI)
        NNFFT(precompute_psi)(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
 	NNFFT(precompute_full_psi)(&my_plan);
  if(my_plan.nnfft_flags & PRE_LIN_PSI)
	 NNFFT(precompute_lin_psi)(&my_plan);
  /** precompute phi_hut, the entries of the matrix D */
/*  if(my_plan.nnfft_flags & PRE_PHI_HUT)
	  NNFFT(precompute_phi_hut)(&my_plan);
*/

NNFFT(precompute_one_psi)(&my_plan);


  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((NFFT_R)rand())/((NFFT_R)RAND_MAX) + _Complex_I*((NFFT_R)rand())/((NFFT_R)RAND_MAX);

  NFFT(vpr_complex)(my_plan.f_hat,my_plan.N_total,"given Fourier coefficients, vector f_hat");

  /** direct trafo and show the result */
  NNFFT(trafo_direct)(&my_plan);
  NFFT(vpr_complex)(my_plan.f,my_plan.M_total,"nndft, vector f");

  /** approx. trafo and show the result */
  NNFFT(trafo)(&my_plan);
  NFFT(vpr_complex)(my_plan.f,my_plan.M_total,"nnfft, vector f");

  /** finalise the one dimensional plan */
  NNFFT(finalize)(&my_plan);
}

static void simple_test_adjoint_nnfft_1d(void)
{
  int j;                                 /**< index for nodes and freqencies   */
  NNFFT(plan) my_plan;                    /**< plan for the nfft                */

  int N[1];
  N[0]=12;

  /** init an one dimensional plan */
  NNFFT(init)(&my_plan, 1, 20, 33, N);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
  }
  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
  }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    NNFFT(precompute_psi)(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
    NNFFT(precompute_full_psi)(&my_plan);

  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    NNFFT(precompute_lin_psi)(&my_plan);

  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    NNFFT(precompute_phi_hut)(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(j=0;j<my_plan.M_total;j++)
    my_plan.f[j] = ((NFFT_R)rand())/((NFFT_R)RAND_MAX) + _Complex_I*((NFFT_R)rand())/((NFFT_R)RAND_MAX);

  NFFT(vpr_complex)(my_plan.f,my_plan.M_total,"given Samples, vector f");

  /** direct trafo and show the result */
  NNFFT(adjoint_direct)(&my_plan);
  NFFT(vpr_complex)(my_plan.f_hat,my_plan.N_total,"adjoint nndft, vector f_hat");

  /** approx. trafo and show the result */
  NNFFT(adjoint)(&my_plan);
  NFFT(vpr_complex)(my_plan.f_hat,my_plan.N_total,"adjoint nnfft, vector f_hat");

  /** finalise the one dimensional plan */
  NNFFT(finalize)(&my_plan);
}

static void simple_test_nnfft_2d(void)
{
  int j,k;                              /**< index for nodes and freqencies   */
  NNFFT(plan) my_plan;                    /**< plan for the nfft                */

  int N[2];
  N[0]=12;
  N[1]=14;

  /** init an one dimensional plan */
  NNFFT(init)(&my_plan, 2,12*14,19, N);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
  {
    my_plan.x[2*j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
    my_plan.x[2*j+1]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
  }

  /** init pseudo random nodes */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_plan.v[2*j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
    my_plan.v[2*j+1]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);
  }

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    NNFFT(precompute_psi)(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
    NNFFT(precompute_full_psi)(&my_plan);

  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    NNFFT(precompute_lin_psi)(&my_plan);

  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    NNFFT(precompute_phi_hut)(&my_plan);

  /** init pseudo random Fourier coefficients and show them */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.f_hat[k] = ((NFFT_R)rand())/((NFFT_R)RAND_MAX) + _Complex_I*((NFFT_R)rand())/((NFFT_R)RAND_MAX);

  NFFT(vpr_complex)(my_plan.f_hat,12,
        "given Fourier coefficients, vector f_hat (first 12 entries)");

  /** direct trafo and show the result */
  NNFFT(trafo_direct)(&my_plan);
  NFFT(vpr_complex)(my_plan.f,my_plan.M_total,"ndft, vector f");

  /** approx. trafo and show the result */
  NNFFT(trafo)(&my_plan);
  NFFT(vpr_complex)(my_plan.f,my_plan.M_total,"nfft, vector f");

  /** finalise the one dimensional plan */
  NNFFT(finalize)(&my_plan);
}

static void simple_test_innfft_1d(void)
{
  int j,k,l,N=8;                        /**< index for nodes, freqencies, iter*/
  NNFFT(plan) my_plan;                   /**< plan for the nnfft               */
  SOLVER(plan_complex) my_iplan;         /**< plan for the inverse nnfft       */

  /** initialise an one dimensional plan */
  NNFFT(init)(&my_plan,1 ,8 ,8 ,&N);

  /** initialise my_iplan */
  SOLVER(init_advanced_complex)(&my_iplan,(NFFT(mv_plan_complex)*)(&my_plan),CGNR);

  /** init pseudo random nodes */
  for(j=0;j<my_plan.M_total;j++)
    my_plan.x[j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);

  /** init pseudo random nodes */
  for(k=0;k<my_plan.N_total;k++)
    my_plan.v[k]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);

  /** precompute psi, the entries of the matrix B */
  if(my_plan.nnfft_flags & PRE_PSI)
    NNFFT(precompute_psi)(&my_plan);

  if(my_plan.nnfft_flags & PRE_FULL_PSI)
      NNFFT(precompute_full_psi)(&my_plan);

  /** precompute phi_hut, the entries of the matrix D */
  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    NNFFT(precompute_phi_hut)(&my_plan);

  /** init pseudo random samples (real) and show them */
  for(j=0;j<my_plan.M_total;j++)
    my_iplan.y[j] = ((NFFT_R)rand())/((NFFT_R)RAND_MAX);

  NFFT(vpr_complex)(my_iplan.y,my_plan.M_total,"given data, vector given_f");

  /** initialise some guess f_hat_0 */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k] = NFFT_K(0.0);

  NFFT(vpr_complex)(my_iplan.f_hat_iter,my_plan.N_total,
        "approximate solution, vector f_hat_iter");

  /** solve the system */
  SOLVER(before_loop_complex)(&my_iplan);

  for(l=0;l<8;l++)
  {
    printf("iteration l=%d\n",l);
    SOLVER(loop_one_step_complex)(&my_iplan);
    NFFT(vpr_complex)(my_iplan.f_hat_iter,my_plan.N_total,
          "approximate solution, vector f_hat_iter");

    NFFT_CSWAP(my_iplan.f_hat_iter,my_plan.f_hat);
    NNFFT(trafo)(&my_plan);
    NFFT(vpr_complex)(my_plan.f,my_plan.M_total,"fitting the data, vector f");
    NFFT_CSWAP(my_iplan.f_hat_iter,my_plan.f_hat);
  }

  SOLVER(finalize_complex)(&my_iplan);
  NNFFT(finalize)(&my_plan);
}

static void measure_time_nnfft_1d(void)
{
  int j,k;                              /**< index for nodes and freqencies   */
  NNFFT(plan) my_plan;                    /**< plan for the nfft                */
  int my_N,N=4;
  double t;
  double t0, t1;

  for(my_N=16; my_N<=16384; my_N*=2)
  {
    NNFFT(init)(&my_plan,1,my_N,my_N,&N);

    for(j=0;j<my_plan.M_total;j++)
      my_plan.x[j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);

    for(j=0;j<my_plan.N_total;j++)
      my_plan.v[j]=((NFFT_R)rand())/((NFFT_R)RAND_MAX)-NFFT_K(0.5);

    if(my_plan.nnfft_flags & PRE_PSI)
      NNFFT(precompute_psi)(&my_plan);

    if(my_plan.nnfft_flags & PRE_FULL_PSI)
        NNFFT(precompute_full_psi)(&my_plan);

    if(my_plan.nnfft_flags & PRE_PHI_HUT)
      NNFFT(precompute_phi_hut)(&my_plan);

    for(k=0;k<my_plan.N_total;k++)
      my_plan.f_hat[k] = ((NFFT_R)rand())/((NFFT_R)RAND_MAX) + _Complex_I*((NFFT_R)rand())/((NFFT_R)RAND_MAX);

    t0 = NFFT(clock_gettime_seconds)();
    NNFFT(trafo_direct)(&my_plan);
    t1 = NFFT(clock_gettime_seconds)();
    t = t1 - t0;
    printf("t_nndft=%e,\t",t);

    t0 = NFFT(clock_gettime_seconds)();
    NNFFT(trafo)(&my_plan);
    t1 = NFFT(clock_gettime_seconds)();
    t = t1 - t0;
    printf("t_nnfft=%e\t",t);

    printf("(N=M=%d)\n",my_N);

    NNFFT(finalize)(&my_plan);
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
