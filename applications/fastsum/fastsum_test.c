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

/*! \file fastsum_test.c
 *  \brief Simple test program for the fast NFFT-based summation algorithm.
 *
 *  \author Markus Fenn
 *  \date 2006
 */
#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
  #include <complex.h>
#endif

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "fastsum.h"
#include "kernels.h"
#include "infft.h"

/**
 * \defgroup applications_fastsum_test fastsum_test
 * \ingroup applications_fastsum
 * \{
 */

int main(int argc, char **argv)
{
  int j, k; /**< indices                 */
  int d; /**< number of dimensions    */
  int N; /**< number of source nodes  */
  int M; /**< number of target nodes  */
  int n; /**< expansion degree        */
  int m; /**< cut-off parameter       */
  int p; /**< degree of smoothness    */
  const char *s; /**< name of kernel          */
  C (*kernel)(R, int, const R *); /**< kernel function         */
  R c; /**< parameter for kernel    */
  fastsum_plan my_fastsum_plan; /**< plan for fast summation */
  C *direct; /**< array for direct computation */
  ticks t0, t1; /**< for time measurement    */
  R time; /**< for time measurement    */
  R error = K(0.0); /**< for error computation   */
  R eps_I; /**< inner boundary          */
  R eps_B; /**< outer boundary          */

  if (argc != 11)
  {
    printf("\nfastsum_test d N M n m p kernel c eps_I eps_B\n\n");
    printf("  d       dimension                 \n");
    printf("  N       number of source nodes    \n");
    printf("  M       number of target nodes    \n");
    printf("  n       expansion degree          \n");
    printf("  m       cut-off parameter         \n");
    printf("  p       degree of smoothness      \n");
    printf("  kernel  kernel function  (e.g., gaussian)\n");
    printf("  c       kernel parameter          \n");
    printf("  eps_I   inner boundary            \n");
    printf("  eps_B   outer boundary            \n\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    d = atoi(argv[1]);
    N = atoi(argv[2]);
    c = K(1.0) / POW((R)(N), K(1.0) / ((R)(d)));
    M = atoi(argv[3]);
    n = atoi(argv[4]);
    m = atoi(argv[5]);
    p = atoi(argv[6]);
    s = argv[7];
    c = (R)(atof(argv[8]));
    eps_I = (R)(atof(argv[9]));
    eps_B = (R)(atof(argv[10]));
    if (strcmp(s, "gaussian") == 0)
      kernel = gaussian;
    else if (strcmp(s, "multiquadric") == 0)
      kernel = multiquadric;
    else if (strcmp(s, "inverse_multiquadric") == 0)
      kernel = inverse_multiquadric;
    else if (strcmp(s, "logarithm") == 0)
      kernel = logarithm;
    else if (strcmp(s, "thinplate_spline") == 0)
      kernel = thinplate_spline;
    else if (strcmp(s, "one_over_square") == 0)
      kernel = one_over_square;
    else if (strcmp(s, "one_over_modulus") == 0)
      kernel = one_over_modulus;
    else if (strcmp(s, "one_over_x") == 0)
      kernel = one_over_x;
    else if (strcmp(s, "inverse_multiquadric3") == 0)
      kernel = inverse_multiquadric3;
    else if (strcmp(s, "sinc_kernel") == 0)
      kernel = sinc_kernel;
    else if (strcmp(s, "cosc") == 0)
      kernel = cosc;
    else if (strcmp(s, "cot") == 0)
      kernel = kcot;
    else if (strcmp(s, "one_over_cube") == 0)
      kernel = one_over_cube;
    else if (strcmp(s, "log_sin") == 0)
      kernel = log_sin;
    else if (strcmp(s, "laplacian_rbf") == 0)
      kernel = laplacian_rbf;
    else
    {
      s = "multiquadric";
      kernel = multiquadric;
    }
  }
  printf(
      "d=%d, N=%d, M=%d, n=%d, m=%d, p=%d, kernel=%s, c=%" __FGS__ ", eps_I=%" __FGS__ ", eps_B=%" __FGS__ " \n",
      d, N, M, n, m, p, s, c, eps_I, eps_B);
#ifdef NF_KUB
  printf("nearfield correction using piecewise cubic Lagrange interpolation\n");
#elif defined(NF_QUADR)
  printf("nearfield correction using piecewise quadratic Lagrange interpolation\n");
#elif defined(NF_LIN)
  printf("nearfield correction using piecewise linear Lagrange interpolation\n");
#endif

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp single
    {
      printf("nthreads=%d\n", omp_get_max_threads());
    }
  }

  FFTW(init_threads)();
#endif

  /** init d-dimensional fastsum plan */
  fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, 0, n, m, p, eps_I,
      eps_B);
  //fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, NEARFIELD_BOXES, n, m, p, eps_I, eps_B);

  if (my_fastsum_plan.flags & NEARFIELD_BOXES)
    printf(
        "determination of nearfield candidates based on partitioning into boxes\n");
  else
    printf("determination of nearfield candidates based on search tree\n");

  /** init source knots in a d-ball with radius 0.25-eps_b/2 */
  k = 0;
  while (k < N)
  {
    R r_max = K(0.25) - my_fastsum_plan.eps_B / K(2.0);
    R r2 = K(0.0);

    for (j = 0; j < d; j++)
      my_fastsum_plan.x[k * d + j] = K(2.0) * r_max * NFFT(drand48)() - r_max;

    for (j = 0; j < d; j++)
      r2 += my_fastsum_plan.x[k * d + j] * my_fastsum_plan.x[k * d + j];

    if (r2 >= r_max * r_max)
      continue;

    k++;
  }

  for (k = 0; k < N; k++)
  {
    /*    R r=(0.25-my_fastsum_plan.eps_B/2.0)*pow((R)rand()/(R)RAND_MAX,1.0/d);
     my_fastsum_plan.x[k*d+0] = r;
     for (j=1; j<d; j++)
     {
     R phi=2.0*KPI*(R)rand()/(R)RAND_MAX;
     my_fastsum_plan.x[k*d+j] = r;
     for (t=0; t<j; t++)
     {
     my_fastsum_plan.x[k*d+t] *= cos(phi);
     }
     my_fastsum_plan.x[k*d+j] *= sin(phi);
     }
     */
    my_fastsum_plan.alpha[k] = NFFT(drand48)() + II * NFFT(drand48)();
  }

  /** init target knots in a d-ball with radius 0.25-eps_b/2 */
  k = 0;
  while (k < M)
  {
    R r_max = K(0.25) - my_fastsum_plan.eps_B / K(2.0);
    R r2 = K(0.0);

    for (j = 0; j < d; j++)
      my_fastsum_plan.y[k * d + j] = K(2.0) * r_max * NFFT(drand48)() - r_max;

    for (j = 0; j < d; j++)
      r2 += my_fastsum_plan.y[k * d + j] * my_fastsum_plan.y[k * d + j];

    if (r2 >= r_max * r_max)
      continue;

    k++;
  }
  /*  for (k=0; k<M; k++)
   {
   R r=(0.25-my_fastsum_plan.eps_B/2.0)*pow((R)rand()/(R)RAND_MAX,1.0/d);
   my_fastsum_plan.y[k*d+0] = r;
   for (j=1; j<d; j++)
   {
   R phi=2.0*KPI*(R)rand()/(R)RAND_MAX;
   my_fastsum_plan.y[k*d+j] = r;
   for (t=0; t<j; t++)
   {
   my_fastsum_plan.y[k*d+t] *= cos(phi);
   }
   my_fastsum_plan.y[k*d+j] *= sin(phi);
   }
   } */

  /** direct computation */
  printf("direct computation: ");
  fflush(NULL);
  t0 = getticks();
  fastsum_exact(&my_fastsum_plan);
  t1 = getticks();
  time = NFFT(elapsed_seconds)(t1, t0);
  printf(__FI__ "sec\n", time);

  /** copy result */
  direct = (C *) NFFT(malloc)((size_t)(my_fastsum_plan.M_total) * (sizeof(C)));
  for (j = 0; j < my_fastsum_plan.M_total; j++)
    direct[j] = my_fastsum_plan.f[j];

  /** precomputation */
  printf("pre-computation:    ");
  fflush(NULL);
  t0 = getticks();
  fastsum_precompute(&my_fastsum_plan);
  t1 = getticks();
  time = NFFT(elapsed_seconds)(t1, t0);
  printf(__FI__ "sec\n", time);

  /** fast computation */
  printf("fast computation:   ");
  fflush(NULL);
  t0 = getticks();
  fastsum_trafo(&my_fastsum_plan);
  t1 = getticks();
  time = NFFT(elapsed_seconds)(t1, t0);
  printf(__FI__ "sec\n", time);

  /** compute max error */
  error = K(0.0);
  for (j = 0; j < my_fastsum_plan.M_total; j++)
  {
    if (CABS(direct[j] - my_fastsum_plan.f[j]) / CABS(direct[j]) > error)
      error = CABS(direct[j] - my_fastsum_plan.f[j]) / CABS(direct[j]);
  }
  printf("max relative error: %" __FES__ "\n", error);

  /** finalise the plan */
  fastsum_finalize(&my_fastsum_plan);

  return EXIT_SUCCESS;
}
/* \} */
