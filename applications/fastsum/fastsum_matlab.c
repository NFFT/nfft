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

/*! \file fastsum_matlab.c
 *  \brief Simple test program for the fast NFFT-based summation algorithm, called by fastsum.m.
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

#include "fastsum.h"
#include "kernels.h"
#include "infft.h"

/**
 * \defgroup applications_fastsum_matlab fastsum_matlab
 * \ingroup applications_fastsum
 * \{
 */

int main(int argc, char **argv)
{
  int j, k, t; /**< indices                 */
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
  FILE *fid1, *fid2;
  R temp;

  if (argc != 11)
  {
    printf("\nfastsum_test d N M n m p kernel c\n\n");
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
    exit(-1);
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
      printf("Unrecognized kernel function!\n");
      exit(EXIT_FAILURE);
    }
  }
  printf(
      "d=%d, N=%d, M=%d, n=%d, m=%d, p=%d, kernel=%s, c=%" __FGS__ ", eps_I=%" __FGS__ ", eps_B=%" __FGS__ " \n",
      d, N, M, n, m, p, s, c, eps_I, eps_B);

  /** init two dimensional fastsum plan */
  fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, 0, n, m, p, eps_I,
      eps_B);
  /*fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, EXACT_NEARFIELD, n, m, p);*/

  /** load source knots and coefficients */
  fid1 = fopen("x.dat", "r");
  fid2 = fopen("alpha.dat", "r");
  for (k = 0; k < N; k++)
  {
    for (t = 0; t < d; t++)
    {
      fscanf(fid1, __FR__, &my_fastsum_plan.x[k * d + t]);
    }
    fscanf(fid2, __FR__, &temp);
    my_fastsum_plan.alpha[k] = temp;
    fscanf(fid2, __FR__, &temp);
    my_fastsum_plan.alpha[k] += temp * II;
  }
  fclose(fid1);
  fclose(fid2);

  /** load target knots */
  fid1 = fopen("y.dat", "r");
  for (j = 0; j < M; j++)
  {
    for (t = 0; t < d; t++)
    {
      fscanf(fid1, __FR__, &my_fastsum_plan.y[j * d + t]);
    }
  }
  fclose(fid1);

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
  printf("max relative error: " __FE__ "\n", error);

  /** write result to file */
  fid1 = fopen("f.dat", "w+");
  fid2 = fopen("f_direct.dat", "w+");
  if (fid1 == NULL)
  {
    printf("Error writing to file f.dat!\n");
    exit(EXIT_FAILURE);
  }
  for (j = 0; j < M; j++)
  {
    temp = CREAL(my_fastsum_plan.f[j]);
    fprintf(fid1, "  % .16" __FES__ "", temp);
    temp = CIMAG(my_fastsum_plan.f[j]);
    fprintf(fid1, "  % .16" __FES__ "\n", temp);

    temp = CREAL(direct[j]);
    fprintf(fid2, "  % .16" __FES__ "", temp);
    temp = CIMAG(direct[j]);
    fprintf(fid2, "  % .16" __FES__ "\n", temp);
  }
  fclose(fid1);
  fclose(fid2);

  /** finalise the plan */
  fastsum_finalize(&my_fastsum_plan);

  return EXIT_SUCCESS;
}
/* \} */
