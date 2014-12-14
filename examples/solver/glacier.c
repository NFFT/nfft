/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id$ */
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

#undef X
#define X(name) NFFT(name)

/**
 * \defgroup examples_solver_glacier Reconstruction of a glacier from \
 scattered data
 * \ingroup examples_solver
 * \{
 */

/** Generalised Sobolev weight */
static R my_weight(R z, R a, R b, R c)
{
  return POW(K(0.25) - z * z, b) / (c + POW(FABS(z), K(2.0) * a));
}

/** Reconstruction routine */
static void glacier(int N, int M)
{
  int j, k, k0, k1, l, my_N[2], my_n[2];
  R tmp_y;
  X(plan) p;
  SOLVER(plan_complex) ip;
  FILE* fp;

  /* initialise p */
  my_N[0] = N;
  my_n[0] = (int)X(next_power_of_2)(N);
  my_N[1] = N;
  my_n[1] = (int)X(next_power_of_2)(N);
  X(init_guru)(&p, 2, my_N, M, my_n, 6,
  PRE_PHI_HUT | PRE_FULL_PSI |
  MALLOC_X | MALLOC_F_HAT | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  /* initialise ip, specific */
  SOLVER(init_advanced_complex)(&ip, (X(mv_plan_complex)*) (&p),
      CGNE | PRECOMPUTE_DAMP);
  fprintf(stderr, "Using the generic solver!");

  /* init nodes */
  fp = fopen("input_data.dat", "r");
  for (j = 0; j < p.M_total; j++)
  {
    fscanf(fp, "%" __FES__ " %" __FES__ " %" __FES__, &p.x[2 * j + 0], &p.x[2 * j + 1], &tmp_y);
    ip.y[j] = tmp_y;
  }
  fclose(fp);

  /* precompute psi */
  if (p.flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&p);

  /* initialise damping factors */
  if (ip.flags & PRECOMPUTE_DAMP)
    for (k0 = 0; k0 < p.N[0]; k0++)
      for (k1 = 0; k1 < p.N[1]; k1++)
        ip.w_hat[k0 * p.N[1] + k1] = my_weight(((R) ((R)(k0) - (R)(p.N[0]) / K(2.0))) / ((R)(p.N[0])),
            K(0.5), K(3.0), K(0.001))
            * my_weight((((R)(k1) - (R)(p.N[1]) / K(2.0))) / ((R)(p.N[1])), K(0.5), K(3.0), K(0.001));

  /* init some guess */
  for (k = 0; k < p.N_total; k++)
    ip.f_hat_iter[k] = K(0.0);

  /* inverse trafo */
  SOLVER(before_loop_complex)(&ip);

  for (l = 0; l < 40; l++)
  {
    fprintf(stderr, "Residual ||r||=%" __FES__ ",\n", SQRT(ip.dot_r_iter));
    SOLVER(loop_one_step_complex)(&ip);
  }

  for (k = 0; k < p.N_total; k++)
    printf("%" __FES__ " %" __FES__ "\n", CREAL(ip.f_hat_iter[k]), CIMAG(ip.f_hat_iter[k]));

  SOLVER(finalize_complex)(&ip);
  X(finalize)(&p);
}

/** Reconstruction routine with cross validation */
static void glacier_cv(int N, int M, int M_cv, unsigned solver_flags)
{
  int j, k, k0, k1, l, my_N[2], my_n[2];
  R tmp_y, r;
  X(plan) p, cp;
  SOLVER(plan_complex) ip;
  C* cp_y;
  FILE* fp;
  int M_re = M - M_cv;

  /* initialise p for reconstruction */
  my_N[0] = N;
  my_n[0] = (int)X(next_power_of_2)(N);
  my_N[1] = N;
  my_n[1] = (int)X(next_power_of_2)(N);
  X(init_guru)(&p, 2, my_N, M_re, my_n, 6,
  PRE_PHI_HUT | PRE_FULL_PSI |
  MALLOC_X | MALLOC_F_HAT | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  /* initialise ip, specific */
  SOLVER(init_advanced_complex)(&ip, (X(mv_plan_complex)*) (&p), solver_flags);

  /* initialise cp for validation */
  cp_y = (C*) X(malloc)((size_t)(M) * sizeof(C));
  X(init_guru)(&cp, 2, my_N, M, my_n, 6,
  PRE_PHI_HUT | PRE_FULL_PSI |
  MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  cp.f_hat = ip.f_hat_iter;

  /* set up data in cp and cp_y */
  fp = fopen("input_data.dat", "r");
  for (j = 0; j < cp.M_total; j++)
  {
    fscanf(fp, "%" __FES__ " %" __FES__ " %" __FES__, &cp.x[2 * j + 0], &cp.x[2 * j + 1], &tmp_y);
    cp_y[j] = tmp_y;
  }
  fclose(fp);

  /* copy part of the data to p and ip */
  for (j = 0; j < p.M_total; j++)
  {
    p.x[2 * j + 0] = cp.x[2 * j + 0];
    p.x[2 * j + 1] = cp.x[2 * j + 1];
    ip.y[j] = tmp_y;
  }

  /* precompute psi */
  if (p.flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&p);

  /* precompute psi */
  if (cp.flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&cp);

  /* initialise damping factors */
  if (ip.flags & PRECOMPUTE_DAMP)
    for (k0 = 0; k0 < p.N[0]; k0++)
      for (k1 = 0; k1 < p.N[1]; k1++)
        ip.w_hat[k0 * p.N[1] + k1] = my_weight((((R)(k0) - (R)(p.N[0]) / K(2.0))) / ((R)(p.N[0])),
            K(0.5), K(3.0), K(0.001))
            * my_weight((((R)(k1) - (R)(p.N[1]) / K(2.0))) / ((R)(p.N[1])), K(0.5), K(3.0), K(0.001));

  /* init some guess */
  for (k = 0; k < p.N_total; k++)
    ip.f_hat_iter[k] = K(0.0);

  /* inverse trafo */
  SOLVER(before_loop_complex)(&ip);
  //  fprintf(stderr,"iteration starts,\t");
  for (l = 0; l < 40; l++)
    SOLVER(loop_one_step_complex)(&ip);

  //fprintf(stderr,"r=%1.2e, ",sqrt(ip.dot_r_iter)/M_re);

  CSWAP(p.f_hat, ip.f_hat_iter);
  X(trafo)(&p);
  CSWAP(p.f_hat, ip.f_hat_iter);
  Y(upd_axpy_complex)(p.f, -1, ip.y, M_re);
  r = SQRT(Y(dot_complex)(p.f, M_re) / Y(dot_complex)(cp_y, M));
  fprintf(stderr, "r=%1.2" __FES__ ", ", r);
  printf("$%1.1" __FES__ "$ & ", r);

  X(trafo)(&cp);
  Y(upd_axpy_complex)(&cp.f[M_re], -1, &cp_y[M_re], M_cv);
  r = SQRT(Y(dot_complex)(&cp.f[M_re], M_cv) / Y(dot_complex)(cp_y, M));
  fprintf(stderr, "r_1=%1.2" __FES__ "\t", r);
  printf("$%1.1" __FES__ "$ & ", r);

  X(finalize)(&cp);
  SOLVER(finalize_complex)(&ip);
  X(finalize)(&p);
}

/** Main routine */
int main(int argc, char **argv)
{
  int M_cv;

  if (argc < 3)
  {
    fprintf(stderr, "Call this program from the Matlab script glacier.m!");
    return EXIT_FAILURE;
  }

  if (argc == 3)
    glacier(atoi(argv[1]), atoi(argv[2]));
  else
    for (M_cv = atoi(argv[3]); M_cv <= atoi(argv[5]); M_cv += atoi(argv[4]))
    {
      fprintf(stderr, "\nM_cv=%d,\t", M_cv);
      printf("$%d$ & ", M_cv);
      fprintf(stderr, "cgne+damp: ");
      glacier_cv(atoi(argv[1]), atoi(argv[2]), M_cv, CGNE | PRECOMPUTE_DAMP);
      //fprintf(stderr,"cgne: ");
      //glacier_cv(atoi(argv[1]),atoi(argv[2]),M_cv,CGNE);
      fprintf(stderr, "cgnr: ");
      glacier_cv(atoi(argv[1]), atoi(argv[2]), M_cv, CGNR);
      fprintf(stderr, "cgnr: ");
      glacier_cv(atoi(argv[1]) / 4, atoi(argv[2]), M_cv, CGNR);
      printf("XXX \\\\\n");
    }

  fprintf(stderr, "\n");

  return EXIT_SUCCESS;
}
/* \} */
