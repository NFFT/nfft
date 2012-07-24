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

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"

/* void simple_test_infft_1d(int N, int M, int iter) */
/* { */
/*   int k,l;                            /\**< index for nodes, freqencies,iter*\/ */
/*   nfft_plan p;                          /\**< plan for the nfft               *\/ */
/*   infft_plan ip;                        /\**< plan for the inverse nfft       *\/ */

/*   /\** initialise an one dimensional plan *\/ */
/*   nfft_init_1d(&p, N, M); */

/*   /\** init pseudo random nodes *\/ */
/*   nfft_vrand_shifted_unit_double(p.x,p.M_total); */

/*   /\** precompute psi, the entries of the matrix B *\/ */
/*   if(p.nfft_flags & PRE_ONE_PSI) */
/*     nfft_precompute_one_psi(&p); */

/*   /\** initialise inverse plan *\/ */
/*   infft_init(&ip,&p); */

/*   /\** init pseudo random samples and show them *\/ */
/*   nfft_vrand_unit_complex(ip.y,p.M_total); */
/*   nfft_vpr_complex(ip.y,p.M_total,"Given data, vector y"); */

/*   /\** initialise some guess f_hat_0 and solve *\/ */
/*   for(k=0;k<p.N_total;k++) */
/*       ip.f_hat_iter[k]=0; */

/*   nfft_vpr_complex(ip.f_hat_iter,p.N_total,"Initial guess, vector f_hat_iter"); */

/*   NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat); */
/*   nfft_trafo(&p); */
/*   nfft_vpr_complex(p.f,p.M_total,"Data fit, vector f"); */
/*   NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat); */

/*   infft_before_loop(&ip); */
/*   printf("\n Residual r=%e\n",ip.dot_r_iter); */

/*   for(l=0;l<iter;l++) */
/*     { */
/*       printf("\n********** Iteration l=%d **********\n",l); */
/*       infft_loop_one_step(&ip); */
/*       nfft_vpr_complex(ip.f_hat_iter,p.N_total, */
/* 		  "Approximate solution, vector f_hat_iter"); */

/*       NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat); */
/*       nfft_trafo(&p); */
/*       nfft_vpr_complex(p.f,p.M_total,"Data fit, vector f"); */
/*       NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat); */

/*       printf("\n Residual r=%e\n",ip.dot_r_iter); */
/*     } */

/*   infft_finalize(&ip); */
/*   nfft_finalize(&p); */
/* } */

/** Simple test routine for the inverse nfft */
static void simple_test_solver_nfft_1d(int N, int M, int iter)
{
  int k,l;                            /**< index for nodes, freqencies,iter*/
  nfft_plan p;                          /**< plan for the nfft               */
  solver_plan_complex ip;                        /**< plan for the inverse nfft       */

  /** initialise an one dimensional plan */
  nfft_init_1d(&p, N, M);

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(p.x,p.M_total);

  /** precompute psi, the entries of the matrix B */
  if(p.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&p);

  /** initialise inverse plan */
  solver_init_complex(&ip,(nfft_mv_plan_complex*)(&p));

  /** init pseudo random samples and show them */
  nfft_vrand_unit_complex(ip.y,p.M_total);
  nfft_vpr_complex(ip.y,p.M_total,"Given data, vector y");

  /** initialise some guess f_hat_0 and solve */
  for(k=0;k<p.N_total;k++)
      ip.f_hat_iter[k]=0;

  nfft_vpr_complex(ip.f_hat_iter,p.N_total,"Initial guess, vector f_hat_iter");

  NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat);
  nfft_trafo(&p);
  nfft_vpr_complex(p.f,p.M_total,"Data fit, vector f");
  NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat);

  solver_before_loop_complex(&ip);
  printf("\n Residual r=%e\n",ip.dot_r_iter);

  for(l=0;l<iter;l++)
    {
      printf("\n********** Iteration l=%d **********\n",l);
      solver_loop_one_step_complex(&ip);
      nfft_vpr_complex(ip.f_hat_iter,p.N_total,
		  "Approximate solution, vector f_hat_iter");

      NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat);
      nfft_trafo(&p);
      nfft_vpr_complex(p.f,p.M_total,"Data fit, vector f");
      NFFT_SWAP_complex(ip.f_hat_iter,p.f_hat);

      printf("\n Residual r=%e\n",ip.dot_r_iter);
    }

  solver_finalize_complex(&ip);
  nfft_finalize(&p);
}

/** Main routine */
int main(void)
{
  printf("\n Computing a one dimensional inverse nfft\n");

  simple_test_solver_nfft_1d(8,4,5);

  return 1;
}
/* \} */
