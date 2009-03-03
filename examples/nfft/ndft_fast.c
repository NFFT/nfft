/*
 * $Id$
 *
 * Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
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

/*! \file ndft_fast.c
 *
 * \brief Testing ndft, Horner-like ndft, and fully precomputed ndft.
 *
 * \author Stefan Kunis
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "util.h"
#include "nfft3.h"

void ndft_horner_trafo(nfft_plan *ths)
{
  int j,k;
  double _Complex *f_hat_k, *f_j;
  double _Complex exp_omega_0;

  for(j=0, f_j=ths->f; j<ths->M_total; j++, f_j++)
    (*f_j) =0;

  for(j=0, f_j=ths->f; j<ths->M_total; j++, f_j++)
    {
      exp_omega_0 = cexp(+2*PI*_Complex_I*ths->x[j]);
      for(k=0, f_hat_k= ths->f_hat; k<ths->N[0]; k++, f_hat_k++)
        {
          (*f_j)+=(*f_hat_k);
          (*f_j)*=exp_omega_0;
	}
      (*f_j)*=cexp(-PI*_Complex_I*ths->N[0]*ths->x[j]);
    }
} /* ndft_horner_trafo */

void ndft_pre_full_trafo(nfft_plan *ths, double _Complex *A)
{
  int j,k;
  double _Complex *f_hat_k, *f_j;
  double _Complex *A_local;

  for(j=0, f_j=ths->f; j<ths->M_total; j++, f_j++)
    (*f_j) =0;

  for(j=0, f_j=ths->f, A_local=A; j<ths->M_total; j++, f_j++)
    for(k=0, f_hat_k= ths->f_hat; k<ths->N[0]; k++, f_hat_k++, A_local++)
      (*f_j) += (*f_hat_k)*(*A_local);
} /* ndft_pre_full_trafo */

void ndft_pre_full_init(nfft_plan *ths, double _Complex *A)
{
  int j,k;
  double _Complex *f_hat_k, *f_j, *A_local;

  for(j=0, f_j=ths->f, A_local=A; j<ths->M_total; j++, f_j++)
    for(k=0, f_hat_k= ths->f_hat; k<ths->N[0]; k++, f_hat_k++, A_local++)
      (*A_local) = cexp(-2*PI*_Complex_I*(k-ths->N[0]/2)*ths->x[j]);

} /* ndft_pre_full_init */

void ndft_time(int N, int M, unsigned test_ndft, unsigned test_pre_full)
{
  int r;
  double t, t_ndft, t_horner, t_pre_full, t_nfft;
  double _Complex *A;

  nfft_plan np;

  printf("%d\t%d\t",N, M);

  nfft_init_1d(&np, N, M);

  /** init pseudo random nodes */
  nfft_vrand_shifted_unit_double(np.x, np.M_total);

  if(test_pre_full)
   {
     A=(double _Complex*)nfft_malloc(N*M*sizeof(double _Complex));
     ndft_pre_full_init(&np, A);
   }

  /** init pseudo random Fourier coefficients */
  nfft_vrand_unit_complex(np.f_hat, np.N_total);

  /** NDFT */
  if(test_ndft)
    {
      t_ndft=0;
      r=0;
      while(t_ndft<0.1)
        {
          r++;
          t=nfft_second();
          ndft_trafo(&np);
          t=nfft_second()-t;
          t_ndft+=t;
        }
      t_ndft/=r;

      printf("%.2e\t",t_ndft);
    }
  else
    printf("nan\t\t");

  /** Horner NDFT */
  t_horner=0;
  r=0;
  while(t_horner<0.1)
    {
      r++;
      t=nfft_second();
      ndft_horner_trafo(&np);
      t=nfft_second()-t;
      t_horner+=t;
    }
  t_horner/=r;

  printf("%.2e\t", t_horner);

  /** Fully precomputed NDFT */
  if(test_pre_full)
    {
      t_pre_full=0;
      r=0;
      while(t_pre_full<0.1)
        {
          r++;
          t=nfft_second();
          ndft_pre_full_trafo(&np,A);
          t=nfft_second()-t;
          t_pre_full+=t;
        }
      t_pre_full/=r;

      printf("%.2e\t", t_pre_full);
    }
  else
    printf("nan\t\t");

  t_nfft=0;
  r=0;
  while(t_nfft<0.1)
    {
      r++;
      t=nfft_second();
      nfft_trafo(&np);
      t=nfft_second()-t;
      t_nfft+=t;
    }
  t_nfft/=r;

  printf("%.2e\n", t_nfft);

  fflush(stdout);

  if(test_pre_full)
    nfft_free(A);

  nfft_finalize(&np);
}

int main(int argc,char **argv)
{
  int l,trial;

  if(argc<=2)
    {
      fprintf(stderr,"ndft_fast type first last trials\n");
      return -1;
    }

  fprintf(stderr,"Testing ndft, Horner-like ndft, fully precomputed ndft.\n");
  fprintf(stderr,"Columns: N, M, t_ndft, t_horner, t_pre_full, t_nfft\n\n");

  /* time vs. N=M */
  if(atoi(argv[1])==0)
    {
      for(l=atoi(argv[2]); l<=atoi(argv[3]); l++)
        for(trial=0; trial<atoi(argv[4]); trial++)
          if(l<=13)
            ndft_time((1U<< l), (1U<< l), 1, 1);
          else
            if(l<=15)
              ndft_time((1U<< l), (1U<< l), 1, 0);
            else
              ndft_time((1U<< l), (1U<< l), 0, 0);
    }

  return 1;
}
