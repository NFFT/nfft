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
#include "wigner.h"

#define DEFAULT_NFFT_CUTOFF    6
#define FPT_THRESHOLD          1000

static fpt_set SO3_fpt_init(int l, unsigned int flags, int kappa);

void nfsoft_init(nfsoft_plan *plan, int N, int M)
{
  nfsoft_init_advanced(plan, N, M, NFSOFT_MALLOC_X | NFSOFT_MALLOC_F
      | NFSOFT_MALLOC_F_HAT);
}

void nfsoft_init_advanced(nfsoft_plan *plan, int N, int M,
    unsigned int nfsoft_flags)
{
  nfsoft_init_guru(plan, N, M, nfsoft_flags, PRE_PHI_HUT | PRE_PSI | MALLOC_X | NFFT_OMP_BLOCKWISE_ADJOINT
      | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
      DEFAULT_NFFT_CUTOFF, FPT_THRESHOLD);
}

void nfsoft_init_guru(nfsoft_plan *plan, int B, int M,
    unsigned int nfsoft_flags, unsigned int nfft_flags, int nfft_cutoff,
    int fpt_kappa)
{
	nfsoft_init_guru_advanced(plan, B, M, nfsoft_flags, nfft_flags, nfft_cutoff, fpt_kappa, 8* B);
}

void nfsoft_init_guru_advanced(nfsoft_plan *plan, int B, int M,
    unsigned int nfsoft_flags, unsigned int nfft_flags, int nfft_cutoff,
    int fpt_kappa, int nn_oversampled)
{
  int N[3];
  int n[3];

  N[0] = 2* B + 2;
  N[1] = 2* B + 2;
  N[2] = 2* B + 2;

  n[0] = nn_oversampled ;
  n[1] = nn_oversampled ;
  n[2] = nn_oversampled ;

  nfft_init_guru(&plan->p_nfft, 3, N, M, n, nfft_cutoff, nfft_flags,
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  if ((plan->p_nfft).flags & PRE_LIN_PSI)
  {
    nfft_precompute_lin_psi(&(plan->p_nfft));
  }

  plan->N_total = B;
  plan->M_total = M;
  plan->fpt_kappa = fpt_kappa;
  plan->flags = nfsoft_flags;

  if (plan->flags & NFSOFT_MALLOC_F_HAT)
  {
    plan->f_hat = (C*) nfft_malloc((B + 1) * (4* (B +1)*(B+1)-1)/3*sizeof(C));
    if (plan->f_hat == NULL ) printf("Allocation failed!\n");
  }

  if (plan->flags & NFSOFT_MALLOC_X)
  {
    plan->x = (R*) nfft_malloc(plan->M_total*3*sizeof(R));
    if (plan->x == NULL ) printf("Allocation failed!\n");
  }
  if (plan->flags & NFSOFT_MALLOC_F)
  {
    plan->f = (C*) nfft_malloc(plan->M_total*sizeof(C));
      if (plan->f == NULL ) printf("Allocation failed!\n");
  }

  plan->wig_coeffs = (C*) nfft_malloc((X(next_power_of_2)(B)+1)*sizeof(C));
  plan->cheby = (C*) nfft_malloc((2*B+2)*sizeof(C));
  plan->aux = (C*) nfft_malloc((2*B+4)*sizeof(C));

  if (plan->wig_coeffs == NULL ) printf("Allocation failed!\n");
  if (plan->cheby == NULL ) printf("Allocation failed!\n");
  if (plan->aux == NULL ) printf("Allocation failed!\n");

  plan->mv_trafo = (void (*) (void* ))nfsoft_trafo;
  plan->mv_adjoint = (void (*) (void* ))nfsoft_adjoint;

  plan->internal_fpt_set = SO3_fpt_init(plan->N_total, plan->flags, plan->fpt_kappa);

}

static void c2e(nfsoft_plan *my_plan, int even)
{
  int j, N;

  /**initialize the bigger plan*/
  N = 2* (my_plan ->N_total+1);

  /** prepare the coefficients for the new plan */
  my_plan->cheby[my_plan->N_total+1] = my_plan->wig_coeffs[0];
  my_plan->cheby[0]=0.0;

  for (j=1;j<my_plan->N_total+1;j++)
  {
    my_plan->cheby[my_plan->N_total+1+j]=0.5* my_plan->wig_coeffs[j];
    my_plan->cheby[my_plan->N_total+1-j]=0.5* my_plan->wig_coeffs[j];
  }

  C *aux= (C*) nfft_malloc((N+2)*sizeof(C));

  for(j=1;j<N;j++)
  aux[j]=my_plan->cheby[j];

  aux[0]=0.;
  aux[N]=0.;

  if (even>0)
  {
    my_plan->cheby[0]=(C) (-1.)/(2.*_Complex_I) * aux[1];
    for (j=1;j<N;j++)
    {
      my_plan->cheby[j]=(1./(2.*_Complex_I)*(aux[j+1]-aux[j-1]));
    }

  }
  nfft_free(aux);
  aux = NULL;
}


static fpt_set SO3_fpt_init(int l, unsigned int flags, int kappa)
{
  fpt_set set = 0;
  int N, t, k_start, k_end, k, m;
  int glo = 0;
  R *alpha, *beta, *gamma;

  /** Read in transfrom length. */
  if (flags & NFSOFT_USE_DPT)
  {
    if (l < 2)
      N = 2;
    else
      N = l;

    t = (int) log2(X(next_power_of_2)(N));

  }
  else
  {
    /** workaround to compute polynomials of degree less than 2*/
    if (l < 2)
      N = 2;
    else
      N = X(next_power_of_2)(l);

    t = (int) log2(N);
  }

  /**memory for the recurrence coefficients*/
  alpha = (R*) nfft_malloc((N + 2) * sizeof(R));
  beta = (R*) nfft_malloc((N + 2) * sizeof(R));
  gamma = (R*) nfft_malloc((N + 2) * sizeof(R));

  /** Initialize DPT. */
  if (flags & NFSOFT_NO_STABILIZATION)
  {
    set = fpt_init((2* N + 1) * (2* N + 1), t, 0U | FPT_NO_STABILIZATION);
  }
  else
  {
    set = fpt_init((2* N + 1) * (2* N + 1), t, 0U);
  }

  for (k = -N; k <= N; k++)
    for (m = -N; m <= N; m++)
    {
      /** Read in start and end indeces */
      k_start = (ABS(k) >= ABS(m)) ? ABS(k) : ABS(m);
      k_end = N;

      SO3_alpha_row(alpha, N, k, m);
      SO3_beta_row(beta, N, k, m);
      SO3_gamma_row(gamma, N, k, m);

      fpt_precompute(set, glo, alpha, beta, gamma, k_start, kappa);
      glo++;
    }

  nfft_free(alpha);
  nfft_free(beta);
  nfft_free(gamma);
  alpha = NULL;
  beta = NULL;
  gamma = NULL;

  return set;
}

static fpt_set SO3_single_fpt_init(int l, int k, int m, unsigned int flags, int kappa)
{
  int N, t, k_start, k_end;
  R *alpha, *beta, *gamma;
  fpt_set set = 0;

  /** Read in transfrom length. */
  if (flags & NFSOFT_USE_DPT)
  {
    if (l < 2)
      N = 2;
    else
      N = l;

    t = (int) log2(X(next_power_of_2)(N));

  }
  else
  {
    /** workaround to compute polynomials of degree less than 2*/
    if (l < 2)
      N = 2;
    else
      N = X(next_power_of_2)(l);

    t = (int) log2(N);
  }

  /**memory for the recurrence coefficients*/
  alpha = (R*) nfft_malloc((N + 2) * sizeof(R));
  beta = (R*) nfft_malloc((N + 2) * sizeof(R));
  gamma = (R*) nfft_malloc((N + 2) * sizeof(R));

  /** Initialize DPT. */
  {
    unsigned int fptflags = 0U 
      | IF(flags & NFSOFT_USE_DPT,FPT_NO_FAST_ALGORITHM,IF(t > 1,FPT_NO_DIRECT_ALGORITHM,0U))
      | IF(flags & NFSOFT_NO_STABILIZATION,FPT_NO_STABILIZATION,0U);
    set = fpt_init(1, t, fptflags);
  }

  /** Read in start and end indices */
  k_start = (ABS(k) >= ABS(m)) ? ABS(k) : ABS(m);
  k_end = N;

  SO3_alpha_row(alpha, N, k, m);
  SO3_beta_row(beta, N, k, m);
  SO3_gamma_row(gamma, N, k, m);
  
  /*{
    int rr;
    for (rr = 0; rr < N + 2; rr++)
      fprintf(stderr, "a[%4d] = %10e b[%4d] = %10e c[%4d] = %10e\n",rr,alpha[rr],rr,beta[rr],rr,gamma[rr]);
  }*/

  fpt_precompute(set, 0, alpha, beta, gamma, k_start, kappa);

  nfft_free(alpha);
  nfft_free(beta);
  nfft_free(gamma);
  alpha = NULL;
  beta = NULL;
  gamma = NULL;

  return set;
}

void SO3_fpt(C *coeffs, fpt_set set, int l, int k, int m, unsigned int flags)
{
  int N;
  /** The Wigner  coefficients */
  C* x;
  /** The Chebyshev coefficients */
  C* y;

  int trafo_nr; /**gives the index of the trafo in the FPT_set*/
  int k_start, k_end, j;
  int function_values = 0;

  /** Read in transfrom length. */
  if (flags & NFSOFT_USE_DPT)
  {
    N = l;
    if (l < 2)
      N = 2;
  }
  else
  {
    if (l < 2)
      N = 2;
    else
      N = X(next_power_of_2)(l);
  }

  /** Read in start and end indeces */
  k_start = (ABS(k) >= ABS(m)) ? ABS(k) : ABS(m);
  k_end = N;
  trafo_nr = (N + k) * (2* N + 1) + (m + N);

  /** Read in Wigner coefficients. */
  x = (C*) nfft_malloc((k_end + 1) * sizeof(C));

  for (j = 0; j <= k_end; j++)
   x[j] = K(0.0);


  for (j = 0; j <= l - k_start; j++)
  {
    x[j + k_start] = coeffs[j];
  }
  for (j = l - k_start + 1; j <= k_end - k_start; j++)
  {
    x[j + k_start] = K(0.0);
  }

  /** Allocate memory for Chebyshev coefficients. */
  y = (C*) nfft_malloc((k_end + 1) * sizeof(C));

  if (flags & NFSOFT_USE_DPT)
  { /** Execute DPT. */
    fpt_trafo_direct(set, trafo_nr, &x[k_start], y, k_end, 0U
        | (function_values ? FPT_FUNCTION_VALUES : 0U));
  }
  else
  { /** compute fpt*/
    fpt_trafo(set, trafo_nr, &x[k_start], y, k_end, 0U
        | (function_values ? FPT_FUNCTION_VALUES : 0U));
  }

  /**write computed coeffs in the plan*/
  for (j = 0; j <= l; j++)
  {
    coeffs[j] = y[j];
  }

  /** Free memory. */

  nfft_free(x);
  nfft_free(y);
  x = NULL;
  y = NULL;
}

void SO3_fpt_transposed(C *coeffs, fpt_set set, int l, int k, int m,
    unsigned int flags)
{
  int N, k_start, k_end, j;
  int trafo_nr; /**gives the index of the trafo in the FPT_set*/
  int function_values = 0;
  /** The Wigner  coefficients */
  C* x;
  /** The Chebyshev coefficients */
  C* y;

  /** Read in transfrom length. */

  if (flags & NFSOFT_USE_DPT)
  {
    N = l;
    if (l < 2)
      N = 2;
  }
  else
  {
    if (l < 2)
      N = 2;
    else
      N = X(next_power_of_2)(l);
  }

  /** Read in start and end indeces */
  k_start = (ABS(k) >= ABS(m)) ? ABS(k) : ABS(m);
  k_end = N;
  trafo_nr = (N + k) * (2* N + 1) + (m + N);

  /** Allocate memory for Chebychev coefficients. */
  y = (C*) nfft_malloc((k_end + 1) * sizeof(C));
  /** Allocate memory for Wigner coefficients. */
  x = (C*) nfft_malloc((k_end + 1) * sizeof(C));

  for (j = 0; j <= l; j++)
  {
    y[j] = coeffs[j];
  }
  for (j = l + 1; j <= k_end; j++)
  {
    y[j] = K(0.0);
  }

  if (flags & NFSOFT_USE_DPT)
  {
    fpt_transposed_direct(set, trafo_nr, &x[k_start], y, k_end, 0U
        | (function_values ? FPT_FUNCTION_VALUES : 0U));
  }
  else
  {
    fpt_transposed(set, trafo_nr, &x[k_start], y, k_end, 0U
        | (function_values ? FPT_FUNCTION_VALUES : 0U));
  }

  for (j = 0; j <= l; j++)
  {
    coeffs[j] = x[j];
  }

  /** Free memory. */
  nfft_free(x);
  nfft_free(y);
  x = NULL;
  y = NULL;
}

void nfsoft_precompute(nfsoft_plan *plan3D)
{
  int j;
  int N = plan3D->N_total;
  int M = plan3D->M_total;

  /** Node-dependent part*/

  for (j = 0; j < M; j++)
  {
    plan3D->p_nfft.x[3* j ] = plan3D->x[3* j + 2];
    plan3D->p_nfft.x[3* j + 1] = plan3D->x[3* j ];
    plan3D->p_nfft.x[3* j + 2] = plan3D->x[3* j + 1];
  }

  for (j = 0; j < 3* plan3D ->p_nfft.M_total; j++)
  {
    plan3D->p_nfft.x[j] = plan3D->p_nfft.x[j] * (1 / (2* KPI ));
  }

  if ((plan3D->p_nfft).flags & FG_PSI)
  {
    nfft_precompute_one_psi(&(plan3D->p_nfft));
  }
  if ((plan3D->p_nfft).flags & PRE_PSI)
  {
    nfft_precompute_one_psi(&(plan3D->p_nfft));
  }

}

void nfsoft_trafo(nfsoft_plan *plan3D)
{
  int i, j, m, k, max, glo1, glo2;

  i = 0;
  glo1 = 0;
  glo2 = 0;

  int N = plan3D->N_total;
  int M = plan3D->M_total;

  /**almost nothing to be done for polynomial degree 0*/
  if (N == 0)
  {
    for (j = 0; j < M; j++)
      plan3D->f[j] = plan3D->f_hat[0];
    return;
  }

  for (j = 0; j < plan3D->p_nfft.N_total; j++)
    plan3D->p_nfft.f_hat[j] = 0.0;

  for (k = -N; k <= N; k++)
  {
    for (m = -N; m <= N; m++)
    {

      max = (ABS(m) > ABS(k) ? ABS(m) : ABS(k));

      for (j = 0; j <= N - max; j++)
      {
        plan3D->wig_coeffs[j] = plan3D->f_hat[glo1];

        if ((plan3D->flags & NFSOFT_NORMALIZED))
        {
          plan3D->wig_coeffs[j] = plan3D->wig_coeffs[j] * (1. / (2. * KPI))
              * SQRT(0.5 * (2. * (max + j) + 1.));
        }

        if ((plan3D->flags & NFSOFT_REPRESENT))
        {
          if ((k < 0) && (k % 2))
          {
            plan3D->wig_coeffs[j] = plan3D->wig_coeffs[j] * (-1);
          }
          if ((m < 0) && (m % 2))
            plan3D->wig_coeffs[j] = plan3D->wig_coeffs[j] * (-1);

          if ((m + k) % 2)
	    plan3D->wig_coeffs[j] = plan3D->wig_coeffs[j] * (-1);

        }

        glo1++;
      }

      for (j = N - max + 1; j < X(next_power_of_2)(N) + 1; j++)
        plan3D->wig_coeffs[j] = 0.0;
      //fprintf(stdout,"\n k= %d, m= %d \n",k,m);
      SO3_fpt(plan3D->wig_coeffs, plan3D->internal_fpt_set, N, k, m, plan3D->flags);

      c2e(plan3D, ABS((k + m) % 2));

      for (i = 1; i <= 2* plan3D ->N_total + 2; i++)
      {
        plan3D->p_nfft.f_hat[NFSOFT_INDEX(k, m, i - N - 1, N) - 1]
            = plan3D->cheby[i - 1];
        //fprintf(stdout,"%f \t", plan3D->nfft_plan.f_hat[NFSOFT_INDEX(k,m,i-N-1,N)-1]);
        //fprintf(stdout,"another index: %d for k=%d,m=%d,l=%d,N=%d \n", NFSOFT_INDEX(k,m,i-N-1,N)-1,k,m,i-N-1,N);
      }

    }
  }

  if (plan3D->flags & NFSOFT_USE_NDFT)
  {
    nfft_trafo_direct(&(plan3D->p_nfft));
  }
  else
  {
    nfft_trafo(&(plan3D->p_nfft));
  }

  for (j = 0; j < plan3D->M_total; j++)
    plan3D->f[j] = plan3D->p_nfft.f[j];

}

static void e2c(nfsoft_plan *my_plan, int even)
{
  int N;
  int j;

  /**initialize the bigger plan*/
  N = 2* (my_plan ->N_total+1);
  //nfft_vpr_complex(my_plan->cheby,N+1,"chebychev");


      if (even>0)
      {
        //my_plan->aux[N-1]= -1/(2*I)* my_plan->cheby[N-2];
        my_plan->aux[0]= 1/(2*_Complex_I)*my_plan->cheby[1];

        for(j=1;j<N-1;j++)
        {
          my_plan->aux[j]=1/(2*_Complex_I)*(my_plan->cheby[j+1]-my_plan->cheby[j-1]);
}
my_plan->aux[N-1]=1/(2*_Complex_I)*(-my_plan->cheby[j-1]);


for(j=0;j<N;j++)
{
my_plan->cheby[j]= my_plan->aux[j];
}
}

my_plan->wig_coeffs[0]=my_plan->cheby[my_plan->N_total+1];

for(j=1;j<=my_plan->N_total;j++)
{
my_plan->wig_coeffs[j]=0.5*(my_plan->cheby[my_plan->N_total+j+1]+my_plan->cheby[my_plan->N_total+1-j]);
}



//nfft_vpr_complex(my_plan->wig_coeffs,my_plan->N_total,"chebychev ");

}

void nfsoft_adjoint(nfsoft_plan *plan3D)
{
  int i, j, m, k, max, glo1, glo2;

  i = 0;
  glo1 = 0;
  glo2 = 0;

  int N = plan3D->N_total;
  int M = plan3D->M_total;

  //nothing much to be done for polynomial degree 0
  if (N == 0)
  {
    plan3D->f_hat[0]=0;
    for (j = 0; j < M; j++)
      plan3D->f_hat[0] += plan3D->f[j];
    return;
  }

  for (j = 0; j < M; j++)
  {
    plan3D->p_nfft.f[j] = plan3D->f[j];
  }

  if (plan3D->flags & NFSOFT_USE_NDFT)
  {
    nfft_adjoint_direct(&(plan3D->p_nfft));
  }
  else
  {
    nfft_adjoint(&(plan3D->p_nfft));
  }

  //nfft_vpr_complex(plan3D->nfft_plan.f_hat,plan3D->nfft_plan.N_total,"all results");

  glo1 = 0;

  for (k = -N; k <= N; k++)
  {
    for (m = -N; m <= N; m++)
    {

      max = (ABS(m) > ABS(k) ? ABS(m) : ABS(k));

      for (i = 1; i < 2* plan3D ->N_total + 3; i++)
      {
        plan3D->cheby[i - 1] = plan3D->p_nfft.f_hat[NFSOFT_INDEX(k, m, i - N
            - 1, N) - 1];
      }

      //fprintf(stdout,"k=%d,m=%d \n",k,m);
      //nfft_vpr_complex(plan3D->cheby,2*plan3D->N_total+2,"euler");
      e2c(plan3D, ABS((k + m) % 2));

      //nfft_vpr_complex(plan3D->wig_coeffs,plan3D->N_total+1,"chebys");
      SO3_fpt_transposed(plan3D->wig_coeffs, plan3D->internal_fpt_set, N, k, m,
          plan3D->flags);
      //nfft_vpr_complex(plan3D->wig_coeffs,plan3D->N_total+1,"wigners");
      //  SO3_fpt_transposed(plan3D->wig_coeffs,N,k,m,plan3D->flags,plan3D->fpt_kappa);


      for (j = max; j <= N; j++)
      {
        if ((plan3D->flags & NFSOFT_REPRESENT))
        {
          if ((k < 0) && (k % 2))
          {
            plan3D->wig_coeffs[j] = -plan3D->wig_coeffs[j];
          }
          if ((m < 0) && (m % 2))
            plan3D->wig_coeffs[j] = -plan3D->wig_coeffs[j];

          if ((m + k) % 2)
            plan3D->wig_coeffs[j] = plan3D->wig_coeffs[j] * (-1);

        }

        plan3D->f_hat[glo1] = plan3D->wig_coeffs[j];

        if ((plan3D->flags & NFSOFT_NORMALIZED))
        {
          plan3D->f_hat[glo1] = plan3D->f_hat[glo1] * (1 / (2. * KPI)) * SQRT(
              0.5 * (2. * (j) + 1.));
        }

        glo1++;
      }

    }
  }
}

void nfsoft_finalize(nfsoft_plan *plan)
{
  /* Finalise the nfft plan. */
  nfft_finalize(&plan->p_nfft);
  nfft_free(plan->wig_coeffs);
  nfft_free(plan->cheby);
  nfft_free(plan->aux);

  fpt_finalize(plan->internal_fpt_set);
  plan->internal_fpt_set = NULL;

  if (plan->flags & NFSOFT_MALLOC_F_HAT)
  {
    //fprintf(stderr,"deallocating f_hat\n");
    nfft_free(plan->f_hat);
  }

  /* De-allocate memory for samples, if neccesary. */
  if (plan->flags & NFSOFT_MALLOC_F)
  {
    //fprintf(stderr,"deallocating f\n");
    nfft_free(plan->f);
  }

  /* De-allocate memory for nodes, if neccesary. */
  if (plan->flags & NFSOFT_MALLOC_X)
  {
    //fprintf(stderr,"deallocating x\n");
    nfft_free(plan->x);
  }
}

int posN(int n, int m, int B)
{
  int pos;

  if (n > -B)
    pos = posN(n - 1, m, B) + B + 1 - MAX(ABS(m), ABS(n - 1));
  else
    pos = 0;
  //(n > -B? pos=posN(n-1,m,B)+B+1-MAX(ABS(m),ABS(n-1)): pos= 0)
  return pos;
}

