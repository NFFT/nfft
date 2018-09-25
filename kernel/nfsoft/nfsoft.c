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
#include "../fpt/fpt.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define DEFAULT_NFFT_CUTOFF    6
#define FPT_THRESHOLD          1000

#define NFSOFT_INDEX_TWO(m,n,l,B) ((B+1)*(B+1)+(B+1)*(B+1)*(m+B)-((m-1)*m*(2*m-1)+(B+1)*(B+2)*(2*B+3))/6)+(posN(n,m,B))+(l-MAX(ABS(m),ABS(n)))

static fpt_set* SO3_fpt_init(int l, unsigned int flags, int kappa, int nthreads);
static int posN(int n, int m, int B);

void nfsoft_init(nfsoft_plan *plan, int N, int M)
{
  nfsoft_init_advanced(plan, N, M, NFSOFT_MALLOC_X | NFSOFT_MALLOC_F
      | NFSOFT_MALLOC_F_HAT);
}

void nfsoft_init_advanced(nfsoft_plan *plan, int N, int M,
    unsigned int nfsoft_flags)
{
  nfsoft_init_guru(plan, N, M, nfsoft_flags, PRE_PHI_HUT | PRE_PSI | MALLOC_X | NFFT_OMP_BLOCKWISE_ADJOINT
      | MALLOC_F_HAT | MALLOC_F | FFTW_INIT,
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

  plan->wig_coeffs = NULL;
  plan->cheby = NULL;
  plan->aux = NULL;

  plan->mv_trafo = (void (*) (void* ))nfsoft_trafo;
  plan->mv_adjoint = (void (*) (void* ))nfsoft_adjoint;

  plan->nthreads = Y(get_num_threads)();

  plan->internal_fpt_set = SO3_fpt_init(plan->N_total, plan->flags, fpt_kappa, plan->nthreads);

}

static void c2e(nfsoft_plan *my_plan, int even, C* wig_coeffs, int k, int m)
{
  int j, N;

  /**initialize the bigger plan*/
  N = 2* (my_plan ->N_total+1);

  /** prepare the coefficients for the new plan */
  C cheby[2*my_plan->N_total + 2];
  cheby[my_plan->N_total+1] = wig_coeffs[0];
  cheby[0]=0.0;

  for (j=1;j<my_plan->N_total+1;j++)
  {
    cheby[my_plan->N_total+1+j]=0.5* wig_coeffs[j];
    cheby[my_plan->N_total+1-j]=0.5* wig_coeffs[j];
  }

  C aux[N+2];

  for(j=1;j<N;j++)
  aux[j]=cheby[j];

  aux[0]=0.;
  aux[N]=0.;

  if (even>0)
  {
    cheby[0]=(C) (-1.)/(2.*_Complex_I) * aux[1];
    for (j=1;j<N;j++)
    {
      cheby[j]=(1./(2.*_Complex_I)*(aux[j+1]-aux[j-1]));
    }

  }

  for (int i = 1; i <= 2* my_plan ->N_total + 2; i++)
  {
    my_plan->p_nfft.f_hat[NFSOFT_INDEX(k, m, i - my_plan->N_total - 1, my_plan->N_total) - 1]
      = cheby[i - 1];
  }

}


static fpt_set* SO3_fpt_init(int l, unsigned int flags, int kappa, int nthreads)
{
  fpt_set *set = (fpt_set*)nfft_malloc(nthreads * sizeof(fpt_set));
  int N, t, k_start, k, m;

  /** Read in transform length. */
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

  /** Initialize DPT. */
  unsigned int fptflags = 0U
        | IF(flags & NFSOFT_USE_DPT,FPT_NO_FAST_ALGORITHM,IF(t > 1,FPT_NO_DIRECT_ALGORITHM,0U))
        | IF(flags & NFSOFT_NO_STABILIZATION,FPT_NO_STABILIZATION,0U);
/*#ifdef _OPENMP
      #pragma omp parallel default(shared) num_threads(nthreads)
      {
        int threadid = omp_get_thread_num();

        if (threadid == 0)
          set[threadid] = fpt_init((2* N + 1) * (2* N + 1), t, fptflags);
        else
          set[threadid] = fpt_init((2* N + 1) * (2* N + 1), t, fptflags | FPT_NO_INIT_FPT_DATA);

        #pragma omp barrier

        if (threadid != 0)
          set[threadid]->dpt = set[0]->dpt;
      }

#else*/
  set[0] = fpt_init((2* N + 1) * (2* N + 1), t, fptflags);
  for (int i=1; i<nthreads; i++)
    {
      set[i] = fpt_init((2* N + 1) * (2* N + 1), t, fptflags | FPT_NO_INIT_FPT_DATA);
      set[i]->dpt = set[0]->dpt;
    }
//#endif

#ifdef _OPENMP
  for (k = -N; k <= N; k++)
    for (m = -N; m <= N; m++)
    {
      /** Read in start and end indices */
      k_start = (ABS(k) >= ABS(m)) ? ABS(k) : ABS(m);

      fpt_precompute_1(set[0], (k+N)*(2*N+1) + m+N,k_start);
    }
  #pragma omp parallel for default(shared) private(k,m,k_start) schedule(dynamic) num_threads(nthreads)
#endif
  for (k = -N; k <= N; k++)
    for (m = -N; m <= N; m++)
    {
      /** Read in start and end indices */
      k_start = (ABS(k) >= ABS(m)) ? ABS(k) : ABS(m);
      // k_end = N;

      R alpha[N+2], beta[N+2], gamma[N+2];
      SO3_alpha_row(alpha, N, k, m);
      SO3_beta_row(beta, N, k, m);
      SO3_gamma_row(gamma, N, k, m);

#ifdef _OPENMP
      fpt_precompute_2(set[omp_get_thread_num()], (k+N)*(2*N+1) + m+N, alpha, beta, gamma, k_start, kappa);
#else
      fpt_precompute(set[0], (k+N)*(2*N+1) + m+N, alpha, beta, gamma, k_start, kappa);
#endif
    }

  return set;
}

static void SO3_fpt(C *coeffs, fpt_set set, int l, int k, int m, unsigned int flags)
{
  int N;
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

  /** Read in start and end indices */
  k_start = (ABS(k) >= ABS(m)) ? ABS(k) : ABS(m);
  k_end = N;
  trafo_nr = (N + k) * (2* N + 1) + (m + N);

  /** Read in Wigner coefficients. */
  C x[k_end + 1];

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

  /** The Chebyshev coefficients. */
  C y[k_end + 1];

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

}

static void SO3_fpt_transposed(C *coeffs, fpt_set set, int l, int k, int m,
    unsigned int flags)
{
  int N, k_start, k_end, j;
  int trafo_nr; /**gives the index of the trafo in the FPT_set*/
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

  /** The Chebychev coefficients. */
  C y[k_end + 1];
  /** The Wigner coefficients. */
  C x[k_end + 1];

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

}

void nfsoft_precompute(nfsoft_plan *plan3D)
{
  int j;
  int M = plan3D->M_total;

  /** Node-dependent part*/

  if (plan3D->x != plan3D->p_nfft.x)
  {/** Resort and copy nodes*/
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
  // glo1 = 0;

  int N = plan3D->N_total;
  int M = plan3D->M_total;

  /**almost nothing to be done for polynomial degree 0*/
  if (N == 0)
  {
    for (int j = 0; j < M; j++)
      plan3D->f[j] = plan3D->f_hat[0];
    return;
  }

  for (int j = 0; j < plan3D->p_nfft.N_total; j++)
    plan3D->p_nfft.f_hat[j] = 0.0;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) num_threads(plan3D->nthreads)
#endif
  for (int k = -N; k <= N; k++)
  {
    C wig_coeffs[(X(next_power_of_2)(N)+1)];
#ifdef _OPENMP
    int threadid = omp_get_thread_num();
#else
    int threadid = 0;
#endif

    for (int m = -N; m <= N; m++)
    {

      int max = (ABS(m) > ABS(k) ? ABS(m) : ABS(k));

      int glo0 = NFSOFT_INDEX_TWO(k,m,max,N);

      for (int j = 0; j <= N - max; j++)
      {
        wig_coeffs[j] = plan3D->f_hat[glo0 + j];

        if ((plan3D->flags & NFSOFT_NORMALIZED))
        {
          wig_coeffs[j] = wig_coeffs[j] * (1. / (2. * KPI))
              * SQRT(0.5 * (2. * (max + j) + 1.));
        }

        if ((plan3D->flags & NFSOFT_REPRESENT))
        {
          if ((k < 0) && (k % 2))
          {
            wig_coeffs[j] = wig_coeffs[j] * (-1);
          }
          if ((m < 0) && (m % 2))
            wig_coeffs[j] = wig_coeffs[j] * (-1);

          if ((m + k) % 2)
	    wig_coeffs[j] = wig_coeffs[j] * (-1);

        }

        // glo1++;
      }

      for (int j = N - max + 1; j < X(next_power_of_2)(N) + 1; j++)
        wig_coeffs[j] = 0.0;
      //fprintf(stdout,"\n k= %d, m= %d \n",k,m);
      SO3_fpt(wig_coeffs, plan3D->internal_fpt_set[threadid], N, k, m, plan3D->flags);

      c2e(plan3D, ABS((k + m) % 2), wig_coeffs, k, m);
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

  if (plan3D->f != plan3D->p_nfft.f)
    for (int j = 0; j < plan3D->M_total; j++)
      plan3D->f[j] = plan3D->p_nfft.f[j];

}

static void e2c(nfsoft_plan *my_plan, int even, C* wig_coeffs, C* cheby)
{
  int N;
  int j;

  /**initialize the bigger plan*/
  N = 2* (my_plan ->N_total+1);
  //nfft_vpr_complex(my_plan->cheby,N+1,"chebychev");

  C aux[N];

  if (even>0)
  {
    //my_plan->aux[N-1]= -1/(2*I)* my_plan->cheby[N-2];
    aux[0]= 1/(2*_Complex_I)*cheby[1];

    for(j=1;j<N-1;j++)
    {
      aux[j]=1/(2*_Complex_I)*(cheby[j+1]-cheby[j-1]);
    }
    aux[N-1]=1/(2*_Complex_I)*(-cheby[j-1]);


    for(j=0;j<N;j++)
    {
    cheby[j]= aux[j];
    }
  }

  wig_coeffs[0]=cheby[my_plan->N_total+1];

  for(j=1;j<=my_plan->N_total;j++)
  {
  wig_coeffs[j]=0.5*(cheby[my_plan->N_total+j+1]+cheby[my_plan->N_total+1-j]);
  }

//nfft_vpr_complex(my_plan->wig_coeffs,my_plan->N_total,"chebychev ");

}

void nfsoft_adjoint(nfsoft_plan *plan3D)
{
  //int glo1 = 0;

  int N = plan3D->N_total;
  int M = plan3D->M_total;

  //nothing much to be done for polynomial degree 0
  if (N == 0)
  {
    plan3D->f_hat[0]=0;
    for (int j = 0; j < M; j++)
      plan3D->f_hat[0] += plan3D->f[j];
    return;
  }

  if (plan3D->p_nfft.f != plan3D->f)
    for (int j = 0; j < M; j++)
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

#ifdef _OPENMP
  #pragma omp parallel for default(shared) num_threads(plan3D->nthreads)
#endif
  for (int k = -N; k <= N; k++)
  {
#ifdef _OPENMP
    int threadid = omp_get_thread_num();
#else
    int threadid = 0;
#endif
    for (int m = -N; m <= N; m++)
    {
      C wig_coeffs[(X(next_power_of_2)(N)+1)];
      C cheby[2*plan3D->N_total + 2];

      int max = (ABS(m) > ABS(k) ? ABS(m) : ABS(k));

      for (int i = 1; i < 2* plan3D ->N_total + 3; i++)
      {
        cheby[i - 1] = plan3D->p_nfft.f_hat[NFSOFT_INDEX(k, m, i - N
            - 1, N) - 1];
      }

      //fprintf(stdout,"k=%d,m=%d \n",k,m);
      //nfft_vpr_complex(plan3D->cheby,2*plan3D->N_total+2,"euler");

      e2c(plan3D, ABS((k + m) % 2), wig_coeffs, cheby);

      //nfft_vpr_complex(plan3D->wig_coeffs,plan3D->N_total+1,"chebys");
      SO3_fpt_transposed(wig_coeffs, plan3D->internal_fpt_set[threadid], N, k, m,
          plan3D->flags);
      //nfft_vpr_complex(plan3D->wig_coeffs,plan3D->N_total+1,"wigners");
      //  SO3_fpt_transposed(plan3D->wig_coeffs,N,k,m,plan3D->flags,plan3D->fpt_kappa);

      int glo0 = NFSOFT_INDEX_TWO(k,m,0,N);

      for (int j = max; j <= N; j++)
      {
        if ((plan3D->flags & NFSOFT_REPRESENT))
        {
          if ((k < 0) && (k % 2))
          {
            wig_coeffs[j] = -wig_coeffs[j];
          }
          if ((m < 0) && (m % 2))
            wig_coeffs[j] = -wig_coeffs[j];

          if ((m + k) % 2)
            wig_coeffs[j] = wig_coeffs[j] * (-1);

        }

        plan3D->f_hat[glo0+j] = wig_coeffs[j];

        if ((plan3D->flags & NFSOFT_NORMALIZED))
        {
          plan3D->f_hat[glo0+j] = plan3D->f_hat[glo0+j] * (1 / (2. * KPI)) * SQRT(
              0.5 * (2. * (j) + 1.));
        }

        //glo1++;
      }

    }
  }
}

void nfsoft_finalize(nfsoft_plan *plan)
{
  /* Finalise the nfft plan. */
  nfft_finalize(&plan->p_nfft);

  for (int i=0; i<plan->nthreads; i++)
    fpt_finalize(plan->internal_fpt_set[i]);
  nfft_free(plan->internal_fpt_set);
  plan->internal_fpt_set = NULL;

  if (plan->flags & NFSOFT_MALLOC_F_HAT)
  {
    //fprintf(stderr,"deallocating f_hat\n");
    nfft_free(plan->f_hat);
  }

  /* De-allocate memory for samples, if necessary. */
  if (plan->flags & NFSOFT_MALLOC_F)
  {
    //fprintf(stderr,"deallocating f\n");
    nfft_free(plan->f);
  }

  /* De-allocate memory for nodes, if necessary. */
  if (plan->flags & NFSOFT_MALLOC_X)
  {
    //fprintf(stderr,"deallocating x\n");
    nfft_free(plan->x);
  }
}

static int posN(int n, int m, int B)
{
  int pos;

  if (n > -B)
    pos = posN(n - 1, m, B) + B + 1 - MAX(ABS(m), ABS(n - 1));
  else
    pos = 0;
  //(n > -B? pos=posN(n-1,m,B)+B+1-MAX(ABS(m),ABS(n-1)): pos= 0)
  return pos;
}

