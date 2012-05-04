/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* Nonequispaced FFT */

/* Authors: D. Potts, S. Kunis 2002-2009, Jens Keiner 2009, Toni Volkmer 2012 */

/* configure header */
#include "config.h"

/* complex datatype (maybe) */
#ifdef HAVE_COMPLEX_H
#include<complex.h>
#endif

/* NFFT headers */
#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define OMP_ASSERT
#ifdef OMP_ASSERT
#endif
#include <assert.h>

/**
 * Sort nodes (index) to get better cache utilization during multiplication
 * with matrix B.
 * The resulting index set is written to ar[2*j+1], the nodes array remains
 * unchanged.
 *
 * \arg n FFTW length (number of oversampled in each dimension)
 * \arg m window length
 * \arg local_x_num number of nodes
 * \arg local_x nodes array
 * \arg ar_x resulting index array
 *
 * \author Toni Volkmer
 */
static void nfft_sort_nodes_for_better_cache_handle(int d,
    const int *n, int m, int local_x_num, const double *local_x, int *ar_x)
{
  int u_j[d], i, j, help, rhigh;
  int *ar_x_temp;
  int nprod;

  for(i = 0; i < local_x_num; i++) {
    ar_x[2*i] = 0;
    ar_x[2*i+1] = i;
    for(j = 0; j < d; j++) {
      help = floor( n[j]*local_x[d*i+j] - m);
      u_j[j] = (help%n[j]+n[j])%n[j];

      ar_x[2*i] += u_j[j];
      if (j+1 < d)
        ar_x[2*i] *= n[j+1];
    }
  }

  for (j = 0, nprod = 1; j < d; j++)
    nprod *= n[j];

  rhigh = ceil(log2(nprod)) - 1;

  ar_x_temp = (int *) nfft_malloc(2*local_x_num*sizeof(int));
  nfft_sort_node_indices_radix_lsdf(local_x_num, ar_x, ar_x_temp, rhigh);
#ifdef OMP_ASSERT
  for (i = 1; i < local_x_num; i++)
    assert(ar_x[2*(i-1)] <= ar_x[2*i]);
#endif
  nfft_free(ar_x_temp);
}

/**
 * Sort nodes (index) to get better cache utilization during multiplication
 * with matrix B.
 * The resulting index set is written to ths->index_x[2*j+1], the nodes array
 * remains unchanged.
 *
 * \arg ths nfft_plan
 */
static void nfft_sort_nodes(const nfft_plan *ths)
{
  if (ths->nfft_flags & NFFT_SORT_NODES)
    nfft_sort_nodes_for_better_cache_handle(ths->d, ths->n, ths->m, ths->M_total, ths->x, ths->index_x);
}

/** direct computation of non equispaced fourier transforms
 *  nfft_trafo_direct, ndft_conjugated, nfft_adjoint_direct, ndft_transposed
 *  require O(M_total N^d) arithemtical operations
 *
 * direct computation of the nfft_trafo_direct and ndft_conjugated, formula (1.1)
 * nfft_trafo_direct:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 * ndft_conjugated:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint_direct and ndft_transposed, formula (1.2)
 * nfft_adjoint_direct:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(+2(pi) k x[j])
 * ndft_transposed:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(-2(pi) k x[j])
 */

/* some macros to initialize arrays before executing a transformation */
#define MACRO_ndft_init_result_trafo memset(f,0,ths->M_total*sizeof(C));
#define MACRO_ndft_init_result_conjugated MACRO_ndft_init_result_trafo
#define MACRO_ndft_init_result_adjoint memset(f_hat,0,ths->N_total*sizeof(C));
#define MACRO_ndft_init_result_transposed MACRO_ndft_init_result_adjoint

/* exponent of complex exponentials */
#define MACRO_ndft_sign_trafo K2PI*ths->x[j*ths->d+t]
#define MACRO_ndft_sign_conjugated -K2PI*ths->x[j*ths->d+t]
#define MACRO_ndft_sign_adjoint K2PI*ths->x[j*ths->d+t]
#define MACRO_ndft_sign_transposed -K2PI*ths->x[j*ths->d+t]

#define MACRO_init_k_N_Omega_x(which_one)                                     \
{                                                                             \
  for (t = 0; t < ths->d; t++)                                                \
  {                                                                           \
    k[t] = -ths->N[t]/2;                                                      \
    x[t] = MACRO_ndft_sign_ ## which_one;                                     \
    Omega[t+1] = k[t]*x[t] + Omega[t];                                        \
  }                                                                           \
  omega = Omega[ths->d];                                                      \
}                                                                             \

#define MACRO_count_k_N_Omega                                                 \
{                                                                             \
  for (t = ths->d-1; (t >= 1) && (k[t] == ths->N[t]/2-1); t--)                \
    k[t]-= ths->N[t]-1;                                                       \
                                                                              \
  k[t]++;                                                                     \
\
  for (t2 = t; t2 < ths->d; t2++)                                             \
    Omega[t2+1] = k[t2]*x[t2] + Omega[t2];                                    \
                                                                              \
  omega = Omega[ths->d];                                                      \
}

#define MACRO_ndft_compute_trafo f[j] += f_hat[k_L]*CEXP(-II*omega);
#define MACRO_ndft_compute_conjugated MACRO_ndft_compute_trafo
#define MACRO_ndft_compute_adjoint f_hat[k_L] += f[j]*CEXP(+II*omega);
#define MACRO_ndft_compute_transposed MACRO_ndft_compute_adjoint

#define MACRO_ndft(which_one)                                                 \
void nfft_ ## which_one ## _direct (nfft_plan *ths)                           \
{                                                                             \
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;                                 \
                                                                              \
  MACRO_ndft_init_result_ ## which_one                                        \
                                                                              \
  if (ths->d == 1)                                                            \
  {                                                                           \
    /* specialize for univariate case, rationale: faster */                   \
    const int t = 0;                                                          \
    int j;                                                                    \
    for (j = 0; j < ths->M_total; j++)                                        \
    {                                                                         \
      int k_L;                                                                \
      for(k_L = 0; k_L < ths->N_total; k_L++)                                 \
      {                                                                       \
        R omega = (k_L - (ths->N_total/2)) * MACRO_ndft_sign_ ## which_one;   \
        MACRO_ndft_compute_ ## which_one;                                     \
      }                                                                       \
    }                                                                         \
  }                                                                           \
  else                                                                        \
  {                                                                           \
    /* multivariate case */                                                   \
    int j;                                                                    \
    for (j = 0; j < ths->M_total; j++)                                        \
    {                                                                         \
      R x[ths->d], omega, Omega[ths->d+1];                                    \
      int t, t2, k_L, k[ths->d];                                              \
      Omega[0] = K(0.0);                                                      \
      MACRO_init_k_N_Omega_x(which_one);                                      \
      for(k_L = 0; k_L < ths->N_total; k_L++)                                 \
      {                                                                       \
        MACRO_ndft_compute_ ## which_one;                                     \
        MACRO_count_k_N_Omega;                                                \
      }                                                                       \
    }                                                                         \
  }                                                                           \
}

//MACRO_ndft(trafo)
void nfft_trafo_direct (nfft_plan *ths)
{
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;

  memset(f,0,ths->M_total*sizeof(C));

  if (ths->d == 1)
  {
    /* specialize for univariate case, rationale: faster */
    int j;
    #pragma omp parallel for default(shared) private(j)
    for (j = 0; j < ths->M_total; j++)
    {
      int k_L;
      for(k_L = 0; k_L < ths->N_total; k_L++)
      {
        R omega = (k_L - ths->N_total/2) * K2PI * ths->x[j];
        f[j] += f_hat[k_L]*cexp(-( 1.0iF)*omega);
      }
    }
  }
  else
  {
    /* multivariate case */
    int j;
    #pragma omp parallel for default(shared) private(j)
    for (j = 0; j < ths->M_total; j++)
    {
      R x[ths->d], omega, Omega[ths->d+1];
      int t, t2, k_L, k[ths->d];
      Omega[0] = ((R) 0.0);
      for (t = 0; t < ths->d; t++)
      {
        k[t] = -ths->N[t]/2;
        x[t] = K2PI*ths->x[j*ths->d+t];
        Omega[t+1] = k[t]*x[t] + Omega[t];
      }
      omega = Omega[ths->d];

      for(k_L = 0; k_L < ths->N_total; k_L++)
      {
        f[j] += f_hat[k_L]*cexp(-( 1.0iF)*omega);
        {
          for (t = ths->d-1; (t >= 1) && (k[t] == ths->N[t]/2-1); t--)
            k[t]-= ths->N[t]-1;

          k[t]++;

          for (t2 = t; t2 < ths->d; t2++)
            Omega[t2+1] = k[t2]*x[t2] + Omega[t2];

          omega = Omega[ths->d];
        }
      }
    }
  }
}

//MACRO_ndft(adjoint)
void nfft_adjoint_direct (nfft_plan *ths)
{
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;

  memset(f_hat,0,ths->N_total*sizeof(C));

  if (ths->d == 1)
  {
    /* specialize for univariate case, rationale: faster */
#ifdef _OPENMP
      int k_L;
      #pragma omp parallel for default(shared) private(k_L)
      for(k_L = 0; k_L < ths->N_total; k_L++)
      {
        int j;
        for (j = 0; j < ths->M_total; j++)
        {
          R omega = (k_L - (ths->N_total/2)) * K2PI * ths->x[j];
          f_hat[k_L] += f[j]*cexp(+( 1.0iF)*omega);
        }
      }
#else
      int j;
      for (j = 0; j < ths->M_total; j++)
      {
        int k_L;
        for(k_L = 0; k_L < ths->N_total; k_L++)
        {
          R omega = (k_L - (ths->N_total/2)) * K2PI * ths->x[j];
          f_hat[k_L] += f[j]*cexp(+( 1.0iF)*omega);
        }
      }
#endif
  }
  else
  {
    /* multivariate case */
    int j, k_L;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j, k_L)
    for(k_L = 0; k_L < ths->N_total; k_L++)
    {
      int k[ths->d], k_temp, t;

      k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        k[t] = k_temp % ths->N[t] - ths->N[t]/2;
        k_temp /= ths->N[t];
      }

      for (j = 0; j < ths->M_total; j++)
      {
        R omega = 0.0;
        for (t = 0; t < ths->d; t++)
          omega += k[t] * K2PI * ths->x[j*ths->d+t];
        f_hat[k_L] += f[j]*cexp(+( 1.0iF)*omega);
      }
    }
#else
    for (j = 0; j < ths->M_total; j++)
    {
      R x[ths->d], omega, Omega[ths->d+1];
      int t, t2, k[ths->d];
      Omega[0] = ((R) 0.0);
      for (t = 0; t < ths->d; t++)
      {
        k[t] = -ths->N[t]/2;
        x[t] = K2PI * ths->x[j*ths->d+t];
        Omega[t+1] = k[t]*x[t] + Omega[t];
      }
      omega = Omega[ths->d];
      for(k_L = 0; k_L < ths->N_total; k_L++)
      {
        f_hat[k_L] += f[j]*cexp(+( 1.0iF)*omega);

        for (t = ths->d-1; (t >= 1) && (k[t] == ths->N[t]/2-1); t--)
          k[t]-= ths->N[t]-1;

        k[t]++;

        for (t2 = t; t2 < ths->d; t2++)
          Omega[t2+1] = k[t2]*x[t2] + Omega[t2];

        omega = Omega[ths->d];
      }
    }
#endif
  }
}

/** fast computation of non equispaced fourier transforms
 *  require O(N^d log(N) + M_total) arithmetical operations
 *
 * fast computation of the nfft_trafo and nfft_conjugated, formula (1.1)
 * nfft_trafo:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 * nfft_conjugated:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint and nfft_transposed, formula (1.2)
 * nfft_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(+2(pi) k x[j])
 * nfft_transposed:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(-2(pi) k x[j])
 */

/** macros and small sub routines for the fast transforms
 */

/** computes 2m+2 indices for the matrix B
 */
static inline void nfft_uo(const nfft_plan *ths, const int j, int *up, int *op,
  const int act_dim)
{
  const R xj = ths->x[j * ths->d + act_dim];
  int c = LRINT(FLOOR(xj * ths->n[act_dim]));

  (*up) = c - (ths->m);
  (*op) = c + 1 + (ths->m);
}

static inline void nfft_uo2(int *u, int *o, const R x, const int n, const int m)
{
  int c = LRINT(FLOOR(x * n));

  *u = (c - m + n) % n;
  *o = (c + 1 + m + n) % n;
}

#define MACRO_nfft_D_compute_A                                                \
{                                                                             \
  g_hat[k_plain[ths->d]] = f_hat[ks_plain[ths->d]] * c_phi_inv_k[ths->d];     \
}

#define MACRO_nfft_D_compute_T                                                \
{                                                                             \
  f_hat[ks_plain[ths->d]] = g_hat[k_plain[ths->d]] * c_phi_inv_k[ths->d];     \
}

#define MACRO_nfft_D_init_result_A memset(g_hat,0,ths->n_total*sizeof(C));

#define MACRO_nfft_D_init_result_T memset(f_hat,0,ths->N_total*sizeof(C));

#define MACRO_with_PRE_PHI_HUT * ths->c_phi_inv[t2][ks[t2]];

#define MACRO_without_PRE_PHI_HUT / (PHI_HUT(ks[t2]-(ths->N[t2]/2),t2));

#define MACRO_init_k_ks                                                       \
{                                                                             \
  for (t = ths->d-1; 0 <= t; t--)                                             \
  {                                                                           \
    kp[t] = k[t] = 0;                                                         \
    ks[t] = ths->N[t]/2;                                                      \
  }                                                                           \
  t++;                                                                        \
}

#define MACRO_update_c_phi_inv_k(which_one) \
{                                                                             \
  for (t2 = t; t2 < ths->d; t2++)                                             \
  {                                                                           \
    c_phi_inv_k[t2+1] = c_phi_inv_k[t2] MACRO_ ##which_one;                   \
    ks_plain[t2+1] = ks_plain[t2]*ths->N[t2] + ks[t2];                        \
    k_plain[t2+1] = k_plain[t2]*ths->n[t2] + k[t2];                           \
  }                                                                           \
}

#define MACRO_count_k_ks                                                      \
{                                                                             \
  for (t = ths->d-1; (t > 0) && (kp[t] == ths->N[t]-1); t--)                  \
  {                                                                           \
    kp[t] = k[t] = 0;                                                         \
    ks[t]= ths->N[t]/2;                                                       \
  }                                                                           \
                                                                              \
  kp[t]++; k[t]++; ks[t]++;                                                   \
  if(kp[t] == ths->N[t]/2)                                                    \
  {                                                                           \
    k[t] = ths->n[t] - ths->N[t]/2;                                           \
    ks[t] = 0;                                                                \
  }                                                                           \
}                                                                             \


/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$D, D^T\f$
 */
#define MACRO_nfft_D(which_one)                                               \
static inline void nfft_D_serial_ ## which_one (nfft_plan *ths)               \
{                                                                             \
  C *f_hat, *g_hat;                     /**< local copy                     */\
  R c_phi_inv_k[ths->d+1];              /**< postfix product of PHI_HUT     */\
  int t, t2;                            /**< index dimensions               */\
  int k_L;                              /**< plain index                    */\
  int kp[ths->d];                       /**< multi index (simple)           */\
  int k[ths->d];                        /**< multi index in g_hat           */\
  int ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/\
  int k_plain[ths->d+1];                /**< postfix plain index            */\
  int ks_plain[ths->d+1];               /**< postfix plain index            */\
                                                                              \
  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat;                             \
  MACRO_nfft_D_init_result_ ## which_one;                                     \
                                                                              \
  c_phi_inv_k[0] = K(1.0);                                                    \
  k_plain[0] = 0;                                                             \
  ks_plain[0] = 0;                                                            \
                                                                              \
  if (ths->nfft_flags & PRE_PHI_HUT)                                          \
  {                                                                           \
    MACRO_init_k_ks;                                                          \
                                                                              \
    for (k_L = 0; k_L < ths->N_total; k_L++)                                  \
    {                                                                         \
      MACRO_update_c_phi_inv_k(with_PRE_PHI_HUT);                             \
      MACRO_nfft_D_compute_ ## which_one;                                     \
      MACRO_count_k_ks;                                                       \
    } /* for(k_L) */                                                          \
  } /* if(PRE_PHI_HUT) */                                                     \
  else                                                                        \
  {                                                                           \
    MACRO_init_k_ks;                                                          \
    for (k_L = 0; k_L < ths->N_total; k_L++)                                  \
    {                                                                         \
      MACRO_update_c_phi_inv_k(without_PRE_PHI_HUT);                          \
      MACRO_nfft_D_compute_ ## which_one;                                     \
      MACRO_count_k_ks;                                                       \
    } /* for(k_L) */                                                          \
  } /* else(PRE_PHI_HUT) */                                                   \
} /* nfft_D */

static void nfft_D_openmp_A(nfft_plan *ths)
{
  C *f_hat, *g_hat;                     /**< local copy                     */
  int k_L;                              /**< plain index                    */

  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat;
  memset(g_hat,0,ths->n_total*sizeof(C));

  if (ths->nfft_flags & PRE_PHI_HUT)
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      int kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      int k[ths->d];                        /**< multi index in g_hat           */
      int ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      double c_phi_inv_k_val = K(1.0);
      int k_plain_val = 0;
      int ks_plain_val = 0;
      int t;
      int k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val *= ths->c_phi_inv[t][ks[t]];
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*ths->n[t] + k[t];
      }

      g_hat[k_plain_val] = f_hat[ks_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* if(PRE_PHI_HUT) */
  else
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      int kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      int k[ths->d];                        /**< multi index in g_hat           */
      int ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      double c_phi_inv_k_val = K(1.0);
      int k_plain_val = 0;
      int ks_plain_val = 0;
      int t;
      int k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val /= (PHI_HUT(ks[t]-(ths->N[t]/2),t));
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*ths->n[t] + k[t];
      }

      g_hat[k_plain_val] = f_hat[ks_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* else(PRE_PHI_HUT) */
}
MACRO_nfft_D(A)
static void nfft_D_A(nfft_plan *ths)
{
#ifdef _OPENMP
  nfft_D_openmp_A(ths);
#else
  nfft_D_serial_A(ths);
#endif
}

static void nfft_D_openmp_T(nfft_plan *ths)
{
  C *f_hat, *g_hat;                     /**< local copy                     */
  int k_L;                              /**< plain index                    */

  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat;
  memset(f_hat,0,ths->N_total*sizeof(C));

  if (ths->nfft_flags & PRE_PHI_HUT)
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      int kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      int k[ths->d];                        /**< multi index in g_hat           */
      int ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      double c_phi_inv_k_val = K(1.0);
      int k_plain_val = 0;
      int ks_plain_val = 0;
      int t;
      int k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val *= ths->c_phi_inv[t][ks[t]];
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*ths->n[t] + k[t];
      }

      f_hat[ks_plain_val] = g_hat[k_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* if(PRE_PHI_HUT) */
  else
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      int kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      int k[ths->d];                        /**< multi index in g_hat           */
      int ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      double c_phi_inv_k_val = K(1.0);
      int k_plain_val = 0;
      int ks_plain_val = 0;
      int t;
      int k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val /= (PHI_HUT(ks[t]-(ths->N[t]/2),t));
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*ths->n[t] + k[t];
      }

      f_hat[ks_plain_val] = g_hat[k_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* else(PRE_PHI_HUT) */
}

MACRO_nfft_D(T)

static void nfft_D_T(nfft_plan *ths)
{
#ifdef _OPENMP
  nfft_D_openmp_T(ths);
#else
  nfft_D_serial_T(ths);
#endif
}

/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$B, B^{\rm T}\f$
 */
#define MACRO_nfft_B_init_result_A  memset(f,0,ths->M_total*sizeof(C));
#define MACRO_nfft_B_init_result_T memset(g,0,ths->n_total*sizeof(C));

#define MACRO_nfft_B_PRE_FULL_PSI_compute_A                                   \
{                                                                             \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]];                            \
}

#define MACRO_nfft_B_PRE_FULL_PSI_compute_T                                   \
{                                                                             \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj);                            \
}

#define MACRO_nfft_B_compute_A                                                \
{                                                                             \
  (*fj) += phi_prod[ths->d] * g[ll_plain[ths->d]];                            \
}

#define MACRO_nfft_B_compute_T                                                \
{                                                                             \
  g[ll_plain[ths->d]] += phi_prod[ths->d] * (*fj);                            \
}

#define MACRO_with_FG_PSI fg_psi[t2][lj[t2]]

#define MACRO_with_PRE_PSI ths->psi[(j*ths->d+t2) * (2*ths->m+2)+lj[t2]]

#define MACRO_without_PRE_PSI  PHI(ths->x[j*ths->d+t2]                        \
  - ((R)l[t2])/((R)ths->n[t2]), t2)

#define MACRO_init_uo_l_lj_t                                                  \
{                                                                             \
  for (t = ths->d-1; t >= 0; t--)                                             \
  {                                                                           \
    nfft_uo(ths,j,&u[t],&o[t],t);                                             \
    l[t] = u[t];                                                              \
    lj[t] = 0;                                                                \
  } /* for(t) */                                                              \
  t++;                                                                        \
}

#define MACRO_update_phi_prod_ll_plain(which_one) {                           \
  for(t2=t; t2<ths->d; t2++)                                                  \
    {                                                                         \
      phi_prod[t2+1] = phi_prod[t2]* MACRO_ ## which_one;                     \
      ll_plain[t2+1] = ll_plain[t2]*ths->n[t2] +(l[t2]+ths->n[t2])%ths->n[t2];\
    } /* for(t2) */                                                           \
}

#define MACRO_count_uo_l_lj_t {                                               \
  for(t = ths->d-1; (t > 0) && (l[t] == o[t]); t--)                           \
    {                                                                         \
      l[t] = u[t];                                                            \
      lj[t] = 0;                                                              \
    } /* for(t) */                                                            \
                                                                              \
  l[t]++;                                                                     \
  lj[t]++;                                                                    \
}

#define MACRO_nfft_B(which_one)                                               \
static inline void nfft_B_serial_ ## which_one (nfft_plan *ths)               \
{                                                                             \
  int lprod; /* 'regular bandwidth' of matrix B  */                           \
  int u[ths->d], o[ths->d]; /* multi band with respect to x_j */              \
  int t, t2; /* index dimensions */                                           \
  int j; /* index nodes */                                                    \
  int l_L, ix; /* index one row of B */                                       \
  int l[ths->d]; /* multi index u<=l<=o */                                    \
  int lj[ths->d]; /* multi index 0<=lj<u+o+1 */                               \
  int ll_plain[ths->d+1]; /* postfix plain index in g */                      \
  R phi_prod[ths->d+1]; /* postfix product of PHI */                          \
  C *f, *g; /* local copy */                                                  \
  C *fj; /* local copy */                                                     \
  R y[ths->d];                                                                \
  R fg_psi[ths->d][2*ths->m+2];                                               \
  R fg_exp_l[ths->d][2*ths->m+2];                                             \
  int l_fg,lj_fg;                                                             \
  R tmpEXP1, tmpEXP2, tmpEXP2sq, tmp1, tmp2, tmp3;                            \
  R ip_w;                                                                     \
  int ip_u;                                                                   \
  int ip_s = ths->K/(ths->m+2);                                               \
                                                                              \
  f = (C*)ths->f; g = (C*)ths->g;                                             \
                                                                              \
  MACRO_nfft_B_init_result_ ## which_one;                                     \
                                                                              \
  if (ths->nfft_flags & PRE_FULL_PSI)                                         \
  {                                                                           \
    for (ix = 0, j = 0, fj = f; j < ths->M_total; j++, fj++)                  \
      for (l_L = 0; l_L < ths->psi_index_f[j]; l_L++, ix++)                   \
        MACRO_nfft_B_PRE_FULL_PSI_compute_ ## which_one;                      \
    return;                                                                   \
  }                                                                           \
                                                                              \
  phi_prod[0] = K(1.0);                                                       \
  ll_plain[0] = 0;                                                            \
                                                                              \
  for (t = 0, lprod = 1; t < ths->d; t++)                                     \
    lprod *= (2*ths->m+2);                                                    \
                                                                              \
  if (ths->nfft_flags & PRE_PSI)                                              \
  {                                                                           \
    for (j = 0, fj = f; j < ths->M_total; j++, fj++)                          \
    {                                                                         \
      MACRO_init_uo_l_lj_t;                                                   \
                                                                              \
      for (l_L = 0; l_L < lprod; l_L++)                                       \
      {                                                                       \
        MACRO_update_phi_prod_ll_plain(with_PRE_PSI);                         \
                                                                              \
        MACRO_nfft_B_compute_ ## which_one;                                   \
                                                                              \
        MACRO_count_uo_l_lj_t;                                                \
      } /* for(l_L) */                                                        \
    } /* for(j) */                                                            \
    return;                                                                   \
  } /* if(PRE_PSI) */                                                         \
                                                                              \
  if (ths->nfft_flags & PRE_FG_PSI)                                           \
  {                                                                           \
    for(t2 = 0; t2 < ths->d; t2++)                                            \
    {                                                                         \
      tmpEXP2 = EXP(K(-1.0)/ths->b[t2]);                                      \
      tmpEXP2sq = tmpEXP2*tmpEXP2;                                            \
      tmp2 = K(1.0);                                                          \
      tmp3 = K(1.0);                                                          \
      fg_exp_l[t2][0] = K(1.0);                                               \
      for(lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++)                          \
      {                                                                       \
        tmp3 = tmp2*tmpEXP2;                                                  \
        tmp2 *= tmpEXP2sq;                                                    \
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;                     \
      }                                                                       \
    }                                                                         \
    for (j = 0, fj = f; j < ths->M_total; j++, fj++)                          \
    {                                                                         \
      MACRO_init_uo_l_lj_t;                                                   \
                                                                              \
      for (t2 = 0; t2 < ths->d; t2++)                                         \
      {                                                                       \
        fg_psi[t2][0] = ths->psi[2*(j*ths->d+t2)];                            \
        tmpEXP1 = ths->psi[2*(j*ths->d+t2)+1];                                \
        tmp1 = K(1.0);                                                        \
        for (l_fg = u[t2]+1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++)       \
        {                                                                     \
          tmp1 *= tmpEXP1;                                                    \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg];         \
        }                                                                     \
      }                                                                       \
                                                                              \
      for (l_L= 0; l_L < lprod; l_L++)                                        \
      {                                                                       \
        MACRO_update_phi_prod_ll_plain(with_FG_PSI);                          \
                                                                              \
        MACRO_nfft_B_compute_ ## which_one;                                   \
                                                                              \
        MACRO_count_uo_l_lj_t;                                                \
      } /* for(l_L) */                                                        \
    } /* for(j) */                                                            \
    return;                                                                   \
  } /* if(PRE_FG_PSI) */                                                      \
                                                                              \
  if (ths->nfft_flags & FG_PSI)                                               \
  {                                                                           \
    for (t2 = 0; t2 < ths->d; t2++)                                           \
    {                                                                         \
      tmpEXP2 = EXP(K(-1.0)/ths->b[t2]);                                      \
      tmpEXP2sq = tmpEXP2*tmpEXP2;                                            \
      tmp2 = K(1.0);                                                          \
      tmp3 = K(1.0);                                                          \
      fg_exp_l[t2][0] = K(1.0);                                               \
      for (lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++)                         \
      {                                                                       \
        tmp3 = tmp2*tmpEXP2;                                                  \
        tmp2 *= tmpEXP2sq;                                                    \
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;                     \
      }                                                                       \
    }                                                                         \
    for (j = 0, fj = f; j < ths->M_total; j++, fj++)                          \
    {                                                                         \
      MACRO_init_uo_l_lj_t;                                                   \
                                                                              \
      for (t2 = 0; t2 < ths->d; t2++)                                         \
      {                                                                       \
        fg_psi[t2][0] = (PHI((ths->x[j*ths->d+t2]-((R)u[t2])/ths->n[t2]),t2));\
                                                                              \
        tmpEXP1 = EXP(K(2.0)*(ths->n[t2]*ths->x[j*ths->d+t2] - u[t2])         \
          /ths->b[t2]);                                                       \
        tmp1 = K(1.0);                                                        \
        for (l_fg = u[t2] + 1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++)     \
        {                                                                     \
          tmp1 *= tmpEXP1;                                                    \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg];         \
        }                                                                     \
      }                                                                       \
                                                                              \
      for (l_L = 0; l_L < lprod; l_L++)                                       \
      {                                                                       \
        MACRO_update_phi_prod_ll_plain(with_FG_PSI);                          \
                                                                              \
        MACRO_nfft_B_compute_ ## which_one;                                   \
                                                                              \
        MACRO_count_uo_l_lj_t;                                                \
      } /* for(l_L) */                                                        \
    } /* for(j) */                                                            \
    return;                                                                   \
  } /* if(FG_PSI) */                                                          \
                                                                              \
  if (ths->nfft_flags & PRE_LIN_PSI)                                          \
  {                                                                           \
    for (j = 0, fj=f; j<ths->M_total; j++, fj++)                              \
    {                                                                         \
      MACRO_init_uo_l_lj_t;                                                   \
                                                                              \
      for (t2 = 0; t2 < ths->d; t2++)                                         \
      {                                                                       \
        y[t2] = ((ths->n[t2]*ths->x[j*ths->d+t2]-(R)u[t2])                    \
          * ((R)ths->K))/(ths->m+2);                                          \
        ip_u  = LRINT(FLOOR(y[t2]));                                          \
        ip_w  = y[t2]-ip_u;                                                   \
        for (l_fg = u[t2], lj_fg = 0; l_fg <= o[t2]; l_fg++, lj_fg++)         \
        {                                                                     \
          fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s)]  \
            * (1-ip_w) + ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s+1)]     \
            * (ip_w);                                                         \
        }                                                                     \
      }                                                                       \
                                                                              \
      for (l_L = 0; l_L < lprod; l_L++)                                       \
      {                                                                       \
        MACRO_update_phi_prod_ll_plain(with_FG_PSI);                          \
                                                                              \
        MACRO_nfft_B_compute_ ## which_one;                                   \
                                                                              \
        MACRO_count_uo_l_lj_t;                                                \
      } /* for(l_L) */                                                        \
    } /* for(j) */                                                            \
    return;                                                                   \
  } /* if(PRE_LIN_PSI) */                                                     \
                                                                              \
  /* no precomputed psi at all */                                             \
  for (j = 0, fj = f; j < ths->M_total; j++, fj++)                            \
  {                                                                           \
    MACRO_init_uo_l_lj_t;                                                     \
                                                                              \
    for (l_L = 0; l_L < lprod; l_L++)                                         \
    {                                                                         \
      MACRO_update_phi_prod_ll_plain(without_PRE_PSI);                        \
                                                                              \
      MACRO_nfft_B_compute_ ## which_one;                                     \
                                                                              \
      MACRO_count_uo_l_lj_t;                                                  \
    } /* for(l_L) */                                                          \
  } /* for(j) */                                                              \
} /* nfft_B */                                                                \

MACRO_nfft_B(A)

static inline void nfft_B_openmp_A (nfft_plan *ths)
{
  int lprod; /* 'regular bandwidth' of matrix B  */
  int k;

  memset(ths->f,0,ths->M_total*sizeof(C));

  for (k = 0, lprod = 1; k < ths->d; k++)
    lprod *= (2*ths->m+2);

  if (ths->nfft_flags & PRE_FULL_PSI)
  {
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < ths->M_total; k++)
    {
      int l;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < lprod; l++)
        ths->f[j] += ths->psi[j*lprod+l] * ths->g[ths->psi_index_g[j*lprod+l]];
    }
    return;
  }

  if (ths->nfft_flags & PRE_PSI)
  {
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < ths->M_total; k++)
    {
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int t, t2; /* index dimensions */
      int l_L; /* index one row of B */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (l_L = 0; l_L < lprod; l_L++)
      {
        MACRO_update_phi_prod_ll_plain(with_PRE_PSI);

        ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]];

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(PRE_PSI) */

  if (ths->nfft_flags & PRE_FG_PSI)
  {
    int t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];

    for(t2 = 0; t2 < ths->d; t2++)
    {
      int lj_fg;
      R tmpEXP2 = EXP(K(-1.0)/ths->b[t2]);
      R tmpEXP2sq = tmpEXP2*tmpEXP2;
      R tmp2 = K(1.0);
      R tmp3 = K(1.0);
      fg_exp_l[t2][0] = K(1.0);
      for(lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++)
      {
        tmp3 = tmp2*tmpEXP2;
        tmp2 *= tmpEXP2sq;
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;
      }
    }

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      R fg_psi[ths->d][2*ths->m+2];
      R tmpEXP1, tmp1;
      int l_fg,lj_fg;
      int l_L;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (t2 = 0; t2 < ths->d; t2++)
      {
        fg_psi[t2][0] = ths->psi[2*(j*ths->d+t2)];
        tmpEXP1 = ths->psi[2*(j*ths->d+t2)+1];
        tmp1 = K(1.0);
        for (l_fg = u[t2]+1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++)
        {
          tmp1 *= tmpEXP1;
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg];
        }
      }

      for (l_L= 0; l_L < lprod; l_L++)
      {
        MACRO_update_phi_prod_ll_plain(with_FG_PSI);

        ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]];

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(PRE_FG_PSI) */

  if (ths->nfft_flags & FG_PSI)
  {
    int t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];

    nfft_sort_nodes(ths);

    for (t2 = 0; t2 < ths->d; t2++)
    {
      int lj_fg;
      R tmpEXP2 = EXP(K(-1.0)/ths->b[t2]);
      R tmpEXP2sq = tmpEXP2*tmpEXP2;
      R tmp2 = K(1.0);
      R tmp3 = K(1.0);
      fg_exp_l[t2][0] = K(1.0);
      for (lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++)
      {
        tmp3 = tmp2*tmpEXP2;
        tmp2 *= tmpEXP2sq;
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;
      }
    }

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      R fg_psi[ths->d][2*ths->m+2];
      R tmpEXP1, tmp1;
      int l_fg,lj_fg;
      int l_L;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (t2 = 0; t2 < ths->d; t2++)
      {
        fg_psi[t2][0] = (PHI((ths->x[j*ths->d+t2]-((R)u[t2])/ths->n[t2]),t2));

        tmpEXP1 = EXP(K(2.0)*(ths->n[t2]*ths->x[j*ths->d+t2] - u[t2])
          /ths->b[t2]);
        tmp1 = K(1.0);
        for (l_fg = u[t2] + 1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++)
        {
          tmp1 *= tmpEXP1;
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg];
        }
      }

      for (l_L = 0; l_L < lprod; l_L++)
      {
        MACRO_update_phi_prod_ll_plain(with_FG_PSI);

        ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]];

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(FG_PSI) */

  if (ths->nfft_flags & PRE_LIN_PSI)
  {
    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k<ths->M_total; k++)
    {
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int t, t2; /* index dimensions */
      int l_L; /* index one row of B */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      R y[ths->d];
      R fg_psi[ths->d][2*ths->m+2];
      int l_fg,lj_fg;
      R ip_w;
      int ip_u;
      int ip_s = ths->K/(ths->m+2);
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (t2 = 0; t2 < ths->d; t2++)
      {
        y[t2] = ((ths->n[t2]*ths->x[j*ths->d+t2]-(R)u[t2])
          * ((R)ths->K))/(ths->m+2);
        ip_u  = LRINT(FLOOR(y[t2]));
        ip_w  = y[t2]-ip_u;
        for (l_fg = u[t2], lj_fg = 0; l_fg <= o[t2]; l_fg++, lj_fg++)
        {
          fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s)]
            * (1-ip_w) + ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s+1)]
            * (ip_w);
        }
      }

      for (l_L = 0; l_L < lprod; l_L++)
      {
        MACRO_update_phi_prod_ll_plain(with_FG_PSI);

        ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]];

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  nfft_sort_nodes(ths);

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < ths->M_total; k++)
  {
    int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
    int t, t2; /* index dimensions */
    int l_L; /* index one row of B */
    int l[ths->d]; /* multi index u<=l<=o */
    int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
    int ll_plain[ths->d+1]; /* postfix plain index in g */
    R phi_prod[ths->d+1]; /* postfix product of PHI */
    int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    phi_prod[0] = K(1.0);
    ll_plain[0] = 0;

    MACRO_init_uo_l_lj_t;

    for (l_L = 0; l_L < lprod; l_L++)
    {
      MACRO_update_phi_prod_ll_plain(without_PRE_PSI);

      ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]];

      MACRO_count_uo_l_lj_t;
    } /* for(l_L) */
  } /* for(j) */
}

static void nfft_B_A(nfft_plan *ths)
{
#ifdef _OPENMP
  nfft_B_openmp_A(ths);
#else
  nfft_B_serial_A(ths);
#endif
}

#ifdef _OPENMP
/**
 * Performs binary search in sorted index array and returns the offset of the
 * specified key.
 * If the key is not unique, the index of the left-most element with the
 * requested key is returned.
 * If the index array does not contain the key, the index of the next-largest
 * element is returned.
 *
 * \arg ar_x sorted index array containing the key at offset 2*k
 * and the nodes index at offset 2*k+1
 * \arg len number of nodes x
 * \arg key the key value
 *
 * \author Toni Volkmer
 */
static inline int index_x_binary_search(const int *ar_x, const int len, const int key)
{
  int left = 0, right = len - 1;

  if (len == 1)
    return 0;

  while (left < right - 1)
  {
    int i = (left + right) / 2;
    if (ar_x[2*i] >= key)
      right = i;
    else if (ar_x[2*i] < key)
      left = i;
  }

  if (ar_x[2*left] < key && left != len-1)
    return left+1;

  return left;
}
#endif

#ifdef _OPENMP
/**
 * Determines the blocks of vector g the current thread is responsible for.
 *
 * \arg my_u0 lowest index (first component) the current threads writes to in g
 * \arg my_o0 highest index (first component) the current threads writes to in g
 * \arg min_u_a lowest (linearized) index u which could lead to writing to g
 * \arg max_u_a highest (linearized) index u which could lead to writing to g
 * \arg min_u_b lowest (linearized) index u which could lead to writing to g
 * \arg max_u_b highest (linearized) index u which could lead to writing to g
 * \arg d dimensionality
 * \arg n FFTW length
 * \arg m window length
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_B_omp3_init(int *my_u0, int *my_o0, int *min_u_a, int *max_u_a, int *min_u_b, int *max_u_b, const int d, const int *n, const int m)
{
  const int n0 = n[0];
  int k;
  int nthreads = omp_get_num_threads();
  int nthreads_used = MIN(nthreads, n0);
  int size_per_thread = n0 / nthreads_used;
  int size_left = n0 - size_per_thread * nthreads_used;
  int size_g[nthreads_used];
  int offset_g[nthreads_used];
  int my_id = omp_get_thread_num();
  int n_prod_rest = 1;

  for (k = 1; k < d; k++)
    n_prod_rest *= n[k];

  *min_u_a = -1;
  *max_u_a = -1;
  *min_u_b = -1;
  *max_u_b = -1;
  *my_u0 = -1;
  *my_o0 = -1;

  if (my_id < nthreads_used)
  {
    const int m22 = 2 * m + 2;

    offset_g[0] = 0;
    for (k = 0; k < nthreads_used; k++)
    {
      if (k > 0)
        offset_g[k] = offset_g[k-1] + size_g[k-1];
      size_g[k] = size_per_thread;
      if (size_left > 0)
      {
        size_g[k]++;
        size_left--;
      }
    }

    *my_u0 = offset_g[my_id];
    *my_o0 = offset_g[my_id] + size_g[my_id] - 1;

    if (nthreads_used > 1)
    {
      *max_u_a = n_prod_rest*(offset_g[my_id] + size_g[my_id]) - 1;
      *min_u_a = n_prod_rest*(offset_g[my_id] - m22 + 1);
    }
    else
    {
      *min_u_a = 0;
      *max_u_a = n_prod_rest * n0 - 1;
    }

    if (*min_u_a < 0)
    {
      *min_u_b = n_prod_rest * (offset_g[my_id] - m22 + 1 + n0);
      *max_u_b = n_prod_rest * n0 - 1;
      *min_u_a = 0;
    }

    if (*min_u_b != -1 && *min_u_b <= *max_u_a)
    {
      *max_u_a = *max_u_b;
      *min_u_b = -1;
      *max_u_b = -1;
    }
    assert(*min_u_a <= *max_u_a);
    assert(*min_u_b <= *max_u_b);
    assert(*min_u_b == -1 || *max_u_a < *min_u_b);
  }
}
#endif

/**
 * Calculates adjoint NFFT for flag PRE_FULL_PSI.
 * Parallel calculation (OpenMP) with atomic operations.
 *
 * \arg lprod stride (2*m+2)^d
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_B_compute_full_psi(
    C *g, const int *psi_index_g, const R *psi, const C *f,
    const int M, const int d, const int *n, const int m, const int nfft_flags, const int *index_x)
{
  int k;
  int lprod, lprod_m1;
  {
    int t;
    for(t = 0, lprod = 1; t < d; t++)
        lprod *= 2 * m + 2;
  }
  lprod_m1 = lprod / (2 * m + 2);

#ifdef _OPENMP
  if (nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel private(k)
    {
      int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
      const int *ar_x = index_x;
      int n_prod_rest = 1;

      for (k = 1; k < d; k++)
        n_prod_rest *= n[k];

      nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, d, n, m);

      if (min_u_a != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_a || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_a);
#endif
        while (k < M)
        {
          int l0, lrest;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_a || u_prod > max_u_a)
            break;

          for (l0 = 0; l0 < 2 * m + 2; l0++)
          {
            const int start_index = psi_index_g[j * lprod + l0 * lprod_m1];

            if (start_index < my_u0 * n_prod_rest || start_index > (my_o0+1) * n_prod_rest - 1)
              continue;

            for (lrest = 0; lrest < lprod_m1; lrest++)
            {
              const int l = l0 * lprod_m1 + lrest;
              g[psi_index_g[j * lprod + l]] += psi[j * lprod + l] * f[j];
            }
          }

          k++;
        }
      }

      if (min_u_b != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_b || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_b);
#endif
        while (k < M)
        {
          int l0, lrest;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_b || u_prod > max_u_b)
            break;

          for (l0 = 0; l0 < 2 * m + 2; l0++)
          {
            const int start_index = psi_index_g[j * lprod + l0 * lprod_m1];

            if (start_index < my_u0 * n_prod_rest || start_index > (my_o0+1) * n_prod_rest - 1)
              continue;
            for (lrest = 0; lrest < lprod_m1; lrest++)
            {
              const int l = l0 * lprod_m1 + lrest;
              g[psi_index_g[j * lprod + l]] += psi[j * lprod + l] * f[j];
            }
          }

          k++;
        }
      }
    } /* omp parallel */
    return;
  } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < M; k++)
  {
    int l;
    int j = (nfft_flags & NFFT_SORT_NODES) ? index_x[2*k+1] : k;

    for (l = 0; l < lprod; l++)
    {
#ifdef _OPENMP
      C val = psi[j * lprod + l] * f[j];
      C *gref = g + psi_index_g[j * lprod + l];
      R *gref_real = (R*) gref;

      #pragma omp atomic
      gref_real[0] += creal(val);

      #pragma omp atomic
      gref_real[1] += cimag(val);
#else
      g[psi_index_g[j * lprod + l]] += psi[j * lprod + l] * f[j];
#endif
    }
  }
}

MACRO_nfft_B(T)

static inline void nfft_B_openmp_T(nfft_plan *ths)
{
  int lprod; /* 'regular bandwidth' of matrix B  */
  int k;

  memset(ths->g,0,ths->n_total*sizeof(C));

  for (k = 0, lprod = 1; k < ths->d; k++)
    lprod *= (2*ths->m+2);

  if (ths->nfft_flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(ths->g, ths->psi_index_g, ths->psi, ths->f, ths->M_total,
            ths->d, ths->n, ths->m, ths->nfft_flags, ths->index_x);
    return;
  }

  if (ths->nfft_flags & PRE_PSI)
  {
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < ths->M_total; k++)
    {
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int t, t2; /* index dimensions */
      int l_L; /* index one row of B */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (l_L = 0; l_L < lprod; l_L++)
      {
        C *lhs;
        R *lhs_real;
        C val;

        MACRO_update_phi_prod_ll_plain(with_PRE_PSI);

        lhs = ths->g + ll_plain[ths->d];
        lhs_real = (R*)lhs;
        val = phi_prod[ths->d] * ths->f[j];

        #pragma omp atomic
        lhs_real[0] += creal(val);

        #pragma omp atomic
        lhs_real[1] += cimag(val);

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(PRE_PSI) */

  if (ths->nfft_flags & PRE_FG_PSI)
  {
    int t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];
    for(t2 = 0; t2 < ths->d; t2++)
    {
      int lj_fg;
      R tmpEXP2 = EXP(K(-1.0)/ths->b[t2]);
      R tmpEXP2sq = tmpEXP2*tmpEXP2;
      R tmp2 = K(1.0);
      R tmp3 = K(1.0);
      fg_exp_l[t2][0] = K(1.0);
      for(lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++)
      {
        tmp3 = tmp2*tmpEXP2;
        tmp2 *= tmpEXP2sq;
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;
      }
    }

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      R fg_psi[ths->d][2*ths->m+2];
      R tmpEXP1, tmp1;
      int l_fg,lj_fg;
      int l_L;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (t2 = 0; t2 < ths->d; t2++)
      {
        fg_psi[t2][0] = ths->psi[2*(j*ths->d+t2)];
        tmpEXP1 = ths->psi[2*(j*ths->d+t2)+1];
        tmp1 = K(1.0);
        for (l_fg = u[t2]+1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++)
        {
          tmp1 *= tmpEXP1;
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg];
        }
      }

      for (l_L= 0; l_L < lprod; l_L++)
      {
        C *lhs;
        R *lhs_real;
        C val;

        MACRO_update_phi_prod_ll_plain(with_FG_PSI);

        lhs = ths->g + ll_plain[ths->d];
        lhs_real = (R*)lhs;
        val = phi_prod[ths->d] * ths->f[j];

        #pragma omp atomic
        lhs_real[0] += creal(val);

        #pragma omp atomic
        lhs_real[1] += cimag(val);

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(PRE_FG_PSI) */

  if (ths->nfft_flags & FG_PSI)
  {
    int t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];

    nfft_sort_nodes(ths);

    for (t2 = 0; t2 < ths->d; t2++)
    {
      int lj_fg;
      R tmpEXP2 = EXP(K(-1.0)/ths->b[t2]);
      R tmpEXP2sq = tmpEXP2*tmpEXP2;
      R tmp2 = K(1.0);
      R tmp3 = K(1.0);
      fg_exp_l[t2][0] = K(1.0);
      for (lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++)
      {
        tmp3 = tmp2*tmpEXP2;
        tmp2 *= tmpEXP2sq;
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;
      }
    }

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      R fg_psi[ths->d][2*ths->m+2];
      R tmpEXP1, tmp1;
      int l_fg,lj_fg;
      int l_L;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (t2 = 0; t2 < ths->d; t2++)
      {
        fg_psi[t2][0] = (PHI((ths->x[j*ths->d+t2]-((R)u[t2])/ths->n[t2]),t2));

        tmpEXP1 = EXP(K(2.0)*(ths->n[t2]*ths->x[j*ths->d+t2] - u[t2])
          /ths->b[t2]);
        tmp1 = K(1.0);
        for (l_fg = u[t2] + 1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++)
        {
          tmp1 *= tmpEXP1;
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg];
        }
      }

      for (l_L = 0; l_L < lprod; l_L++)
      {
        C *lhs;
        R *lhs_real;
        C val;

        MACRO_update_phi_prod_ll_plain(with_FG_PSI);

        lhs = ths->g + ll_plain[ths->d];
        lhs_real = (R*)lhs;
        val = phi_prod[ths->d] * ths->f[j];

        #pragma omp atomic
        lhs_real[0] += creal(val);

        #pragma omp atomic
        lhs_real[1] += cimag(val);

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(FG_PSI) */

  if (ths->nfft_flags & PRE_LIN_PSI)
  {
    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k<ths->M_total; k++)
    {
      int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
      int t, t2; /* index dimensions */
      int l_L; /* index one row of B */
      int l[ths->d]; /* multi index u<=l<=o */
      int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
      int ll_plain[ths->d+1]; /* postfix plain index in g */
      R phi_prod[ths->d+1]; /* postfix product of PHI */
      R y[ths->d];
      R fg_psi[ths->d][2*ths->m+2];
      int l_fg,lj_fg;
      R ip_w;
      int ip_u;
      int ip_s = ths->K/(ths->m+2);
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      phi_prod[0] = K(1.0);
      ll_plain[0] = 0;

      MACRO_init_uo_l_lj_t;

      for (t2 = 0; t2 < ths->d; t2++)
      {
        y[t2] = ((ths->n[t2]*ths->x[j*ths->d+t2]-(R)u[t2])
          * ((R)ths->K))/(ths->m+2);
        ip_u  = LRINT(FLOOR(y[t2]));
        ip_w  = y[t2]-ip_u;
        for (l_fg = u[t2], lj_fg = 0; l_fg <= o[t2]; l_fg++, lj_fg++)
        {
          fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s)]
            * (1-ip_w) + ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s+1)]
            * (ip_w);
        }
      }

      for (l_L = 0; l_L < lprod; l_L++)
      {
        C *lhs;
        R *lhs_real;
        C val;

        MACRO_update_phi_prod_ll_plain(with_FG_PSI);

        lhs = ths->g + ll_plain[ths->d];
        lhs_real = (R*)lhs;
        val = phi_prod[ths->d] * ths->f[j];

        #pragma omp atomic
        lhs_real[0] += creal(val);

        #pragma omp atomic
        lhs_real[1] += cimag(val);

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */
    } /* for(j) */
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  nfft_sort_nodes(ths);

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < ths->M_total; k++)
  {
    int u[ths->d], o[ths->d]; /* multi band with respect to x_j */
    int t, t2; /* index dimensions */
    int l_L; /* index one row of B */
    int l[ths->d]; /* multi index u<=l<=o */
    int lj[ths->d]; /* multi index 0<=lj<u+o+1 */
    int ll_plain[ths->d+1]; /* postfix plain index in g */
    R phi_prod[ths->d+1]; /* postfix product of PHI */
    int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    phi_prod[0] = K(1.0);
    ll_plain[0] = 0;

    MACRO_init_uo_l_lj_t;

    for (l_L = 0; l_L < lprod; l_L++)
    {
      C *lhs;
      R *lhs_real;
      C val;

      MACRO_update_phi_prod_ll_plain(without_PRE_PSI);

      lhs = ths->g + ll_plain[ths->d];
      lhs_real = (R*)lhs;
      val = phi_prod[ths->d] * ths->f[j];

      #pragma omp atomic
      lhs_real[0] += creal(val);

      #pragma omp atomic
      lhs_real[1] += cimag(val);

      MACRO_count_uo_l_lj_t;
    } /* for(l_L) */
  } /* for(j) */
}

static void nfft_B_T(nfft_plan *ths)
{
#ifdef _OPENMP
  nfft_B_openmp_T(ths);
#else
  nfft_B_serial_T(ths);
#endif
}

/* ## specialized version for d=1  ########################################### */

static void nfft_1d_init_fg_exp_l(R *fg_exp_l, const int m, const R b)
{
  const int tmp2 = 2*m+2;
  int l;
  R fg_exp_b0, fg_exp_b1, fg_exp_b2, fg_exp_b0_sq;

  fg_exp_b0 = EXP(K(-1.0)/b);
  fg_exp_b0_sq = fg_exp_b0*fg_exp_b0;
  fg_exp_b1 = fg_exp_b2 =fg_exp_l[0] = K(1.0);

  for (l = 1; l < tmp2; l++)
  {
    fg_exp_b2 = fg_exp_b1*fg_exp_b0;
    fg_exp_b1 *= fg_exp_b0_sq;
    fg_exp_l[l] = fg_exp_l[l-1]*fg_exp_b2;
  }
}


static void nfft_trafo_1d_compute(C *fj, const C *g,const R *psij_const,
  const R *xj, const int n, const int m)
{
  int u, o, l;
  const C *gj;
  const R *psij;
  psij = psij_const;

  nfft_uo2(&u, &o, *xj, n, m);

  if (u < o)
  {
    for (l = 1, gj = g + u, (*fj) = (*psij++) * (*gj++); l <= 2*m+1; l++)
      (*fj) += (*psij++) * (*gj++);
  }
  else
  {
    for (l = 1, gj = g + u, (*fj) = (*psij++) * (*gj++); l < 2*m+1 - o; l++)
      (*fj) += (*psij++) * (*gj++);
    for (l = 0, gj = g; l <= o; l++)
      (*fj) += (*psij++) * (*gj++);
  }
}

static void nfft_adjoint_1d_compute(const C *fj, C *g,const R *psij_const,
  const R *xj, const int n, const int m)
{
  int u,o,l;
  C *gj;
  const R *psij;
  psij=psij_const;

  nfft_uo2(&u,&o,*xj, n, m);

  if(u<o)
  {
    for (l = 0, gj = g+u; l <= 2*m+1; l++)
      (*gj++) += (*psij++) * (*fj);
  }
  else
  {
    for (l = 0, gj = g+u; l < 2*m+1-o; l++)
      (*gj++) += (*psij++) * (*fj);
    for (l = 0, gj = g; l <= o; l++)
      (*gj++) += (*psij++) * (*fj);
  }
}

#ifdef _OPENMP
/* adjoint NFFT one-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_1d_compute_omp(const C *fj, C *g,const R *psij_const,
  const R *xj, const int n, const int m)
{
  int u,o,l;
  C *gj;
  int index_temp[2*m+2];

  nfft_uo2(&u,&o,*xj, n, m);

  for (l=0; l<=2*m+1; l++)
    index_temp[l] = (l+u)%n;

  for (l = 0, gj = g+u; l <= 2*m+1; l++)
  {
    int i = index_temp[l];
    C *lhs = g+i;
    R *lhs_real = (R*)lhs;
    C val = psij_const[l] * (*fj);
    #pragma omp atomic
    lhs_real[0] += creal(val);

    #pragma omp atomic
    lhs_real[1] += cimag(val);
  }
}
#endif

#ifdef _OPENMP
/**
 * Adjoint NFFT one-dimensional case updating only a specified range of
 * vector g.
 *
 * \arg fg input coefficient f[j]
 * \arg g output vector g
 * \arg psij_const vector of window function values
 * \arg xj node x[j]
 * \arg n FFTW length (number oversampled Fourier coefficients)
 * \arg m window length
 * \arg my_u0 lowest index the current thread writes to in g
 * \arg my_o0 highest index the current thread writes to in g
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_1d_compute_omp3(const C *fj, C *g,const R *psij_const,
  const R *xj, const int n, const int m, const int my_u0, const int my_o0)
{
  int ar_u,ar_o,l;

  nfft_uo2(&ar_u,&ar_o,*xj, n, m);

  if(ar_u<ar_o)
  {
    int u = MAX(my_u0,ar_u);
    int o = MIN(my_o0,ar_o);
    int offset_psij = u-ar_u;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o-u <= 2*m+1);
    assert(offset_psij+o-u <= 2*m+1);
#endif

    for (l = 0; l <= o-u; l++)
      g[u+l] += psij_const[offset_psij+l] * (*fj);
  }
  else
  {
    int u = MAX(my_u0,ar_u);
    int o = my_o0;
    int offset_psij = u-ar_u;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o-u <= 2*m+1);
    assert(offset_psij+o-u <= 2*m+1);
#endif

    for (l = 0; l <= o-u; l++)
      g[u+l] += psij_const[offset_psij+l] * (*fj);

    u = my_u0;
    o = MIN(my_o0,ar_o);
    offset_psij += my_u0-ar_u+n;
#ifdef OMP_ASSERT
if (u<=o)
{
  assert(o-u <= 2*m+1);
  if (offset_psij+o-u > 2*m+1)
  {
    fprintf(stderr, "ERR: %d %d %d %d %d %d %d\n", ar_u, ar_o, my_u0, my_o0, u, o, offset_psij);
  }
  assert(offset_psij+o-u <= 2*m+1);
}
#endif
    for (l = 0; l <= o-u; l++)
      g[u+l] += psij_const[offset_psij+l] * (*fj);
  }
}
#endif

static void nfft_trafo_1d_B(nfft_plan *ths)
{
  const int n = ths->n[0], M = ths->M_total, m = ths->m, m2p2 = 2*m+2;
  const C *g = (C*)ths->g;

  if (ths->nfft_flags & PRE_FULL_PSI)
  {
    int k;
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int l;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < m2p2; l++)
        ths->f[j] += ths->psi[j*m2p2+l] * g[ths->psi_index_g[j*m2p2+l]];
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if (ths->nfft_flags & PRE_PSI)
  {
    int k;
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      nfft_trafo_1d_compute(&ths->f[j], g, ths->psi + j * (2 * m + 2),
        &ths->x[j], n, m);
    }
    return;
  } /* if(PRE_PSI) */

  if (ths->nfft_flags & PRE_FG_PSI)
  {
    int k;
    R fg_exp_l[m2p2];

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      const R fg_psij0 = ths->psi[2 * j], fg_psij1 = ths->psi[2 * j + 1];
      R fg_psij2 = K(1.0);
      R psij_const[m2p2];
      int l;

      psij_const[0] = fg_psij0;

      for (l = 1; l < m2p2; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
      }

      nfft_trafo_1d_compute(&ths->f[j], g, psij_const, &ths->x[j], n, m);
    }

    return;
  } /* if(PRE_FG_PSI) */

  if (ths->nfft_flags & FG_PSI)
  {
    int k;
    R fg_exp_l[m2p2];

    nfft_sort_nodes(ths);

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      int u, o, l;
      R fg_psij0, fg_psij1, fg_psij2;
      R psij_const[m2p2];

      nfft_uo(ths, (int)j, &u, &o, 0);
      fg_psij0 = (PHI(ths->x[j]-((R)u)/n,0));
      fg_psij1 = EXP(K(2.0) * (n * ths->x[j] - u) / ths->b[0]);
      fg_psij2  = K(1.0);

      psij_const[0] = fg_psij0;

      for (l = 1; l < m2p2; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
      }

      nfft_trafo_1d_compute(&ths->f[j], g, psij_const, &ths->x[j], n, m);
    }
    return;
  } /* if(FG_PSI) */

  if (ths->nfft_flags & PRE_LIN_PSI)
  {
    const int K = ths->K, ip_s = K / (m + 2);
    int k;

    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u, o, l;
      R ip_y, ip_w;
      int ip_u;
      R psij_const[m2p2];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      nfft_uo(ths, (int)j, &u, &o, 0);

      ip_y = FABS(n * ths->x[j] - u) * ((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y - ip_u;

      for (l = 0; l < m2p2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w)
          + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w);

      nfft_trafo_1d_compute(&ths->f[j], g, psij_const, &ths->x[j], n, m);
    }
    return;
  } /* if(PRE_LIN_PSI) */
  else
  {
    /* no precomputed psi at all */
    int k;

    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      R psij_const[m2p2];
      int u, o, l;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      nfft_uo(ths, (int)j, &u, &o, 0);

      for (l = 0; l < m2p2; l++)
        psij_const[l] = (PHI(ths->x[j]-((R)((u+l)))/n,0));

      nfft_trafo_1d_compute(&ths->f[j], g, psij_const, &ths->x[j], n, m);
    }
  }
}

#ifdef _OPENMP
/**
 * Determines the blocks of vector g the current thread is responsible for.
 *
 * \arg my_u0 lowest index the current threads writes to in g
 * \arg my_o0 highest index the current threads writes to in g
 * \arg min_u_a lowest index u which could lead to writing to g
 * \arg max_u_a highest index o which could lead to writing to g
 * \arg min_u_b lowest index u which could lead to writing to g
 * \arg max_u_b highest index o which could lead to writing to g
 * \arg n0 FFTW length first component
 * \arg m window length
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_1d_B_omp3_init_(int *my_u0, int *my_o0, int *min_u_a, int *max_u_a, int *min_u_b, int *max_u_b, const int n0, const int m)
{
  int k;
  int nthreads = omp_get_num_threads();
  int nthreads_used = MIN(nthreads, n0);
  int size_per_thread = n0 / nthreads_used;
  int size_left = n0 - size_per_thread * nthreads_used;
  int size_g[nthreads_used];
  int offset_g[nthreads_used];
  int my_id = omp_get_thread_num();

  *min_u_a = -1;
  *max_u_a = -1;
  *min_u_b = -1;
  *max_u_b = -1;
  *my_u0 = -1;
  *my_o0 = -1;

  if (my_id < nthreads_used)
  {
    const int m22 = 2 * m + 2;

    offset_g[0] = 0;
    for (k = 0; k < nthreads_used; k++)
    {
      if (k > 0)
        offset_g[k] = offset_g[k-1] + size_g[k-1];
      size_g[k] = size_per_thread;
      if (size_left > 0)
      {
        size_g[k]++;
	size_left--;
      }
    }

    *my_u0 = offset_g[my_id];
    *my_o0 = offset_g[my_id] + size_g[my_id] - 1;

    if (nthreads_used > 1)
    {
      *max_u_a = offset_g[my_id] + size_g[my_id] - 1;
      *min_u_a = offset_g[my_id] - m22 + 1;
    }
    else
    {
      *min_u_a = 0;
      *max_u_a = n0 - 1;
    }

    if (*min_u_a < 0)
    {
      *min_u_b = offset_g[my_id] - m22 + 1 + n0;
      *max_u_b = n0 - 1;
      *min_u_a = 0;
    }

    if (*min_u_b != -1 && *min_u_b <= *max_u_a)
    {
      *max_u_a = *max_u_b;
      *min_u_b = -1;
      *max_u_b = -1;
    }
    assert(*min_u_a <= *max_u_a);
    assert(*min_u_b <= *max_u_b);
    assert(*min_u_b == -1 || *max_u_a < *min_u_b);
  }
}
#endif

static void nfft_adjoint_1d_B(nfft_plan *ths)
{
  const int n = ths->n[0], M = ths->M_total, m = ths->m;
  int k;
  C *g = (C*)ths->g;

  memset(g,0,ths->n_total*sizeof(C));

  if (ths->nfft_flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, ths->f, M,
        1, ths->n, m, ths->nfft_flags, ths->index_x);
    return;
  } /* if(PRE_FULL_PSI) */

  if (ths->nfft_flags & PRE_PSI)
  {
#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 1, &n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, ths->psi + j * (2 * m + 2), ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, ths->psi + j * (2 * m + 2), ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp(ths->f + j, g, ths->psi + j * (2 * m + 2), ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute(ths->f + j, g, ths->psi + j * (2 * m + 2), ths->x + j, n, m);
#endif
    }

    return;
  } /* if(PRE_PSI) */

  if (ths->nfft_flags & PRE_FG_PSI)
  {
    R fg_exp_l[2 * m + 2];

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;
        R psij_const[2 * m + 2];

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 1, &n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R fg_psij0 = ths->psi[2 * j];
            R fg_psij1 = ths->psi[2 * j + 1];
            R fg_psij2 = K(1.0);

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            psij_const[0] = fg_psij0;
            for (l = 1; l <= 2 * m + 1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
            }

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R fg_psij0 = ths->psi[2 * j];
            R fg_psij1 = ths->psi[2 * j + 1];
            R fg_psij2 = K(1.0);

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            psij_const[0] = fg_psij0;
            for (l = 1; l <= 2 * m + 1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
            }

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */

      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif


    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      R psij_const[2 * m + 2];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      int l;
      R fg_psij0 = ths->psi[2 * j];
      R fg_psij1 = ths->psi[2 * j + 1];
      R fg_psij2 = K(1.0);

      psij_const[0] = fg_psij0;
      for (l = 1; l <= 2 * m + 1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
      }

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp(ths->f + j, g, psij_const, ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if (ths->nfft_flags & FG_PSI)
  {
    R fg_exp_l[2 * m + 2];

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

    nfft_sort_nodes(ths);

#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;
        R psij_const[2 * m + 2];
        R fg_psij0, fg_psij1, fg_psij2;

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 1, &n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_uo(ths, j, &u, &o, 0);
            fg_psij0 = (PHI(ths->x[j]-((R)u)/n,0));
            fg_psij1 = EXP(K(2.0) * (n * (ths->x[j]) - u) / ths->b[0]);
            fg_psij2 = K(1.0);
            psij_const[0] = fg_psij0;
            for (l = 1; l <= 2 * m + 1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
            }

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_uo(ths, j, &u, &o, 0);
            fg_psij0 = (PHI(ths->x[j]-((R)u)/n,0));
            fg_psij1 = EXP(K(2.0) * (n * (ths->x[j]) - u) / ths->b[0]);
            fg_psij2 = K(1.0);
            psij_const[0] = fg_psij0;
            for (l = 1; l <= 2 * m + 1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
            }

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */

      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u,o,l;
      R psij_const[2 * m + 2];
      R fg_psij0, fg_psij1, fg_psij2;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      nfft_uo(ths, j, &u, &o, 0);
      fg_psij0 = (PHI(ths->x[j]-((R)u)/n,0));
      fg_psij1 = EXP(K(2.0) * (n * (ths->x[j]) - u) / ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for (l = 1; l <= 2 * m + 1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
      }

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp(ths->f + j, g, psij_const, ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if (ths->nfft_flags & PRE_LIN_PSI)
  {
    const int K = ths->K;
    const int ip_s = K / (m + 2);

    nfft_sort_nodes(ths);

#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;
        R psij_const[2 * m + 2];
        int ip_u;
        R ip_y, ip_w;

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 1, &n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_uo(ths, j, &u, &o, 0);

            ip_y = FABS(n * ths->x[j] - u) * ((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y - ip_u;
            for (l = 0; l < 2 * m + 2; l++)
              psij_const[l]
                  = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w)
                      + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w);

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_uo(ths, j, &u, &o, 0);

            ip_y = FABS(n * ths->x[j] - u) * ((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y - ip_u;
            for (l = 0; l < 2 * m + 2; l++)
              psij_const[l]
                  = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w)
                      + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w);

            nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma openmp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u,o,l;
      int ip_u;
      R ip_y, ip_w;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[2 * m + 2];

      nfft_uo(ths, j, &u, &o, 0);

      ip_y = FABS(n * ths->x[j] - u) * ((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y - ip_u;
      for (l = 0; l < 2 * m + 2; l++)
        psij_const[l]
            = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w)
                + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w);

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp(ths->f + j, g, psij_const, ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  nfft_sort_nodes(ths);

#ifdef _OPENMP
  if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel private(k)
    {
      int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
      int *ar_x = ths->index_x;
      R psij_const[2 * m + 2];

      nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 1, &n, m);

      if (min_u_a != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_a || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_a);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_a || u_prod > max_u_a)
            break;

          nfft_uo(ths, j, &u, &o, 0);

          for (l = 0; l <= 2 * m + 1; l++)
            psij_const[l] = (PHI(ths->x[j]-((R)((u+l)))/n,0));

          nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

          k++;
        }
      }

      if (min_u_b != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_b || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_b);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_b || u_prod > max_u_b)
            break;

          nfft_uo(ths, j, &u, &o, 0);

          for (l = 0; l <= 2 * m + 1; l++)
            psij_const[l] = (PHI(ths->x[j]-((R)((u+l)))/n,0));

          nfft_adjoint_1d_compute_omp3(ths->f + j, g, psij_const, ths->x + j, n, m, my_u0, my_o0);

          k++;
        }
      }
    } /* omp parallel */

    return;
  } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < M; k++)
  {
    int u,o,l;
    R psij_const[2 * m + 2];
    int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    nfft_uo(ths, j, &u, &o, 0);

    for (l = 0; l <= 2 * m + 1; l++)
      psij_const[l] = (PHI(ths->x[j]-((R)((u+l)))/n,0));

#ifdef _OPENMP
    nfft_adjoint_1d_compute_omp(ths->f + j, g, psij_const, ths->x + j, n, m);
#else
    nfft_adjoint_1d_compute(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
  }
}

void nfft_trafo_1d(nfft_plan *ths)
{
  const int N = ths->N[0], N2 = N/2, n = ths->n[0];
  C *f_hat1 = (C*)ths->f_hat, *f_hat2 = (C*)&ths->f_hat[N2];

  ths->g_hat = ths->g1;
  ths->g = ths->g2;

  {
    C *g_hat1 = (C*)&ths->g_hat[n-N/2], *g_hat2 = (C*)ths->g_hat;
    R *c_phi_inv1, *c_phi_inv2;

    TIC(0)
#ifdef _OPENMP
    {
      int k;
      #pragma omp parallel for default(shared) private(k)
      for (k = 0; k < ths->n_total; k++)
        ths->g_hat[k] = 0.0;
    }
#else
    memset(ths->g_hat, 0, ths->n_total*sizeof(C));
#endif
    if(ths->nfft_flags & PRE_PHI_HUT)
    {
      int k;
      c_phi_inv1 = ths->c_phi_inv[0];
      c_phi_inv2 = &ths->c_phi_inv[0][N2];

      #pragma omp parallel for default(shared) private(k)
      for (k = 0; k < N2; k++)
      {
        g_hat1[k] = f_hat1[k] * c_phi_inv1[k];
        g_hat2[k] = f_hat2[k] * c_phi_inv2[k];
      }
    }
    else
    {
      int k;
      #pragma omp parallel for default(shared) private(k)
      for (k = 0; k < N2; k++)
      {
        g_hat1[k] = f_hat1[k] / (PHI_HUT(k-N2,0));
        g_hat2[k] = f_hat2[k] / (PHI_HUT(k,0));
      }
    }
    TOC(0)

    TIC_FFTW(1)
    fftw_execute(ths->my_fftw_plan1);
    TOC_FFTW(1);

    TIC(2);
    nfft_trafo_1d_B(ths);
    TOC(2);
  }
}

void nfft_adjoint_1d(nfft_plan *ths)
{
  int n,N;
  C *g_hat1,*g_hat2,*f_hat1,*f_hat2;
  R *c_phi_inv1, *c_phi_inv2;

  N=ths->N[0];
  n=ths->n[0];

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  f_hat1=(C*)ths->f_hat;
  f_hat2=(C*)&ths->f_hat[N/2];
  g_hat1=(C*)&ths->g_hat[n-N/2];
  g_hat2=(C*)ths->g_hat;

  TIC(2)
  nfft_adjoint_1d_B(ths);
  TOC(2)

  TIC_FFTW(1)
  fftw_execute(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      int k;
      c_phi_inv1=ths->c_phi_inv[0];
      c_phi_inv2=&ths->c_phi_inv[0][N/2];
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k)
      for (k = 0; k < N/2; k++)
      {
        f_hat1[k] = g_hat1[k] * c_phi_inv1[k];
        f_hat2[k] = g_hat2[k] * c_phi_inv2[k];
      }
#else
      for(k=0;k<N/2;k++)
  {
    (*f_hat1++) = (*g_hat1++) * (*c_phi_inv1++);
    (*f_hat2++) = (*g_hat2++) * (*c_phi_inv2++);
  }
#endif
    }
  else
  {
    int k;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < N/2; k++)
    {
      f_hat1[k] = g_hat1[k] / (PHI_HUT(k-N/2,0));
      f_hat2[k] = g_hat2[k] / (PHI_HUT(k,0));
    }
#else
    for(k=0;k<N/2;k++)
      {
  (*f_hat1++) = (*g_hat1++) / (PHI_HUT(k-N/2,0));
  (*f_hat2++) = (*g_hat2++) / (PHI_HUT(k,0));
      }
#endif
  }
  TOC(0)
}


/* ############################################################ SPECIFIC VERSIONS FOR d=2 */

static void nfft_2d_init_fg_exp_l(R *fg_exp_l, const int m, const R b)
{
  int l;
  R fg_exp_b0, fg_exp_b1, fg_exp_b2, fg_exp_b0_sq;

  fg_exp_b0 = EXP(K(-1.0)/b);
  fg_exp_b0_sq = fg_exp_b0*fg_exp_b0;
  fg_exp_b1 = K(1.0);
  fg_exp_b2 = K(1.0);
  fg_exp_l[0] = K(1.0);
  for(l=1; l <= 2*m+1; l++)
    {
      fg_exp_b2 = fg_exp_b1*fg_exp_b0;
      fg_exp_b1 *= fg_exp_b0_sq;
      fg_exp_l[l] = fg_exp_l[l-1]*fg_exp_b2;
    }
}

static void nfft_trafo_2d_compute(C *fj, const C *g,
          const R *psij_const0, const R *psij_const1,
          const R *xj0, const R *xj1,
          const int n0, const int n1, const int m)
{
  int u0,o0,l0,u1,o1,l1;
  const C *gj;
  const R *psij0,*psij1;

  psij0=psij_const0;
  psij1=psij_const1;

  nfft_uo2(&u0,&o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);

  *fj=0;

  if(u0<o0)
      if(u1<o1)
    for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<=2*m+1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
      else
    for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<2*m+1-o1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
        gj=g+(u0+l0)*n1;
        for(l1=0; l1<=o1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
  else
      if(u1<o1)
      {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<=2*m+1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
    for(l0=0; l0<=o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+l0*n1+u1;
        for(l1=0; l1<=2*m+1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
      }
      else
      {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<2*m+1-o1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
        gj=g+(u0+l0)*n1;
        for(l1=0; l1<=o1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
    for(l0=0; l0<=o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+l0*n1+u1;
        for(l1=0; l1<2*m+1-o1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
        gj=g+l0*n1;
        for(l1=0; l1<=o1; l1++)
      (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
      }
}

#ifdef _OPENMP
/* adjoint NFFT two-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_2d_compute_omp(const C *fj, C *g,
            const R *psij_const0, const R *psij_const1,
            const R *xj0, const R *xj1,
            const int n0, const int n1, const int m)
{
  int u0,o0,l0,u1,o1,l1;
  const int lprod = (2*m+2) * (2*m+2); 

  unsigned long int index_temp0[2*m+2];
  unsigned long int index_temp1[2*m+2];

  nfft_uo2(&u0,&o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);

  for (l0=0; l0<=2*m+1; l0++)
    index_temp0[l0] = (u0+l0)%n0;

  for (l1=0; l1<=2*m+1; l1++)
    index_temp1[l1] = (u1+l1)%n1;

  for(l0=0; l0<=2*m+1; l0++)
  {
    for(l1=0; l1<=2*m+1; l1++)
    {
      unsigned long int i = index_temp0[l0] * n1 + index_temp1[l1];
      C *lhs = g+i;
      R *lhs_real = (R*)lhs;
      C val = psij_const0[l0] * psij_const1[l1] * (*fj);

      #pragma omp atomic
      lhs_real[0] += creal(val);

      #pragma omp atomic
      lhs_real[1] += cimag(val);
    }
  }
}
#endif

#ifdef _OPENMP
/** 
 * Adjoint NFFT two-dimensional case updating only a specified range of
 * vector g.
 *
 * \arg fg input coefficient f[j]
 * \arg g output vector g
 * \arg psij_const0 vector of window function values first component
 * \arg psij_const1 vector of window function values second component
 * \arg xj0 node x[2*j]
 * \arg xj1 node x[2*j+1]
 * \arg n0 FFTW length (number oversampled Fourier coefficients) first comp.
 * \arg n1 FFTW length (number oversampled Fourier coefficients) second comp.
 * \arg m window length
 * \arg my_u0 lowest index (first component) the current thread writes to
 * \arg my_o0 highest index (second component) the current thread writes to
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_2d_compute_omp3(const C *fj, C *g,
            const R *psij_const0, const R *psij_const1,
            const R *xj0, const R *xj1,
            const int n0, const int n1, const int m,
	    const int my_u0, const int my_o0)
{
  int ar_u0,ar_o0,l0,u1,o1,l1;
  const int lprod = (2*m+2) * (2*m+2); 
  unsigned long int index_temp1[2*m+2];

  nfft_uo2(&ar_u0,&ar_o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);

  for (l1=0; l1<=2*m+1; l1++)
    index_temp1[l1] = (u1+l1)%n1;

  if(ar_u0<ar_o0)
  {
    int u0 = MAX(my_u0,ar_u0);
    int o0 = MIN(my_o0,ar_o0);
    int offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      unsigned long int i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
        g[i0 + index_temp1[l1]] += val0 * psij_const1[l1] * (*fj);
    }
  }
  else
  {
    int u0 = MAX(my_u0,ar_u0);
    int o0 = my_o0;
    int offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      unsigned long int i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
        g[i0 + index_temp1[l1]] += val0 * psij_const1[l1] * (*fj);
    }

    u0 = my_u0;
    o0 = MIN(my_o0,ar_o0);
    offset_psij += my_u0-ar_u0+n0;
#ifdef OMP_ASSERT
if (u0<=o0)
{
  assert(o0-u0 <= 2*m+1);
  assert(offset_psij+o0-u0 <= 2*m+1);
}
#endif
    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      unsigned long int i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
        g[i0 + index_temp1[l1]] += val0 * psij_const1[l1] * (*fj);
    }
  }
}
#endif

static void nfft_adjoint_2d_compute(const C *fj, C *g,
            const R *psij_const0, const R *psij_const1,
            const R *xj0, const R *xj1,
            const int n0, const int n1, const int m)
{
  int u0,o0,l0,u1,o1,l1;
  C *gj;
  const R *psij0,*psij1;

  psij0=psij_const0;
  psij1=psij_const1;

  nfft_uo2(&u0,&o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);

  if(u0<o0)
      if(u1<o1)
    for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<=2*m+1; l1++)
    (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
      else
    for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<2*m+1-o1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
        gj=g+(u0+l0)*n1;
        for(l1=0; l1<=o1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
  else
      if(u1<o1)
      {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<=2*m+1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
    for(l0=0; l0<=o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+l0*n1+u1;
        for(l1=0; l1<=2*m+1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
      }
      else
      {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+(u0+l0)*n1+u1;
        for(l1=0; l1<2*m+1-o1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
        gj=g+(u0+l0)*n1;
        for(l1=0; l1<=o1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
    for(l0=0; l0<=o0; l0++,psij0++)
    {
        psij1=psij_const1;
        gj=g+l0*n1+u1;
        for(l1=0; l1<2*m+1-o1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
        gj=g+l0*n1;
        for(l1=0; l1<=o1; l1++)
      (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
      }
}

static void nfft_trafo_2d_B(nfft_plan *ths)
{
  const C *g = (C*)ths->g;
  const int N0 = ths->N[0];
  const int n0 = ths->n[0];
  const int N1 = ths->N[1];
  const int n1 = ths->n[1];
  const int M = ths->M_total;
  const int m = ths->m;

  int k;

  if(ths->nfft_flags & PRE_FULL_PSI)
  {
    const int lprod = (2*m+2) * (2*m+2);
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int l;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < lprod; l++)
        ths->f[j] += ths->psi[j*lprod+l] * g[ths->psi_index_g[j*lprod+l]];
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
  {
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      nfft_trafo_2d_compute(ths->f+j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m);
    }

      return;
  } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+2)];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      int l;
      R fg_psij0 = ths->psi[2*j*2];
      R fg_psij1 = ths->psi[2*j*2+1];
      R fg_psij2 = K(1.0);

      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      fg_psij0 = ths->psi[2*(j*2+1)];
      fg_psij1 = ths->psi[2*(j*2+1)+1];
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      nfft_trafo_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u, o, l;
      R fg_psij0, fg_psij1, fg_psij2;
      R psij_const[2*(2*m+2)];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      nfft_uo(ths,j,&u,&o,0);
      fg_psij0 = (PHI(ths->x[2*j]-((R)u)/n0,0));
      fg_psij1 = EXP(K(2.0)*(n0*(ths->x[2*j]) - u)/ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      nfft_uo(ths,j,&u,&o,1);
      fg_psij0 = (PHI(ths->x[2*j+1]-((R)u)/n1,1));
      fg_psij1 = EXP(K(2.0)*(n1*(ths->x[2*j+1]) - u)/ths->b[1]);
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      nfft_trafo_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
    }

    return;
  } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
  {
    const int K = ths->K, ip_s = K / (m + 2);

    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u, o, l;
      R ip_y, ip_w;
      int ip_u;
      R psij_const[2*(2*m+2)];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      nfft_uo(ths,j,&u,&o,0);
      ip_y = FABS(n0*ths->x[2*j] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_uo(ths,j,&u,&o,1);
      ip_y = FABS(n1*ths->x[2*j+1] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_trafo_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
    }
      return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */

  nfft_sort_nodes(ths);

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < M; k++)
  {
    R psij_const[2*(2*m+2)];
    int u, o, l;
    int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    nfft_uo(ths,j,&u,&o,0);
    for(l=0;l<=2*m+1;l++)
      psij_const[l]=(PHI(ths->x[2*j]-((R)((u+l)))/n0,0));

    nfft_uo(ths,j,&u,&o,1);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*m+2+l]=(PHI(ths->x[2*j+1]-((R)((u+l)))/n1,1));

    nfft_trafo_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
  }
}

#ifdef _OPENMP
/**
 * Determines the blocks of vector g the current thread is responsible for.
 *
 * \arg my_u0 lowest index (first component) the current threads writes to in g
 * \arg my_o0 highest index (first component) the current threads writes to in g
 * \arg min_u_a lowest (linear 2d) index u which could lead to writing to g
 * \arg max_u_a highest (linear 2d) index o which could lead to writing to g
 * \arg min_u_b lowest (linear 2d) index u which could lead to writing to g
 * \arg max_u_b highest (linear 2d) index o which could lead to writing to g
 * \arg n0 FFTW length first component
 * \arg n1 FFTW length second component
 * \arg m window length
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_2d_B_omp3_init_(int *my_u0, int *my_o0, int *min_u_a, int *max_u_a, int *min_u_b, int *max_u_b, const int n0, const int n1, const int m)
{
  int k;
  int nthreads = omp_get_num_threads();
  int nthreads_used = MIN(nthreads, n0);
  int size_per_thread = n0 / nthreads_used;
  int size_left = n0 - size_per_thread * nthreads_used;
  int size_g[nthreads_used];
  int offset_g[nthreads_used];
  int my_id = omp_get_thread_num();

  *min_u_a = -1;
  *max_u_a = -1;
  *min_u_b = -1;
  *max_u_b = -1;
  *my_u0 = -1;
  *my_o0 = -1;

  if (my_id < nthreads_used)
  {
    const int m22 = 2 * m + 2;

    offset_g[0] = 0;
    for (k = 0; k < nthreads_used; k++)
    {
      if (k > 0)
        offset_g[k] = offset_g[k-1] + size_g[k-1];
      size_g[k] = size_per_thread;
      if (size_left > 0)
      {
        size_g[k]++;
	size_left--;
      }
    }

    *my_u0 = offset_g[my_id];
    *my_o0 = offset_g[my_id] + size_g[my_id] - 1;

    if (nthreads_used > 1)
    {
      *max_u_a = n1*(offset_g[my_id] + size_g[my_id]) - 1;
      *min_u_a = n1*(offset_g[my_id] - m22 + 1);
    }
    else
    {
      *min_u_a = 0;
      *max_u_a = n1 * n0 - 1;
    }

    if (*min_u_a < 0)
    {
      *min_u_b = n1 * (offset_g[my_id] - m22 + 1 + n0);
      *max_u_b = n1 * n0 - 1;
      *min_u_a = 0;
    }

    if (*min_u_b != -1 && *min_u_b <= *max_u_a)
    {
      *max_u_a = *max_u_b;
      *min_u_b = -1;
      *max_u_b = -1;
    }
    assert(*min_u_a <= *max_u_a);
    assert(*min_u_b <= *max_u_b);
    assert(*min_u_b == -1 || *max_u_a < *min_u_b);
  }
}
#endif

static void nfft_adjoint_2d_B(nfft_plan *ths)
{
  const int N0 = ths->N[0];
  const int n0 = ths->n[0];
  const int N1 = ths->N[1];
  const int n1 = ths->n[1];
  const int M = ths->M_total;
  const int m = ths->m;
  C* g = (C*) ths->g;
  int k;

  memset(g,0,ths->n_total*sizeof(C));

  if(ths->nfft_flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, ths->f, M,
        2, ths->n, m, ths->nfft_flags, ths->index_x);
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
  {
#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 2, ths->n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_adjoint_2d_compute_omp3(ths->f+j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          int k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_adjoint_2d_compute_omp3(ths->f+j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp(ths->f+j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute(ths->f+j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

#ifdef _OPENMP
  if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel private(k)
    {
      int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
      int *ar_x = ths->index_x;
      R psij_const[2*(2*m+2)];
      R fg_psij0, fg_psij1, fg_psij2;

      nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 2, ths->n, m);

      if (min_u_a != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_a || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_a);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];
          R fg_psij0 = ths->psi[2*j*2];
          R fg_psij1 = ths->psi[2*j*2+1];
          R fg_psij2 = K(1.0);

          if (u_prod < min_u_a || u_prod > max_u_a)
            break;

          psij_const[0] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
          }

          fg_psij0 = ths->psi[2*(j*2+1)];
          fg_psij1 = ths->psi[2*(j*2+1)+1];
          fg_psij2 = K(1.0);
          psij_const[2*m+2] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
          }

          nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

          k++;
        }
      }

      if (min_u_b != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_b || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_b);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];
          R fg_psij0 = ths->psi[2*j*2];
          R fg_psij1 = ths->psi[2*j*2+1];
          R fg_psij2 = K(1.0);

          if (u_prod < min_u_b || u_prod > max_u_b)
            break;

          psij_const[0] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
          }

          fg_psij0 = ths->psi[2*(j*2+1)];
          fg_psij1 = ths->psi[2*(j*2+1)+1];
          fg_psij2 = K(1.0);
          psij_const[2*m+2] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
          }

          nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

          k++;
        }
      }
    } /* omp parallel */
    return;
  } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif


    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+2)];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      int l;
      R fg_psij0 = ths->psi[2*j*2];
      R fg_psij1 = ths->psi[2*j*2+1];
      R fg_psij2 = K(1.0);

      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      fg_psij0 = ths->psi[2*(j*2+1)];
      fg_psij1 = ths->psi[2*(j*2+1)+1];
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

    nfft_sort_nodes(ths);

#ifdef _OPENMP
  if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel private(k)
    {
      int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
      int *ar_x = ths->index_x;
      R psij_const[2*(2*m+2)];
      R fg_psij0, fg_psij1, fg_psij2;

      nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 2, ths->n, m);

      if (min_u_a != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_a || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_a);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_a || u_prod > max_u_a)
            break;

          nfft_uo(ths,j,&u,&o,0);
          fg_psij0 = (PHI(ths->x[2*j]-((R)u)/n0,0));
          fg_psij1 = EXP(K(2.0)*(n0*(ths->x[2*j]) - u)/ths->b[0]);
          fg_psij2 = K(1.0);
          psij_const[0] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
          }

          nfft_uo(ths,j,&u,&o,1);
          fg_psij0 = (PHI(ths->x[2*j+1]-((R)u)/n1,1));
          fg_psij1 = EXP(K(2.0)*(n1*(ths->x[2*j+1]) - u)/ths->b[1]);
          fg_psij2 = K(1.0);
          psij_const[2*m+2] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
          }

          nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

          k++;
        }
      }

      if (min_u_b != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_b || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_b);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_b || u_prod > max_u_b)
            break;

          nfft_uo(ths,j,&u,&o,0);
          fg_psij0 = (PHI(ths->x[2*j]-((R)u)/n0,0));
          fg_psij1 = EXP(K(2.0)*(n0*(ths->x[2*j]) - u)/ths->b[0]);
          fg_psij2 = K(1.0);
          psij_const[0] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
          }

          nfft_uo(ths,j,&u,&o,1);
          fg_psij0 = (PHI(ths->x[2*j+1]-((R)u)/n1,1));
          fg_psij1 = EXP(K(2.0)*(n1*(ths->x[2*j+1]) - u)/ths->b[1]);
          fg_psij2 = K(1.0);
          psij_const[2*m+2] = fg_psij0;
          for(l=1; l<=2*m+1; l++)
          {
            fg_psij2 *= fg_psij1;
            psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
          }

          nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

          k++;
        }
      }
    } /* omp parallel */
    return;
  } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u, o, l;
      R fg_psij0, fg_psij1, fg_psij2;
      R psij_const[2*(2*m+2)];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      nfft_uo(ths,j,&u,&o,0);
      fg_psij0 = (PHI(ths->x[2*j]-((R)u)/n0,0));
      fg_psij1 = EXP(K(2.0)*(n0*(ths->x[2*j]) - u)/ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      nfft_uo(ths,j,&u,&o,1);
      fg_psij0 = (PHI(ths->x[2*j+1]-((R)u)/n1,1));
      fg_psij1 = EXP(K(2.0)*(n1*(ths->x[2*j+1]) - u)/ths->b[1]);
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
  {
    const int K = ths->K;
    const int ip_s = K / (m + 2);

    nfft_sort_nodes(ths);

#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;
        R psij_const[2*(2*m+2)];

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 2, ths->n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int ip_u;
            R ip_y, ip_w;
            int j = ar_x[2*k+1];

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_uo(ths,j,&u,&o,0);
            ip_y = FABS(n0*(ths->x[2*j]) - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_uo(ths,j,&u,&o,1);
            ip_y = FABS(n1*(ths->x[2*j+1]) - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int ip_u;
            R ip_y, ip_w;
            int j = ar_x[2*k+1];

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_uo(ths,j,&u,&o,0);
            ip_y = FABS(n0*(ths->x[2*j]) - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_uo(ths,j,&u,&o,1);
            ip_y = FABS(n1*(ths->x[2*j+1]) - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma openmp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u,o,l;
      int ip_u;
      R ip_y, ip_w;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[2*(2*m+2)];

      nfft_uo(ths,j,&u,&o,0);
      ip_y = FABS(n0*(ths->x[2*j]) - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_uo(ths,j,&u,&o,1);
      ip_y = FABS(n1*(ths->x[2*j+1]) - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
  }
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  nfft_sort_nodes(ths);

#ifdef _OPENMP
  if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel private(k)
    {
      int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
      int *ar_x = ths->index_x;
      R psij_const[2*(2*m+2)];

      nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 2, ths->n, m);

      if (min_u_a != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_a || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_a);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_a || u_prod > max_u_a)
            break;

          nfft_uo(ths,j,&u,&o,0);
          for(l=0;l<=2*m+1;l++)
            psij_const[l]=(PHI(ths->x[2*j]-((R)((u+l)))/n0,0));

          nfft_uo(ths,j,&u,&o,1);
          for(l=0;l<=2*m+1;l++)
            psij_const[2*m+2+l]=(PHI(ths->x[2*j+1]-((R)((u+l)))/n1,1));

          nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

          k++;
        }
      }

      if (min_u_b != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_b || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_b);
#endif
        while (k < M)
        {
          int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];

          if (u_prod < min_u_b || u_prod > max_u_b)
            break;

          nfft_uo(ths,j,&u,&o,0);
          for(l=0;l<=2*m+1;l++)
            psij_const[l]=(PHI(ths->x[2*j]-((R)((u+l)))/n0,0));

          nfft_uo(ths,j,&u,&o,1);
          for(l=0;l<=2*m+1;l++)
            psij_const[2*m+2+l]=(PHI(ths->x[2*j+1]-((R)((u+l)))/n1,1));

          nfft_adjoint_2d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

          k++;
        }
      }
    } /* omp parallel */
    return;
  } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif


  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < M; k++)
  {
    int u,o,l;
    R psij_const[2*(2*m+2)];
    int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    nfft_uo(ths,j,&u,&o,0);
    for(l=0;l<=2*m+1;l++)
      psij_const[l]=(PHI(ths->x[2*j]-((R)((u+l)))/n0,0));

    nfft_uo(ths,j,&u,&o,1);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*m+2+l]=(PHI(ths->x[2*j+1]-((R)((u+l)))/n1,1));

#ifdef _OPENMP
    nfft_adjoint_2d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
    nfft_adjoint_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
  }
}


void nfft_trafo_2d(nfft_plan *ths)
{
  int k0,k1,n0,n1,N0,N1;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  R ck01, ck02, ck11, ck12;
  C *g_hat11,*f_hat11,*g_hat21,*f_hat21,*g_hat12,*f_hat12,*g_hat22,*f_hat22;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(0)
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k0)
  for (k0 = 0; k0 < ths->n_total; k0++)
    ths->g_hat[k0] = 0.0;
#else
  memset(ths->g_hat,0,ths->n_total*sizeof(C));
#endif
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,c_phi_inv11,c_phi_inv12,g_hat11,f_hat11,g_hat21,f_hat21,g_hat12,f_hat12,g_hat22,f_hat22,ck11,ck12)
      for(k0=0;k0<N0/2;k0++)
      {
        ck01=c_phi_inv01[k0];
        ck02=c_phi_inv02[k0];

        c_phi_inv11=ths->c_phi_inv[1];
        c_phi_inv12=&ths->c_phi_inv[1][N1/2];

        g_hat11=g_hat + (n0-(N0/2)+k0)*n1+n1-(N1/2);
        f_hat11=f_hat + k0*N1;
        g_hat21=g_hat + k0*n1+n1-(N1/2);
        f_hat21=f_hat + ((N0/2)+k0)*N1;
        g_hat12=g_hat + (n0-(N0/2)+k0)*n1;
        f_hat12=f_hat + k0*N1+(N1/2);
        g_hat22=g_hat + k0*n1;
        f_hat22=f_hat + ((N0/2)+k0)*N1+(N1/2);

        for(k1=0;k1<N1/2;k1++)
        {
        ck11=c_phi_inv11[k1];
        ck12=c_phi_inv12[k1];

        g_hat11[k1] = f_hat11[k1] * ck01 * ck11;
        g_hat21[k1] = f_hat21[k1] * ck02 * ck11;
        g_hat12[k1] = f_hat12[k1] * ck01 * ck12;
        g_hat22[k1] = f_hat22[k1] * ck02 * ck12;
      }
  }
#else
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=(*c_phi_inv01++);
    ck02=(*c_phi_inv02++);

    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];

    g_hat11=g_hat + (n0-(N0/2)+k0)*n1+n1-(N1/2);
    f_hat11=f_hat + k0*N1;
          g_hat21=g_hat + k0*n1+n1-(N1/2);
          f_hat21=f_hat + ((N0/2)+k0)*N1;
          g_hat12=g_hat + (n0-(N0/2)+k0)*n1;
          f_hat12=f_hat + k0*N1+(N1/2);
    g_hat22=g_hat + k0*n1;
    f_hat22=f_hat + ((N0/2)+k0)*N1+(N1/2);
    for(k1=0;k1<N1/2;k1++)
      {
        ck11=(*c_phi_inv11++);
        ck12=(*c_phi_inv12++);

        (*g_hat11++) = (*f_hat11++) * ck01 * ck11;
        (*g_hat21++) = (*f_hat21++) * ck02 * ck11;
        (*g_hat12++) = (*f_hat12++) * ck01 * ck12;
        (*g_hat22++) = (*f_hat22++) * ck02 * ck12;
      }
  }
#endif
    }
  else
    #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,ck11,ck12)
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(k1,1));
      g_hat[(n0-N0/2+k0)*n1+n1-N1/2+k1] = f_hat[k0*N1+k1]             * ck01 * ck11;
      g_hat[k0*n1+n1-N1/2+k1]           = f_hat[(N0/2+k0)*N1+k1]      * ck02 * ck11;
      g_hat[(n0-N0/2+k0)*n1+k1]         = f_hat[k0*N1+N1/2+k1]        * ck01 * ck12;
      g_hat[k0*n1+k1]                   = f_hat[(N0/2+k0)*N1+N1/2+k1] * ck02 * ck12;
    }
      }

  TOC(0)

  TIC_FFTW(1)
  fftw_execute(ths->my_fftw_plan1);
  TOC_FFTW(1);

  TIC(2);
  nfft_trafo_2d_B(ths);
  TOC(2);
}

void nfft_adjoint_2d(nfft_plan *ths)
{
  int k0,k1,n0,n1,N0,N1;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  R ck01, ck02, ck11, ck12;
  C *g_hat11,*f_hat11,*g_hat21,*f_hat21,*g_hat12,*f_hat12,*g_hat22,*f_hat22;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(2);
  nfft_adjoint_2d_B(ths);
  TOC(2);

  TIC_FFTW(1)
  fftw_execute(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,c_phi_inv11,c_phi_inv12,g_hat11,f_hat11,g_hat21,f_hat21,g_hat12,f_hat12,g_hat22,f_hat22,ck11,ck12)
      for(k0=0;k0<N0/2;k0++)
      {
        ck01=c_phi_inv01[k0];
        ck02=c_phi_inv02[k0];

        c_phi_inv11=ths->c_phi_inv[1];
        c_phi_inv12=&ths->c_phi_inv[1][N1/2];

        g_hat11=g_hat + (n0-(N0/2)+k0)*n1+n1-(N1/2);
        f_hat11=f_hat + k0*N1;
        g_hat21=g_hat + k0*n1+n1-(N1/2);
        f_hat21=f_hat + ((N0/2)+k0)*N1;
        g_hat12=g_hat + (n0-(N0/2)+k0)*n1;
        f_hat12=f_hat + k0*N1+(N1/2);
        g_hat22=g_hat + k0*n1;
        f_hat22=f_hat + ((N0/2)+k0)*N1+(N1/2);

        for(k1=0;k1<N1/2;k1++)
        {
          ck11=c_phi_inv11[k1];
          ck12=c_phi_inv12[k1];

          f_hat11[k1] = g_hat11[k1] * ck01 * ck11;
          f_hat21[k1] = g_hat21[k1] * ck02 * ck11;
          f_hat12[k1] = g_hat12[k1] * ck01 * ck12;
          f_hat22[k1] = g_hat22[k1] * ck02 * ck12;
        }
      }
#else
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=(*c_phi_inv01++);
    ck02=(*c_phi_inv02++);

    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];
    g_hat11=g_hat + (n0-(N0/2)+k0)*n1+n1-(N1/2);
    f_hat11=f_hat + k0*N1;
          g_hat21=g_hat + k0*n1+n1-(N1/2);
          f_hat21=f_hat + ((N0/2)+k0)*N1;
          g_hat12=g_hat + (n0-(N0/2)+k0)*n1;
          f_hat12=f_hat + k0*N1+(N1/2);
    g_hat22=g_hat + k0*n1;
    f_hat22=f_hat + ((N0/2)+k0)*N1+(N1/2);
    for(k1=0;k1<N1/2;k1++)
      {
        ck11=(*c_phi_inv11++);
        ck12=(*c_phi_inv12++);

        (*f_hat11++) = (*g_hat11++) * ck01 * ck11;
        (*f_hat21++) = (*g_hat21++) * ck02 * ck11;
        (*f_hat12++) = (*g_hat12++) * ck01 * ck12;
        (*f_hat22++) = (*g_hat22++) * ck02 * ck12;
      }
  }
#endif
    }
  else
    #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,ck11,ck12)
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(k1,1));
      f_hat[k0*N1+k1]             = g_hat[(n0-N0/2+k0)*n1+n1-N1/2+k1] * ck01 * ck11;
      f_hat[(N0/2+k0)*N1+k1]      = g_hat[k0*n1+n1-N1/2+k1]           * ck02 * ck11;
      f_hat[k0*N1+N1/2+k1]        = g_hat[(n0-N0/2+k0)*n1+k1]         * ck01 * ck12;
      f_hat[(N0/2+k0)*N1+N1/2+k1] = g_hat[k0*n1+k1]                   * ck02 * ck12;
    }
      }
  TOC(0)
}

/* ############################################################ SPECIFIC VERSIONS FOR d=3 */

static void nfft_3d_init_fg_exp_l(R *fg_exp_l, const int m, const R b)
{
  int l;
  R fg_exp_b0, fg_exp_b1, fg_exp_b2, fg_exp_b0_sq;

  fg_exp_b0 = EXP(-1.0/b);
  fg_exp_b0_sq = fg_exp_b0*fg_exp_b0;
  fg_exp_b1 = K(1.0);
  fg_exp_b2 = K(1.0);
  fg_exp_l[0] = K(1.0);
  for(l=1; l <= 2*m+1; l++)
    {
      fg_exp_b2 = fg_exp_b1*fg_exp_b0;
      fg_exp_b1 *= fg_exp_b0_sq;
      fg_exp_l[l] = fg_exp_l[l-1]*fg_exp_b2;
    }
}

static void nfft_trafo_3d_compute(C *fj, const C *g,
          const R *psij_const0, const R *psij_const1, const R *psij_const2,
          const R *xj0, const R *xj1, const R *xj2,
          const int n0, const int n1, const int n2, const int m)
{
  int u0,o0,l0,u1,o1,l1,u2,o2,l2;
  const C *gj;
  const R *psij0,*psij1,*psij2;

  psij0=psij_const0;
  psij1=psij_const1;
  psij2=psij_const2;

  nfft_uo2(&u0,&o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);
  nfft_uo2(&u2,&o2,*xj2, n2, m);

  *fj=0;

  if(u0<o0)
    if(u1<o1)
      if(u2<o2)
  for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
      psij1=psij_const1;
      for(l1=0; l1<=2*m+1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
    for(l2=0; l2<=2*m+1; l2++)
      (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
        }
    }
      else/* asserts (u2>o2)*/
  for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
      psij1=psij_const1;
      for(l1=0; l1<=2*m+1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
    for(l2=0; l2<2*m+1-o2; l2++)
      (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    gj=g+((u0+l0)*n1+(u1+l1))*n2;
    for(l2=0; l2<=o2; l2++)
      (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
        }
    }
    else/* asserts (u1>o1)*/
      if(u2<o2)
  for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
      psij1=psij_const1;
      for(l1=0; l1<2*m+1-o1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
    for(l2=0; l2<=2*m+1; l2++)
      (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
        }
      for(l1=0; l1<=o1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+l1)*n2+u2;
    for(l2=0; l2<=2*m+1; l2++)
      (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
        }
    }
      else/* asserts (u2>o2) */
  {
    for(l0=0; l0<=2*m+1; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+((u0+l0)*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+l1)*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+((u0+l0)*n1+l1)*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }
  }
  else/* asserts (u0>o0) */
    if(u1<o1)
      if(u2<o2)
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }

    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }
  }
      else/* asserts (u2>o2) */
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+((u0+l0)*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }

    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+(l0*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }
  }
    else/* asserts (u1>o1) */
      if(u2<o2)
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+l1)*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }
    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+l1)*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }
  }
      else/* asserts (u2>o2) */
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+((u0+l0)*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+l1)*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+((u0+l0)*n1+l1)*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }

    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+(l0*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+l1)*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
      gj=g+(l0*n1+l1)*n2;
      for(l2=0; l2<=o2; l2++)
        (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
    }
      }
  }
}

#ifdef _OPENMP
/** 
 * Adjoint NFFT three-dimensional case updating only a specified range of
 * vector g.
 *
 * \arg fg input coefficient f[j]
 * \arg g output vector g
 * \arg psij_const0 vector of window function values first component
 * \arg psij_const1 vector of window function values second component
 * \arg psij_const2 vector of window function values third component
 * \arg xj0 node x[3*j]
 * \arg xj1 node x[3*j+1]
 * \arg xj2 node x[3*j+2]
 * \arg n0 FFTW length (number oversampled Fourier coefficients) first comp.
 * \arg n1 FFTW length (number oversampled Fourier coefficients) second comp.
 * \arg n2 FFTW length (number oversampled Fourier coefficients) third comp.
 * \arg m window length
 * \arg my_u0 lowest index (first component) the current thread writes to
 * \arg my_o0 highest index (second component) the current thread writes to
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_3d_compute_omp3(const C *fj, C *g,
            const R *psij_const0, const R *psij_const1, const R *psij_const2,
            const R *xj0, const R *xj1, const R *xj2,
            const int n0, const int n1, const int n2, const int m,
	    const int my_u0, const int my_o0)
{
  int ar_u0,ar_o0,l0,u1,o1,l1,u2,o2,l2;
  const int lprod = (2*m+2) * (2*m+2) * (2*m+2);

  unsigned long int index_temp1[2*m+2];
  unsigned long int index_temp2[2*m+2];

  nfft_uo2(&ar_u0,&ar_o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);
  nfft_uo2(&u2,&o2,*xj2, n2, m);

  for (l1=0; l1<=2*m+1; l1++)
    index_temp1[l1] = (u1+l1)%n1;

  for (l2=0; l2<=2*m+1; l2++)
    index_temp2[l2] = (u2+l2)%n2;

  if(ar_u0<ar_o0)
  {
    int u0 = MAX(my_u0,ar_u0);
    int o0 = MIN(my_o0,ar_o0);
    int offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      const unsigned long int i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
      {
        const unsigned long int i1 = (i0 + index_temp1[l1]) * n2;
        const C val1 = psij_const1[l1];

        for(l2=0; l2<=2*m+1; l2++)
          g[i1 + index_temp2[l2]] += val0 * val1 * psij_const2[l2] * (*fj);
      }
    }  
  }
  else
  {
    int u0 = MAX(my_u0,ar_u0);
    int o0 = my_o0;
    int offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      unsigned long int i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
      {
        const unsigned long int i1 = (i0 + index_temp1[l1]) * n2;
        const C val1 = psij_const1[l1];

        for(l2=0; l2<=2*m+1; l2++)
          g[i1 + index_temp2[l2]] += val0 * val1 * psij_const2[l2] * (*fj);
      }
    }

    u0 = my_u0;
    o0 = MIN(my_o0,ar_o0);
    offset_psij += my_u0-ar_u0+n0;
#ifdef OMP_ASSERT
if (u0<=o0)
{
  assert(o0-u0 <= 2*m+1);
  assert(offset_psij+o0-u0 <= 2*m+1);
}
#endif
    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      unsigned long int i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
      {
        const unsigned long int i1 = (i0 + index_temp1[l1]) * n2;
        const C val1 = psij_const1[l1];

        for(l2=0; l2<=2*m+1; l2++)
          g[i1 + index_temp2[l2]] += val0 * val1 * psij_const2[l2] * (*fj);
      }
    }
  }
}
#endif

#ifdef _OPENMP
/* adjoint NFFT three-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_3d_compute_omp(const C *fj, C *g,
            const R *psij_const0, const R *psij_const1, const R *psij_const2,
            const R *xj0, const R *xj1, const R *xj2,
            const int n0, const int n1, const int n2, const int m)
{
  int u0,o0,l0,u1,o1,l1,u2,o2,l2;
  const int lprod = (2*m+2) * (2*m+2) * (2*m+2);

  unsigned long int index_temp0[2*m+2];
  unsigned long int index_temp1[2*m+2];
  unsigned long int index_temp2[2*m+2];

  nfft_uo2(&u0,&o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);
  nfft_uo2(&u2,&o2,*xj2, n2, m);

  for (l0=0; l0<=2*m+1; l0++)
    index_temp0[l0] = (u0+l0)%n0;

  for (l1=0; l1<=2*m+1; l1++)
    index_temp1[l1] = (u1+l1)%n1;

  for (l2=0; l2<=2*m+1; l2++)
    index_temp2[l2] = (u2+l2)%n2;

  for(l0=0; l0<=2*m+1; l0++)
  {
    for(l1=0; l1<=2*m+1; l1++)
    {
      for(l2=0; l2<=2*m+1; l2++)
      {
        unsigned long int i = (index_temp0[l0] * n1 + index_temp1[l1]) * n2 + index_temp2[l2];
        C *lhs = g+i;
        R *lhs_real = (R*)lhs;
        C val = psij_const0[l0] * psij_const1[l1] * psij_const2[l2] * (*fj);

        #pragma omp atomic
        lhs_real[0] += creal(val);

        #pragma omp atomic
        lhs_real[1] += cimag(val);
      }
    }
  }
}
#endif

static void nfft_adjoint_3d_compute(const C *fj, C *g,
            const R *psij_const0, const R *psij_const1, const R *psij_const2,
            const R *xj0, const R *xj1, const R *xj2,
            const int n0, const int n1, const int n2, const int m)
{
  int u0,o0,l0,u1,o1,l1,u2,o2,l2;
  C *gj;
  const R *psij0,*psij1,*psij2;

  psij0=psij_const0;
  psij1=psij_const1;
  psij2=psij_const2;

  nfft_uo2(&u0,&o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);
  nfft_uo2(&u2,&o2,*xj2, n2, m);

  if(u0<o0)
    if(u1<o1)
      if(u2<o2)
  for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
      psij1=psij_const1;
      for(l1=0; l1<=2*m+1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
    for(l2=0; l2<=2*m+1; l2++)
      (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
        }
    }
      else/* asserts (u2>o2)*/
  for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
      psij1=psij_const1;
      for(l1=0; l1<=2*m+1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
    for(l2=0; l2<2*m+1-o2; l2++)
      (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    gj=g+((u0+l0)*n1+(u1+l1))*n2;
    for(l2=0; l2<=o2; l2++)
      (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
        }
    }
    else/* asserts (u1>o1)*/
      if(u2<o2)
  for(l0=0; l0<=2*m+1; l0++,psij0++)
    {
      psij1=psij_const1;
      for(l1=0; l1<2*m+1-o1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
    for(l2=0; l2<=2*m+1; l2++)
      (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
        }
      for(l1=0; l1<=o1; l1++,psij1++)
        {
    psij2=psij_const2;
    gj=g+((u0+l0)*n1+l1)*n2+u2;
    for(l2=0; l2<=2*m+1; l2++)
      (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
        }
    }
      else/* asserts (u2>o2) */
  {
    for(l0=0; l0<=2*m+1; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+((u0+l0)*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+l1)*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+((u0+l0)*n1+l1)*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }
  }
  else/* asserts (u0>o0) */
    if(u1<o1)
      if(u2<o2)
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }

    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }
  }
      else/* asserts (u2>o2) */
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+((u0+l0)*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }

    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<=2*m+1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+(l0*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }
  }
    else/* asserts (u1>o1) */
      if(u2<o2)
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+l1)*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }
    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+l1)*n2+u2;
      for(l2=0; l2<=2*m+1; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }
  }
      else/* asserts (u2>o2) */
  {
    for(l0=0; l0<2*m+1-o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+((u0+l0)*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+((u0+l0)*n1+l1)*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+((u0+l0)*n1+l1)*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }

    for(l0=0; l0<=o0; l0++,psij0++)
      {
        psij1=psij_const1;
        for(l1=0; l1<2*m+1-o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+(u1+l1))*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+(l0*n1+(u1+l1))*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
        for(l1=0; l1<=o1; l1++,psij1++)
    {
      psij2=psij_const2;
      gj=g+(l0*n1+l1)*n2+u2;
      for(l2=0; l2<2*m+1-o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
      gj=g+(l0*n1+l1)*n2;
      for(l2=0; l2<=o2; l2++)
        (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
    }
      }
  }
}


static void nfft_trafo_3d_B(nfft_plan *ths)
{
  const int N0 = ths->N[0];
  const int n0 = ths->n[0];
  const int N1 = ths->N[1];
  const int n1 = ths->n[1];
  const int N2 = ths->N[2];
  const int n2 = ths->n[2];
  const int M = ths->M_total;
  const int m = ths->m;

  const C* g = (C*) ths->g;

  int k;

  if(ths->nfft_flags & PRE_FULL_PSI)
  {
    const int lprod = (2*m+2) * (2*m+2) * (2*m+2);
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int l;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < lprod; l++)
        ths->f[j] += ths->psi[j*lprod+l] * g[ths->psi_index_g[j*lprod+l]];
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
  {
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      nfft_trafo_3d_compute(ths->f+j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      int l;
      R psij_const[3*(2*m+2)];
      R fg_psij0 = ths->psi[2*j*3];
      R fg_psij1 = ths->psi[2*j*3+1];
      R fg_psij2 = K(1.0);

      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      fg_psij0 = ths->psi[2*(j*3+1)];
      fg_psij1 = ths->psi[2*(j*3+1)+1];
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      fg_psij0 = ths->psi[2*(j*3+2)];
      fg_psij1 = ths->psi[2*(j*3+2)+1];
      fg_psij2 = K(1.0);
      psij_const[2*(2*m+2)] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
      }

      nfft_trafo_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      int u, o, l;
      R psij_const[3*(2*m+2)];
      R fg_psij0, fg_psij1, fg_psij2;

      nfft_uo(ths,j,&u,&o,0);
      fg_psij0 = (PHI(ths->x[3*j]-((R)u)/n0,0));
      fg_psij1 = EXP(K(2.0)*(n0*(ths->x[3*j]) - u)/ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      nfft_uo(ths,j,&u,&o,1);
      fg_psij0 = (PHI(ths->x[3*j+1]-((R)u)/n1,1));
      fg_psij1 = EXP(K(2.0)*(n1*(ths->x[3*j+1]) - u)/ths->b[1]);
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      nfft_uo(ths,j,&u,&o,2);
      fg_psij0 = (PHI(ths->x[3*j+2]-((R)u)/n2,2));
      fg_psij1 = EXP(K(2.0)*(n2*(ths->x[3*j+2]) - u)/ths->b[2]);
      fg_psij2 = K(1.0);
      psij_const[2*(2*m+2)] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
      }

      nfft_trafo_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
    }

    return;
  } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
  {
    const int K = ths->K, ip_s = K / (m + 2);

    nfft_sort_nodes(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u, o, l;
      R ip_y, ip_w;
      int ip_u;
      R psij_const[3*(2*m+2)];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      nfft_uo(ths,j,&u,&o,0);
      ip_y = FABS(n0*ths->x[3*j+0] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_uo(ths,j,&u,&o,1);
      ip_y = FABS(n1*ths->x[3*j+1] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_uo(ths,j,&u,&o,2);
      ip_y = FABS(n2*ths->x[3*j+2] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_trafo_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */

  nfft_sort_nodes(ths);

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < M; k++)
  {
    R psij_const[3*(2*m+2)];
    int u, o, l;
    int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    nfft_uo(ths,j,&u,&o,0);
    for(l=0;l<=2*m+1;l++)
      psij_const[l]=(PHI(ths->x[3*j]-((R)((u+l)))/n0,0));

    nfft_uo(ths,j,&u,&o,1);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*m+2+l]=(PHI(ths->x[3*j+1]-((R)((u+l)))/n1,1));

    nfft_uo(ths,j,&u,&o,2);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*(2*m+2)+l]=(PHI(ths->x[3*j+2]-((R)((u+l)))/n2,2));

    nfft_trafo_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
  }
}

#ifdef _OPENMP
/**
 * Determines the blocks of vector g the current thread is responsible for.
 *
 * \arg my_u0 lowest index (first component) the current threads writes to in g
 * \arg my_o0 highest index (first component) the current threads writes to in g
 * \arg min_u_a lowest (linear 3d) index u which could lead to writing to g
 * \arg max_u_a highest (linear 3d) index o which could lead to writing to g
 * \arg min_u_b lowest (linear 3d) index u which could lead to writing to g
 * \arg max_u_b highest (linear 3d) index o which could lead to writing to g
 * \arg n0 FFTW length first component
 * \arg n1 FFTW length second component
 * \arg n2 FFTW length third component
 * \arg m window length
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_3d_B_omp3_init_(int *my_u0, int *my_o0, int *min_u_a, int *max_u_a, int *min_u_b, int *max_u_b, const int n0, const int n1, const int n2, const int m)
{
  int k;
  int nthreads = omp_get_num_threads();
  int nthreads_used = MIN(nthreads, n0);
  int size_per_thread = n0 / nthreads_used;
  int size_left = n0 - size_per_thread * nthreads_used;
  int size_g[nthreads_used];
  int offset_g[nthreads_used];
  int my_id = omp_get_thread_num();

  *min_u_a = -1;
  *max_u_a = -1;
  *min_u_b = -1;
  *max_u_b = -1;
  *my_u0 = -1;
  *my_o0 = -1;

  if (my_id < nthreads_used)
  {
    const int m22 = 2 * m + 2;

    offset_g[0] = 0;
    for (k = 0; k < nthreads_used; k++)
    {
      if (k > 0)
        offset_g[k] = offset_g[k-1] + size_g[k-1];
      size_g[k] = size_per_thread;
      if (size_left > 0)
      {
        size_g[k]++;
	size_left--;
      }
    }

    *my_u0 = offset_g[my_id];
    *my_o0 = offset_g[my_id] + size_g[my_id] - 1;

    if (nthreads_used > 1)
    {
      *max_u_a = n2*n1*(offset_g[my_id] + size_g[my_id]) - 1;
      *min_u_a = n2*n1*(offset_g[my_id] - m22 + 1);
    }
    else
    {
      *min_u_a = 0;
      *max_u_a = n2 * n1 * n0 - 1;
    }

    if (*min_u_a < 0)
    {
      *min_u_b = n2 * n1 * (offset_g[my_id] - m22 + 1 + n0);
      *max_u_b = n2 * n1 * n0 - 1;
      *min_u_a = 0;
    }

    if (*min_u_b != -1 && *min_u_b <= *max_u_a)
    {
      *max_u_a = *max_u_b;
      *min_u_b = -1;
      *max_u_b = -1;
    }
    assert(*min_u_a <= *max_u_a);
    assert(*min_u_b <= *max_u_b);
    assert(*min_u_b == -1 || *max_u_a < *min_u_b);
  }
}
#endif

static void nfft_adjoint_3d_B(nfft_plan *ths)
{
  int k;
  const int N0 = ths->N[0];
  const int n0 = ths->n[0];
  const int N1 = ths->N[1];
  const int n1 = ths->n[1];
  const int N2 = ths->N[2];
  const int n2 = ths->n[2];
  const int M = ths->M_total;
  const int m = ths->m;

  C* g = (C*) ths->g;

  memset(g,0,ths->n_total*sizeof(C));

  if(ths->nfft_flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, ths->f, M,
        3, ths->n, m, ths->nfft_flags, ths->index_x);
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
  {
#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 3, ths->n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          int k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp(ths->f+j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute(ths->f+j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 3, ths->n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R psij_const[3*(2*m+2)];
            R fg_psij0 = ths->psi[2*j*3];
            R fg_psij1 = ths->psi[2*j*3+1];
            R fg_psij2 = K(1.0);

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            psij_const[0] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
            }

            fg_psij0 = ths->psi[2*(j*3+1)];
            fg_psij1 = ths->psi[2*(j*3+1)+1];
            fg_psij2 = K(1.0);
            psij_const[2*m+2] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
            }

            fg_psij0 = ths->psi[2*(j*3+2)];
            fg_psij1 = ths->psi[2*(j*3+2)+1];
            fg_psij2 = K(1.0);
            psij_const[2*(2*m+2)] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
            }

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          int k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R psij_const[3*(2*m+2)];
            R fg_psij0 = ths->psi[2*j*3];
            R fg_psij1 = ths->psi[2*j*3+1];
            R fg_psij2 = K(1.0);

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            psij_const[0] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
            }

            fg_psij0 = ths->psi[2*(j*3+1)];
            fg_psij1 = ths->psi[2*(j*3+1)+1];
            fg_psij2 = K(1.0);
            psij_const[2*m+2] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
            }

            fg_psij0 = ths->psi[2*(j*3+2)];
            fg_psij1 = ths->psi[2*(j*3+2)+1];
            fg_psij2 = K(1.0);
            psij_const[2*(2*m+2)] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
            }

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      R psij_const[3*(2*m+2)];
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      int l;
      R fg_psij0 = ths->psi[2*j*3];
      R fg_psij1 = ths->psi[2*j*3+1];
      R fg_psij2 = K(1.0);

      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      fg_psij0 = ths->psi[2*(j*3+1)];
      fg_psij1 = ths->psi[2*(j*3+1)+1];
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      fg_psij0 = ths->psi[2*(j*3+2)];
      fg_psij1 = ths->psi[2*(j*3+2)+1];
      fg_psij2 = K(1.0);
      psij_const[2*(2*m+2)] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
      }

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

    nfft_sort_nodes(ths);

#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 3, ths->n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R psij_const[3*(2*m+2)];
            R fg_psij0, fg_psij1, fg_psij2;

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_uo(ths,j,&u,&o,0);
            fg_psij0 = (PHI(ths->x[3*j]-((R)u)/n0,0));
            fg_psij1 = EXP(K(2.0)*(n0*(ths->x[3*j]) - u)/ths->b[0]);
            fg_psij2 = K(1.0);
            psij_const[0] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
            }

            nfft_uo(ths,j,&u,&o,1);
            fg_psij0 = (PHI(ths->x[3*j+1]-((R)u)/n1,1));
            fg_psij1 = EXP(K(2.0)*(n1*(ths->x[3*j+1]) - u)/ths->b[1]);
            fg_psij2 = K(1.0);
            psij_const[2*m+2] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
            }

            nfft_uo(ths,j,&u,&o,2);
            fg_psij0 = (PHI(ths->x[3*j+2]-((R)u)/n2,2));
            fg_psij1 = EXP(K(2.0)*(n2*(ths->x[3*j+2]) - u)/ths->b[2]);
            fg_psij2 = K(1.0);
            psij_const[2*(2*m+2)] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
            }

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          int k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R psij_const[3*(2*m+2)];
            R fg_psij0, fg_psij1, fg_psij2;

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_uo(ths,j,&u,&o,0);
            fg_psij0 = (PHI(ths->x[3*j]-((R)u)/n0,0));
            fg_psij1 = EXP(K(2.0)*(n0*(ths->x[3*j]) - u)/ths->b[0]);
            fg_psij2 = K(1.0);
            psij_const[0] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
            }

            nfft_uo(ths,j,&u,&o,1);
            fg_psij0 = (PHI(ths->x[3*j+1]-((R)u)/n1,1));
            fg_psij1 = EXP(K(2.0)*(n1*(ths->x[3*j+1]) - u)/ths->b[1]);
            fg_psij2 = K(1.0);
            psij_const[2*m+2] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
            }

            nfft_uo(ths,j,&u,&o,2);
            fg_psij0 = (PHI(ths->x[3*j+2]-((R)u)/n2,2));
            fg_psij1 = EXP(K(2.0)*(n2*(ths->x[3*j+2]) - u)/ths->b[2]);
            fg_psij2 = K(1.0);
            psij_const[2*(2*m+2)] = fg_psij0;
            for(l=1; l<=2*m+1; l++)
            {
              fg_psij2 *= fg_psij1;
              psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
            }

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma openmp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u,o,l;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[3*(2*m+2)];
      R fg_psij0, fg_psij1, fg_psij2;

      nfft_uo(ths,j,&u,&o,0);
      fg_psij0 = (PHI(ths->x[3*j]-((R)u)/n0,0));
      fg_psij1 = EXP(K(2.0)*(n0*(ths->x[3*j]) - u)/ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      nfft_uo(ths,j,&u,&o,1);
      fg_psij0 = (PHI(ths->x[3*j+1]-((R)u)/n1,1));
      fg_psij1 = EXP(K(2.0)*(n1*(ths->x[3*j+1]) - u)/ths->b[1]);
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      nfft_uo(ths,j,&u,&o,2);
      fg_psij0 = (PHI(ths->x[3*j+2]-((R)u)/n2,2));
      fg_psij1 = EXP(K(2.0)*(n2*(ths->x[3*j+2]) - u)/ths->b[2]);
      fg_psij2 = K(1.0);
      psij_const[2*(2*m+2)] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
      }

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
  {
    const int K = ths->K;
    const int ip_s = K / (m + 2);

    nfft_sort_nodes(ths);

#ifdef _OPENMP
    if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    {
      #pragma omp parallel private(k)
      {
        int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
        int *ar_x = ths->index_x;
        int ip_u;
        R ip_y, ip_w;
        R psij_const[3*(2*m+2)];

        nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 3, ths->n, m);

        if (min_u_a != -1)
        {
          k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_a || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_a);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R psij_const[3*(2*m+2)];

            if (u_prod < min_u_a || u_prod > max_u_a)
              break;

            nfft_uo(ths,j,&u,&o,0);
            ip_y = FABS(n0*ths->x[3*j+0] - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_uo(ths,j,&u,&o,1);
            ip_y = FABS(n1*ths->x[3*j+1] - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_uo(ths,j,&u,&o,2);
            ip_y = FABS(n2*ths->x[3*j+2] - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }

        if (min_u_b != -1)
        {
          int k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
          assert(ar_x[2*k] >= min_u_b || k == M-1);
          if (k > 0)
            assert(ar_x[2*k-2] < min_u_b);
#endif
          while (k < M)
          {
            int u, o, l;
            int u_prod = ar_x[2*k];
            int j = ar_x[2*k+1];
            R psij_const[3*(2*m+2)];

            if (u_prod < min_u_b || u_prod > max_u_b)
              break;

            nfft_uo(ths,j,&u,&o,0);
            ip_y = FABS(n0*ths->x[3*j+0] - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_uo(ths,j,&u,&o,1);
            ip_y = FABS(n1*ths->x[3*j+1] - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_uo(ths,j,&u,&o,2);
            ip_y = FABS(n2*ths->x[3*j+2] - u)*((R)ip_s);
            ip_u = LRINT(FLOOR(ip_y));
            ip_w = ip_y-ip_u;
            for(l=0; l < 2*m+2; l++)
              psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
                ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

            nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

            k++;
          }
        }
      } /* omp parallel */
      return;
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

    #pragma openmp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      int u,o,l;
      int ip_u;
      R ip_y, ip_w;
      int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[3*(2*m+2)];

      nfft_uo(ths,j,&u,&o,0);
      ip_y = FABS(n0*ths->x[3*j+0] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_uo(ths,j,&u,&o,1);
      ip_y = FABS(n1*ths->x[3*j+1] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_uo(ths,j,&u,&o,2);
      ip_y = FABS(n2*ths->x[3*j+2] - u)*((R)ip_s);
      ip_u = LRINT(FLOOR(ip_y));
      ip_w = ip_y-ip_u;
      for(l=0; l < 2*m+2; l++)
        psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  nfft_sort_nodes(ths);

#ifdef _OPENMP
  if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel private(k)
    {
      int my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
      int *ar_x = ths->index_x;
      R psij_const[3*(2*m+2)];

      nfft_adjoint_B_omp3_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, 3, ths->n, m);

      if (min_u_a != -1)
      {
        k = index_x_binary_search(ar_x, M, min_u_a);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_a || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_a);
#endif
        while (k < M)
        {
	  int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];
          R psij_const[3*(2*m+2)];

          if (u_prod < min_u_a || u_prod > max_u_a)
            break;

          nfft_uo(ths,j,&u,&o,0);
          for(l=0;l<=2*m+1;l++)
            psij_const[l]=(PHI(ths->x[3*j]-((R)((u+l)))/n0,0));

          nfft_uo(ths,j,&u,&o,1);
          for(l=0;l<=2*m+1;l++)
            psij_const[2*m+2+l]=(PHI(ths->x[3*j+1]-((R)((u+l)))/n1,1));

          nfft_uo(ths,j,&u,&o,2);
          for(l=0;l<=2*m+1;l++)
            psij_const[2*(2*m+2)+l]=(PHI(ths->x[3*j+2]-((R)((u+l)))/n2,2));

          nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

          k++;
        }
      }

      if (min_u_b != -1)
      {
        int k = index_x_binary_search(ar_x, M, min_u_b);
#ifdef OMP_ASSERT
        assert(ar_x[2*k] >= min_u_b || k == M-1);
        if (k > 0)
          assert(ar_x[2*k-2] < min_u_b);
#endif
        while (k < M)
        {
	  int u, o, l;
          int u_prod = ar_x[2*k];
          int j = ar_x[2*k+1];
          R psij_const[3*(2*m+2)];

          if (u_prod < min_u_b || u_prod > max_u_b)
            break;

          nfft_uo(ths,j,&u,&o,0);
          for(l=0;l<=2*m+1;l++)
            psij_const[l]=(PHI(ths->x[3*j]-((R)((u+l)))/n0,0));

          nfft_uo(ths,j,&u,&o,1);
          for(l=0;l<=2*m+1;l++)
            psij_const[2*m+2+l]=(PHI(ths->x[3*j+1]-((R)((u+l)))/n1,1));

          nfft_uo(ths,j,&u,&o,2);
          for(l=0;l<=2*m+1;l++)
            psij_const[2*(2*m+2)+l]=(PHI(ths->x[3*j+2]-((R)((u+l)))/n2,2));

          nfft_adjoint_3d_compute_omp3(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m, my_u0, my_o0);

          k++;
        }
      }
    } /* omp parallel */
    return;
  } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */
#endif

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < M; k++)
  {
    int u,o,l;
    R psij_const[3*(2*m+2)];
    int j = (ths->nfft_flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    nfft_uo(ths,j,&u,&o,0);
    for(l=0;l<=2*m+1;l++)
      psij_const[l]=(PHI(ths->x[3*j]-((R)((u+l)))/n0,0));

    nfft_uo(ths,j,&u,&o,1);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*m+2+l]=(PHI(ths->x[3*j+1]-((R)((u+l)))/n1,1));

    nfft_uo(ths,j,&u,&o,2);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*(2*m+2)+l]=(PHI(ths->x[3*j+2]-((R)((u+l)))/n2,2));

#ifdef _OPENMP
    nfft_adjoint_3d_compute_omp(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
    nfft_adjoint_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
  }
}


void nfft_trafo_3d(nfft_plan *ths)
{
  int k0,k1,k2,n0,n1,n2,N0,N1,N2;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12, *c_phi_inv21, *c_phi_inv22;
  R ck01, ck02, ck11, ck12, ck21, ck22;
  C *g_hat111,*f_hat111,*g_hat211,*f_hat211,*g_hat121,*f_hat121,*g_hat221,*f_hat221;
  C *g_hat112,*f_hat112,*g_hat212,*f_hat212,*g_hat122,*f_hat122,*g_hat222,*f_hat222;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  N2=ths->N[2];
  n0=ths->n[0];
  n1=ths->n[1];
  n2=ths->n[2];

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(0)
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k0)
  for (k0 = 0; k0 < ths->n_total; k0++)
    ths->g_hat[k0] = 0.0;
#else
  memset(ths->g_hat,0,ths->n_total*sizeof(C));
#endif

  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,c_phi_inv11,c_phi_inv12,ck11,ck12,c_phi_inv21,c_phi_inv22,g_hat111,f_hat111,g_hat211,f_hat211,g_hat121,f_hat121,g_hat221,f_hat221,g_hat112,f_hat112,g_hat212,f_hat212,g_hat122,f_hat122,g_hat222,f_hat222,ck21,ck22)
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=c_phi_inv01[k0];
    ck02=c_phi_inv02[k0];
    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];

    for(k1=0;k1<N1/2;k1++)
      {
        ck11=c_phi_inv11[k1];
        ck12=c_phi_inv12[k1];
        c_phi_inv21=ths->c_phi_inv[2];
        c_phi_inv22=&ths->c_phi_inv[2][N2/2];

        g_hat111=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat111=f_hat + (k0*N1+k1)*N2;
        g_hat211=g_hat + (k0*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat211=f_hat + (((N0/2)+k0)*N1+k1)*N2;
        g_hat121=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2+n2-(N2/2);
        f_hat121=f_hat + (k0*N1+(N1/2)+k1)*N2;
        g_hat221=g_hat + (k0*n1+k1)*n2+n2-(N2/2);
        f_hat221=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2;

        g_hat112=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2;
        f_hat112=f_hat + (k0*N1+k1)*N2+(N2/2);
        g_hat212=g_hat + (k0*n1+n1-(N1/2)+k1)*n2;
        f_hat212=f_hat + (((N0/2)+k0)*N1+k1)*N2+(N2/2);
        g_hat122=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2;
        f_hat122=f_hat + (k0*N1+N1/2+k1)*N2+(N2/2);
        g_hat222=g_hat + (k0*n1+k1)*n2;
        f_hat222=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2+(N2/2);

        for(k2=0;k2<N2/2;k2++)
    {
      ck21=c_phi_inv21[k2];
      ck22=c_phi_inv22[k2];

      g_hat111[k2] = f_hat111[k2] * ck01 * ck11 * ck21;
      g_hat211[k2] = f_hat211[k2] * ck02 * ck11 * ck21;
      g_hat121[k2] = f_hat121[k2] * ck01 * ck12 * ck21;
      g_hat221[k2] = f_hat221[k2] * ck02 * ck12 * ck21;

      g_hat112[k2] = f_hat112[k2] * ck01 * ck11 * ck22;
      g_hat212[k2] = f_hat212[k2] * ck02 * ck11 * ck22;
      g_hat122[k2] = f_hat122[k2] * ck01 * ck12 * ck22;
      g_hat222[k2] = f_hat222[k2] * ck02 * ck12 * ck22;
    }
      }
  }
#else
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=(*c_phi_inv01++);
    ck02=(*c_phi_inv02++);
    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];

    for(k1=0;k1<N1/2;k1++)
      {
        ck11=(*c_phi_inv11++);
        ck12=(*c_phi_inv12++);
        c_phi_inv21=ths->c_phi_inv[2];
        c_phi_inv22=&ths->c_phi_inv[2][N2/2];

        g_hat111=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat111=f_hat + (k0*N1+k1)*N2;
        g_hat211=g_hat + (k0*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat211=f_hat + (((N0/2)+k0)*N1+k1)*N2;
        g_hat121=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2+n2-(N2/2);
        f_hat121=f_hat + (k0*N1+(N1/2)+k1)*N2;
        g_hat221=g_hat + (k0*n1+k1)*n2+n2-(N2/2);
        f_hat221=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2;

        g_hat112=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2;
        f_hat112=f_hat + (k0*N1+k1)*N2+(N2/2);
        g_hat212=g_hat + (k0*n1+n1-(N1/2)+k1)*n2;
        f_hat212=f_hat + (((N0/2)+k0)*N1+k1)*N2+(N2/2);
        g_hat122=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2;
        f_hat122=f_hat + (k0*N1+N1/2+k1)*N2+(N2/2);
        g_hat222=g_hat + (k0*n1+k1)*n2;
        f_hat222=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2+(N2/2);

        for(k2=0;k2<N2/2;k2++)
    {
      ck21=(*c_phi_inv21++);
      ck22=(*c_phi_inv22++);

      (*g_hat111++) = (*f_hat111++) * ck01 * ck11 * ck21;
      (*g_hat211++) = (*f_hat211++) * ck02 * ck11 * ck21;
      (*g_hat121++) = (*f_hat121++) * ck01 * ck12 * ck21;
      (*g_hat221++) = (*f_hat221++) * ck02 * ck12 * ck21;

      (*g_hat112++) = (*f_hat112++) * ck01 * ck11 * ck22;
      (*g_hat212++) = (*f_hat212++) * ck02 * ck11 * ck22;
      (*g_hat122++) = (*f_hat122++) * ck01 * ck12 * ck22;
      (*g_hat222++) = (*f_hat222++) * ck02 * ck12 * ck22;
    }
      }
  }
#endif
    }
  else
    #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,ck11,ck12,ck21,ck22)
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(k1,1));

      for(k2=0;k2<N2/2;k2++)
        {
    ck21=K(1.0)/(PHI_HUT(k2-N2/2,2));
    ck22=K(1.0)/(PHI_HUT(k2,2));

    g_hat[((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2+n2-N2/2+k2] = f_hat[(k0*N1+k1)*N2+k2]                  * ck01 * ck11 * ck21;
    g_hat[(k0*n1+n1-N1/2+k1)*n2+n2-N2/2+k2]           = f_hat[((N0/2+k0)*N1+k1)*N2+k2]           * ck02 * ck11 * ck21;
    g_hat[((n0-N0/2+k0)*n1+k1)*n2+n2-N2/2+k2]         = f_hat[(k0*N1+N1/2+k1)*N2+k2]             * ck01 * ck12 * ck21;
    g_hat[(k0*n1+k1)*n2+n2-N2/2+k2]                   = f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+k2]      * ck02 * ck12 * ck21;

    g_hat[((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2+k2]         = f_hat[(k0*N1+k1)*N2+N2/2+k2]             * ck01 * ck11 * ck22;
    g_hat[(k0*n1+n1-N1/2+k1)*n2+k2]                   = f_hat[((N0/2+k0)*N1+k1)*N2+N2/2+k2]      * ck02 * ck11 * ck22;
    g_hat[((n0-N0/2+k0)*n1+k1)*n2+k2]                 = f_hat[(k0*N1+N1/2+k1)*N2+N2/2+k2]        * ck01 * ck12 * ck22;
    g_hat[(k0*n1+k1)*n2+k2]                           = f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+N2/2+k2] * ck02 * ck12 * ck22;
        }
    }
      }

  TOC(0)

  TIC_FFTW(1)
  fftw_execute(ths->my_fftw_plan1);
  TOC_FFTW(1);

  TIC(2);
  nfft_trafo_3d_B(ths);
  TOC(2);
}

void nfft_adjoint_3d(nfft_plan *ths)
{
  int k0,k1,k2,n0,n1,n2,N0,N1,N2;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12, *c_phi_inv21, *c_phi_inv22;
  R ck01, ck02, ck11, ck12, ck21, ck22;
  C *g_hat111,*f_hat111,*g_hat211,*f_hat211,*g_hat121,*f_hat121,*g_hat221,*f_hat221;
  C *g_hat112,*f_hat112,*g_hat212,*f_hat212,*g_hat122,*f_hat122,*g_hat222,*f_hat222;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  N2=ths->N[2];
  n0=ths->n[0];
  n1=ths->n[1];
  n2=ths->n[2];

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(2);
  nfft_adjoint_3d_B(ths);
  TOC(2);

  TIC_FFTW(1)
  fftw_execute(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,c_phi_inv11,c_phi_inv12,ck11,ck12,c_phi_inv21,c_phi_inv22,g_hat111,f_hat111,g_hat211,f_hat211,g_hat121,f_hat121,g_hat221,f_hat221,g_hat112,f_hat112,g_hat212,f_hat212,g_hat122,f_hat122,g_hat222,f_hat222,ck21,ck22)
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=c_phi_inv01[k0];
    ck02=c_phi_inv02[k0];
    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];

    for(k1=0;k1<N1/2;k1++)
      {
        ck11=c_phi_inv11[k1];
        ck12=c_phi_inv12[k1];
        c_phi_inv21=ths->c_phi_inv[2];
        c_phi_inv22=&ths->c_phi_inv[2][N2/2];

        g_hat111=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat111=f_hat + (k0*N1+k1)*N2;
        g_hat211=g_hat + (k0*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat211=f_hat + (((N0/2)+k0)*N1+k1)*N2;
        g_hat121=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2+n2-(N2/2);
        f_hat121=f_hat + (k0*N1+(N1/2)+k1)*N2;
        g_hat221=g_hat + (k0*n1+k1)*n2+n2-(N2/2);
        f_hat221=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2;

        g_hat112=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2;
        f_hat112=f_hat + (k0*N1+k1)*N2+(N2/2);
        g_hat212=g_hat + (k0*n1+n1-(N1/2)+k1)*n2;
        f_hat212=f_hat + (((N0/2)+k0)*N1+k1)*N2+(N2/2);
        g_hat122=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2;
        f_hat122=f_hat + (k0*N1+(N1/2)+k1)*N2+(N2/2);
        g_hat222=g_hat + (k0*n1+k1)*n2;
        f_hat222=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2+(N2/2);

        for(k2=0;k2<N2/2;k2++)
    {
      ck21=c_phi_inv21[k2];
      ck22=c_phi_inv22[k2];

      f_hat111[k2] = g_hat111[k2] * ck01 * ck11 * ck21;
      f_hat211[k2] = g_hat211[k2] * ck02 * ck11 * ck21;
      f_hat121[k2] = g_hat121[k2] * ck01 * ck12 * ck21;
      f_hat221[k2] = g_hat221[k2] * ck02 * ck12 * ck21;

      f_hat112[k2] = g_hat112[k2] * ck01 * ck11 * ck22;
      f_hat212[k2] = g_hat212[k2] * ck02 * ck11 * ck22;
      f_hat122[k2] = g_hat122[k2] * ck01 * ck12 * ck22;
      f_hat222[k2] = g_hat222[k2] * ck02 * ck12 * ck22;
    }
      }
  }
#else
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=(*c_phi_inv01++);
    ck02=(*c_phi_inv02++);
    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];

    for(k1=0;k1<N1/2;k1++)
      {
        ck11=(*c_phi_inv11++);
        ck12=(*c_phi_inv12++);
        c_phi_inv21=ths->c_phi_inv[2];
        c_phi_inv22=&ths->c_phi_inv[2][N2/2];

        g_hat111=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat111=f_hat + (k0*N1+k1)*N2;
        g_hat211=g_hat + (k0*n1+n1-(N1/2)+k1)*n2+n2-(N2/2);
        f_hat211=f_hat + (((N0/2)+k0)*N1+k1)*N2;
        g_hat121=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2+n2-(N2/2);
        f_hat121=f_hat + (k0*N1+(N1/2)+k1)*N2;
        g_hat221=g_hat + (k0*n1+k1)*n2+n2-(N2/2);
        f_hat221=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2;

        g_hat112=g_hat + ((n0-(N0/2)+k0)*n1+n1-(N1/2)+k1)*n2;
        f_hat112=f_hat + (k0*N1+k1)*N2+(N2/2);
        g_hat212=g_hat + (k0*n1+n1-(N1/2)+k1)*n2;
        f_hat212=f_hat + (((N0/2)+k0)*N1+k1)*N2+(N2/2);
        g_hat122=g_hat + ((n0-(N0/2)+k0)*n1+k1)*n2;
        f_hat122=f_hat + (k0*N1+(N1/2)+k1)*N2+(N2/2);
        g_hat222=g_hat + (k0*n1+k1)*n2;
        f_hat222=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2+(N2/2);

        for(k2=0;k2<N2/2;k2++)
    {
      ck21=(*c_phi_inv21++);
      ck22=(*c_phi_inv22++);

      (*f_hat111++) = (*g_hat111++) * ck01 * ck11 * ck21;
      (*f_hat211++) = (*g_hat211++) * ck02 * ck11 * ck21;
      (*f_hat121++) = (*g_hat121++) * ck01 * ck12 * ck21;
      (*f_hat221++) = (*g_hat221++) * ck02 * ck12 * ck21;

      (*f_hat112++) = (*g_hat112++) * ck01 * ck11 * ck22;
      (*f_hat212++) = (*g_hat212++) * ck02 * ck11 * ck22;
      (*f_hat122++) = (*g_hat122++) * ck01 * ck12 * ck22;
      (*f_hat222++) = (*g_hat222++) * ck02 * ck12 * ck22;
    }
      }
  }
#endif
    }
  else
    #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,ck11,ck12,ck21,ck22)
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(k1,1));

      for(k2=0;k2<N2/2;k2++)
        {
    ck21=K(1.0)/(PHI_HUT(k2-N2/2,2));
    ck22=K(1.0)/(PHI_HUT(k2,2));

    f_hat[(k0*N1+k1)*N2+k2]                  = g_hat[((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2+n2-N2/2+k2] * ck01 * ck11 * ck21;
    f_hat[((N0/2+k0)*N1+k1)*N2+k2]           = g_hat[(k0*n1+n1-N1/2+k1)*n2+n2-N2/2+k2]           * ck02 * ck11 * ck21;
    f_hat[(k0*N1+N1/2+k1)*N2+k2]             = g_hat[((n0-N0/2+k0)*n1+k1)*n2+n2-N2/2+k2]         * ck01 * ck12 * ck21;
    f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+k2]      = g_hat[(k0*n1+k1)*n2+n2-N2/2+k2]                   * ck02 * ck12 * ck21;

    f_hat[(k0*N1+k1)*N2+N2/2+k2]             = g_hat[((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2+k2]         * ck01 * ck11 * ck22;
    f_hat[((N0/2+k0)*N1+k1)*N2+N2/2+k2]      = g_hat[(k0*n1+n1-N1/2+k1)*n2+k2]                   * ck02 * ck11 * ck22;
    f_hat[(k0*N1+N1/2+k1)*N2+N2/2+k2]        = g_hat[((n0-N0/2+k0)*n1+k1)*n2+k2]                 * ck01 * ck12 * ck22;
    f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+N2/2+k2] = g_hat[(k0*n1+k1)*n2+k2]                           * ck02 * ck12 * ck22;
        }
    }
      }

  TOC(0)
}

/** user routines
 */
void nfft_trafo(nfft_plan *ths)
{
  switch(ths->d)
    {
    case 1: nfft_trafo_1d(ths); break;
    case 2: nfft_trafo_2d(ths); break;
    case 3: nfft_trafo_3d(ths); break;
    default:
    /* use ths->my_fftw_plan1 */
    ths->g_hat=ths->g1;
    ths->g=ths->g2;

    /** form \f$ \hat g_k = \frac{\hat f_k}{c_k\left(\phi\right)} \text{ for }
     *  k \in I_N \f$
     */
    TIC(0)
    nfft_D_A(ths);
    TOC(0)

    /** compute by d-variate discrete Fourier transform
     *  \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
     *  \text{ for } l \in I_n \f$
     */
    TIC_FFTW(1)
    fftw_execute(ths->my_fftw_plan1);
    TOC_FFTW(1)

    /** set \f$ f_j =\sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
     *  \text{ for } j=0,\hdots,M_total-1 \f$
     */
    TIC(2)
    nfft_B_A(ths);
    TOC(2)
    }
} /* nfft_trafo */

void nfft_adjoint(nfft_plan *ths)
{
  switch(ths->d)
    {
    case 1: nfft_adjoint_1d(ths); break;
    case 2: nfft_adjoint_2d(ths); break;
    case 3: nfft_adjoint_3d(ths); break;
    default:
      /* use ths->my_fftw_plan2 */
      ths->g_hat=ths->g1;
      ths->g=ths->g2;
      
      /** set \f$ g_l = \sum_{j=0}^{M_total-1} f_j \psi\left(x_j-\frac{l}{n}\right)
       *  \text{ for } l \in I_n,m(x_j) \f$
       */
      TIC(2)
      nfft_B_T(ths);
      TOC(2)

      /** compute by d-variate discrete Fourier transform
       *  \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{+2\pi {\rm i} \frac{kl}{n}}
       *  \text{ for }  k \in I_N\f$
       */
      TIC_FFTW(1)
      fftw_execute(ths->my_fftw_plan2);
      TOC_FFTW(1)

      /** form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
       *  k \in I_N \f$
       */
      TIC(0)
      nfft_D_T(ths);
      TOC(0)
    }
} /* nfft_adjoint */


/** initialisation of direct transform
 */
static void nfft_precompute_phi_hut(nfft_plan *ths)
{
  int ks[ths->d];                       /**< index over all frequencies      */
  int t;                                /**< index over all dimensions       */

  ths->c_phi_inv = (R**) nfft_malloc(ths->d*sizeof(R*));

  for(t=0; t<ths->d; t++)
    {
      ths->c_phi_inv[t]= (R*)nfft_malloc(ths->N[t]*sizeof(R));
      for(ks[t]=0; ks[t]<ths->N[t]; ks[t]++)
  ths->c_phi_inv[t][ks[t]]= K(1.0)/(PHI_HUT(ks[t]-ths->N[t]/2,t));
    }
} /* nfft_phi_hut */

/** create a lookup table, but NOT for each node
 *  good idea K=2^xx
 *  TODO: estimate K, call from init
 *  assumes an EVEN window function
 */
void nfft_precompute_lin_psi(nfft_plan *ths)
{
  int t;                                /**< index over all dimensions       */
  int j;                                /**< index over all nodes            */
  R step;                          /**< step size in [0,(m+2)/n]        */

  for (t=0; t<ths->d; t++)
    {
      step=((R)(ths->m+2))/(((double)ths->K)*ths->n[t]);
      for(j=0;j<=ths->K;j++)
  {
    ths->psi[(ths->K+1)*t + j] = PHI(j*step,t);
  } /* for(j) */
    } /* for(t) */
}

static void nfft_precompute_fg_psi(nfft_plan *ths)
{
  int t;                                /**< index over all dimensions       */
  int u, o;                             /**< depends on x_j                  */

  nfft_sort_nodes(ths);

  for (t=0; t<ths->d; t++)
  {
    int j;
    #pragma omp parallel for default(shared) private(j,u,o)
    for (j = 0; j < ths->M_total; j++)
      {
  nfft_uo(ths,j,&u,&o,t);

        ths->psi[2*(j*ths->d+t)]=
            (PHI((ths->x[j*ths->d+t]-((R)u)/ths->n[t]),t));

        ths->psi[2*(j*ths->d+t)+1]=
            EXP(K(2.0)*(ths->n[t]*ths->x[j*ths->d+t] - u) / ths->b[t]);
      } /* for(j) */
  }
  /* for(t) */
} /* nfft_precompute_fg_psi */

void nfft_precompute_psi(nfft_plan *ths)
{
  int t;                                /**< index over all dimensions       */
  int l;                                /**< index u<=l<=o                   */
  int lj;                               /**< index 0<=lj<u+o+1               */
  int u, o;                             /**< depends on x_j                  */

  nfft_sort_nodes(ths);

  for (t=0; t<ths->d; t++)
  {
    int j;
    #pragma omp parallel for default(shared) private(j,l,lj,u,o)
    for (j = 0; j < ths->M_total; j++)
      {
  nfft_uo(ths,j,&u,&o,t);

  for(l=u, lj=0; l <= o; l++, lj++)
    ths->psi[(j*ths->d+t)*(2*ths->m+2)+lj]=
      (PHI((ths->x[j*ths->d+t]-((R)l)/ths->n[t]),t));
      } /* for(j) */
  }
  /* for(t) */
} /* nfft_precompute_psi */

#ifdef _OPENMP
static void nfft_precompute_full_psi_omp(nfft_plan *ths)
{
  int j;                                /**< index over all nodes            */
  int lprod;                            /**< 'bandwidth' of matrix B         */

  {
    int t;
    for(t=0,lprod = 1; t<ths->d; t++)
        lprod *= 2*ths->m+2;
  }

  #pragma omp parallel for default(shared) private(j)
  for(j=0; j<ths->M_total; j++)
    {
      int t,t2;                             /**< index over all dimensions       */
      int l_L;                              /**< plain index 0<=l_L<lprod        */
      int l[ths->d];                        /**< multi index u<=l<=o             */
      int lj[ths->d];                       /**< multi index 0<=lj<u+o+1         */
      int ll_plain[ths->d+1];               /**< postfix plain index             */

      int u[ths->d], o[ths->d];             /**< depends on x_j                  */

      R phi_prod[ths->d+1];
      int ix = j*lprod;

      phi_prod[0]=1;
      ll_plain[0]=0;

      MACRO_init_uo_l_lj_t;

      for(l_L=0; l_L<lprod; l_L++, ix++)
      {
        MACRO_update_phi_prod_ll_plain(without_PRE_PSI);

        ths->psi_index_g[ix]=ll_plain[ths->d];
        ths->psi[ix]=phi_prod[ths->d];

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */

      ths->psi_index_f[j]=lprod;
    } /* for(j) */
}
#endif

void nfft_precompute_full_psi(nfft_plan *ths)
{
#ifdef _OPENMP
  nfft_sort_nodes(ths);

  nfft_precompute_full_psi_omp(ths);
#else
  int t,t2;                             /**< index over all dimensions       */
  int j;                                /**< index over all nodes            */
  int l_L;                              /**< plain index 0<=l_L<lprod        */
  int l[ths->d];                        /**< multi index u<=l<=o             */
  int lj[ths->d];                       /**< multi index 0<=lj<u+o+1         */
  int ll_plain[ths->d+1];               /**< postfix plain index             */
  int lprod;                            /**< 'bandwidth' of matrix B         */
  int u[ths->d], o[ths->d];             /**< depends on x_j                  */

  R phi_prod[ths->d+1];

  int ix,ix_old;

  nfft_sort_nodes(ths);

  phi_prod[0]=1;
  ll_plain[0]=0;

  for(t=0,lprod = 1; t<ths->d; t++)
      lprod *= 2*ths->m+2;

  for(j=0,ix=0,ix_old=0; j<ths->M_total; j++)
    {
      MACRO_init_uo_l_lj_t;

      for(l_L=0; l_L<lprod; l_L++, ix++)
  {
    MACRO_update_phi_prod_ll_plain(without_PRE_PSI);

    ths->psi_index_g[ix]=ll_plain[ths->d];
    ths->psi[ix]=phi_prod[ths->d];

    MACRO_count_uo_l_lj_t;
  } /* for(l_L) */


      ths->psi_index_f[j]=ix-ix_old;
      ix_old=ix;
    } /* for(j) */
#endif
}

void nfft_precompute_one_psi(nfft_plan *ths)
{
  if(ths->nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(ths);
  if(ths->nfft_flags & PRE_FG_PSI)
    nfft_precompute_fg_psi(ths);
  if(ths->nfft_flags & PRE_PSI)
    nfft_precompute_psi(ths);
  if(ths->nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(ths);
}


static void nfft_init_help(nfft_plan *ths)
{
  int t;                                /**< index over all dimensions       */
  int lprod;                            /**< 'bandwidth' of matrix B         */

  if (ths->nfft_flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    ths->nfft_flags |= NFFT_SORT_NODES;

  ths->N_total=nfft_prod_int(ths->N, ths->d);
  ths->n_total=nfft_prod_int(ths->n, ths->d);

  ths->sigma = (R*) nfft_malloc(ths->d*sizeof(R));
  for(t = 0;t < ths->d; t++)
    ths->sigma[t] = ((R)ths->n[t])/ths->N[t];

  WINDOW_HELP_INIT;

  if(ths->nfft_flags & MALLOC_X)
    ths->x = (R*)nfft_malloc(ths->d*ths->M_total*sizeof(R));

  if(ths->nfft_flags & MALLOC_F_HAT)
    ths->f_hat = (fftw_complex*)nfft_malloc(ths->N_total*sizeof(C));

  if(ths->nfft_flags & MALLOC_F)
    ths->f = (fftw_complex*)nfft_malloc(ths->M_total*sizeof(C));

  if(ths->nfft_flags & PRE_PHI_HUT)
    nfft_precompute_phi_hut(ths);

  if(ths->nfft_flags & PRE_LIN_PSI)
  {
      ths->K=(1U<< 10)*(ths->m+2);
      ths->psi = (R*) nfft_malloc((ths->K+1)*ths->d*sizeof(R));
  }

  if(ths->nfft_flags & PRE_FG_PSI)
    ths->psi = (R*) nfft_malloc(ths->M_total*ths->d*2*sizeof(R));

  if(ths->nfft_flags & PRE_PSI)
    ths->psi = (R*) nfft_malloc(ths->M_total*ths->d*
             (2*ths->m+2)*sizeof(R));

  if(ths->nfft_flags & PRE_FULL_PSI)
  {
      for(t=0,lprod = 1; t<ths->d; t++)
    lprod *= 2*ths->m+2;

      ths->psi = (R*) nfft_malloc(ths->M_total*lprod*sizeof(R));

      ths->psi_index_f = (int*) nfft_malloc(ths->M_total*sizeof(int));
      ths->psi_index_g = (int*) nfft_malloc(ths->M_total*lprod*sizeof(int));
  }

  if(ths->nfft_flags & FFTW_INIT)
  {
#ifdef _OPENMP
    int nthreads = nfft_get_omp_num_threads();
#endif

    ths->g1=(fftw_complex*)nfft_malloc(ths->n_total*sizeof(C));

    if(ths->nfft_flags & FFT_OUT_OF_PLACE)
      ths->g2 = (fftw_complex*) nfft_malloc(ths->n_total*sizeof(C));
    else
      ths->g2 = ths->g1;

#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    fftw_plan_with_nthreads(nthreads);
#endif
    ths->my_fftw_plan1 = fftw_plan_dft(ths->d, ths->n, ths->g1, ths->g2, FFTW_FORWARD, ths->fftw_flags);
    ths->my_fftw_plan2 = fftw_plan_dft(ths->d, ths->n, ths->g2, ths->g1,
      FFTW_BACKWARD, ths->fftw_flags);
#ifdef _OPENMP
}
#endif
  }

  if(ths->nfft_flags & NFFT_SORT_NODES)
    ths->index_x = (int*) nfft_malloc(sizeof(int)*2*ths->M_total);
  else
    ths->index_x = NULL;

  ths->mv_trafo = (void (*) (void* ))nfft_trafo;
  ths->mv_adjoint = (void (*) (void* ))nfft_adjoint;
}

void nfft_init(nfft_plan *ths, int d, int *N, int M_total)
{
  int t;                                /**< index over all dimensions       */

  ths->d = d;

  ths->N=(int*) nfft_malloc(d*sizeof(int));

  for (t = 0;t < d; t++)
    ths->N[t] = N[t];

  ths->M_total = M_total;

  ths->n = (int*) nfft_malloc(d*sizeof(int));
  for (t = 0;t < d; t++)
    ths->n[t] = 2*X(next_power_of_2)(ths->N[t]);

  ths->m = WINDOW_HELP_ESTIMATE_m;

  if (d > 1)
  {
#ifdef _OPENMP
    ths->nfft_flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE | NFFT_SORT_NODES |
		      NFFT_OMP_BLOCKWISE_ADJOINT;
#else
    ths->nfft_flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE | NFFT_SORT_NODES;
#endif
  }
  else
    ths->nfft_flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE;

  ths->fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;

  nfft_init_help(ths);
}

void nfft_init_guru(nfft_plan *ths, int d, int *N, int M_total, int *n,
      int m, unsigned nfft_flags, unsigned fftw_flags)
{
  int t;                                /**< index over all dimensions       */

  ths->d =d;
  ths->N= (int*) nfft_malloc(ths->d*sizeof(int));
  for(t=0; t<d; t++)
    ths->N[t]= N[t];
  ths->M_total= M_total;
  ths->n= (int*) nfft_malloc(ths->d*sizeof(int));
  for(t=0; t<d; t++)
    ths->n[t]= n[t];
  ths->m= m;
  ths->nfft_flags= nfft_flags;
  ths->fftw_flags= fftw_flags;

  nfft_init_help(ths);
}

void nfft_init_1d(nfft_plan *ths, int N1, int M_total)
{
  int N[1];

  N[0]=N1;

  nfft_init(ths, 1, N, M_total);
}

void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M_total)
{
  int N[2];

  N[0]=N1;
  N[1]=N2;
  nfft_init(ths,2,N,M_total);
}

void nfft_init_3d(nfft_plan *ths, int N1, int N2, int N3, int M_total)
{
  int N[3];

  N[0]=N1;
  N[1]=N2;
  N[2]=N3;
  nfft_init(ths,3,N,M_total);
}

const char* nfft_check(nfft_plan *ths)
{
  int j;

  for(j=0;j<ths->M_total*ths->d;j++)
    if((ths->x[j]<-K(0.5)) || (ths->x[j]>= K(0.5)))
      return "ths->x out of range [-0.5,0.5)";

  for(j=0;j<ths->d;j++)
  {
    if(ths->sigma[j]<=1)
      return "nfft_check: oversampling factor too small";

    if(ths->N[j]<=ths->m)
      return "Polynomial degree N is smaller than cut-off m";

    if(ths->N[j]%2==1)
      return "polynomial degree N has to be even";
  }
  return 0;
}

void nfft_finalize(nfft_plan *ths)
{
  int t; /* index over dimensions */

  if(ths->nfft_flags & NFFT_SORT_NODES)
    nfft_free(ths->index_x);

  if(ths->nfft_flags & FFTW_INIT)
  {
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    fftw_destroy_plan(ths->my_fftw_plan2);
    fftw_destroy_plan(ths->my_fftw_plan1);
}

    if(ths->nfft_flags & FFT_OUT_OF_PLACE)
      nfft_free(ths->g2);

    nfft_free(ths->g1);
  }

  if(ths->nfft_flags & PRE_FULL_PSI)
    {
      nfft_free(ths->psi_index_g);
      nfft_free(ths->psi_index_f);
      nfft_free(ths->psi);
    }

  if(ths->nfft_flags & PRE_PSI)
    nfft_free(ths->psi);

  if(ths->nfft_flags & PRE_FG_PSI)
    nfft_free(ths->psi);

  if(ths->nfft_flags & PRE_LIN_PSI)
    nfft_free(ths->psi);

  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      for(t=0; t<ths->d; t++)
        nfft_free(ths->c_phi_inv[t]);
      nfft_free(ths->c_phi_inv);
    }

  if(ths->nfft_flags & MALLOC_F)
    nfft_free(ths->f);

  if(ths->nfft_flags & MALLOC_F_HAT)
    nfft_free(ths->f_hat);

  if(ths->nfft_flags & MALLOC_X)
  nfft_free(ths->x);

  WINDOW_HELP_FINALIZE;

  nfft_free(ths->sigma);
  nfft_free(ths->n);
  nfft_free(ths->N);
}
