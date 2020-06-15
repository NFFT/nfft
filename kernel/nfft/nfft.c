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

/* Nonequispaced FFT */

/* Authors: D. Potts, S. Kunis 2002-2009, Jens Keiner 2009, Toni Volkmer 2012 */

/* configure header */
#include "config.h"

/* complex datatype (maybe) */
#ifdef HAVE_COMPLEX_H
#include<complex.h>
#endif

/* NFFT headers */
#include "nfft3.h"
#include "infft.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef OMP_ASSERT
#include <assert.h>
#endif

#undef X
#define X(name) NFFT(name)

/** Compute aggregated product of integer array. */
static inline INT intprod(const INT *vec, const INT a, const INT d)
{
  INT t, p;

  p = 1;
  for (t = 0; t < d; t++)
    p *= vec[t] - a;

  return p;
}

/* handy shortcuts */
#define BASE(x) CEXP(x)

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
static inline void sort0(const INT d, const INT *n, const INT m,
    const INT local_x_num, const R *local_x, INT *ar_x)
{
  INT u_j[d], i, j, help, rhigh;
  INT *ar_x_temp;
  INT nprod;

  for (i = 0; i < local_x_num; i++)
  {
    ar_x[2 * i] = 0;
    ar_x[2 *i + 1] = i;
    for (j = 0; j < d; j++)
    {
      help = (INT) LRINT(FLOOR((R)(n[j]) * local_x[d * i + j] - (R)(m)));
      u_j[j] = (help % n[j] + n[j]) % n[j];

      ar_x[2 * i] += u_j[j];
      if (j + 1 < d)
        ar_x[2 * i] *= n[j + 1];
    }
  }

  for (j = 0, nprod = 1; j < d; j++)
    nprod *= n[j];

  rhigh = (INT) LRINT(CEIL(LOG2((R)nprod))) - 1;

  ar_x_temp = (INT*) Y(malloc)(2 * (size_t)(local_x_num) * sizeof(INT));
  Y(sort_node_indices_radix_lsdf)(local_x_num, ar_x, ar_x_temp, rhigh);
#ifdef OMP_ASSERT
  for (i = 1; i < local_x_num; i++)
    assert(ar_x[2 * (i - 1)] <= ar_x[2 * i]);
#endif
  Y(free)(ar_x_temp);
}

/**
 * Sort nodes (index) to get better cache utilization during multiplication
 * with matrix B.
 * The resulting index set is written to ths->index_x[2*j+1], the nodes array
 * remains unchanged.
 *
 * \arg ths nfft_plan
 */
static inline void sort(const X(plan) *ths)
{
  if (ths->flags & NFFT_SORT_NODES)
    sort0(ths->d, ths->n, ths->m, ths->M_total, ths->x, ths->index_x);
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
void X(trafo_direct)(const X(plan) *ths)
{
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;

  memset(f, 0, (size_t)(ths->M_total) * sizeof(C));

  if (ths->d == 1)
  {
    /* specialize for univariate case, rationale: faster */
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j)
#endif
    for (j = 0; j < ths->M_total; j++)
    {
      INT k_L;
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        R omega = K2PI * ((R)(k_L - ths->N_total/2)) * ths->x[j];
        f[j] += f_hat[k_L] * BASE(-II * omega);
      }
    }
  }
  else
  {
    /* multivariate case */
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j)
#endif
    for (j = 0; j < ths->M_total; j++)
    {
      R x[ths->d], omega, Omega[ths->d + 1];
      INT t, t2, k_L, k[ths->d];
      Omega[0] = K(0.0);
      for (t = 0; t < ths->d; t++)
      {
        k[t] = -ths->N[t]/2;
        x[t] = K2PI * ths->x[j * ths->d + t];
        Omega[t+1] = ((R)k[t]) * x[t] + Omega[t];
      }
      omega = Omega[ths->d];

      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        f[j] += f_hat[k_L] * BASE(-II * omega);
        {
          for (t = ths->d - 1; (t >= 1) && (k[t] == ths->N[t]/2 - 1); t--)
            k[t]-= ths->N[t]-1;

          k[t]++;

          for (t2 = t; t2 < ths->d; t2++)
            Omega[t2+1] = ((R)k[t2]) * x[t2] + Omega[t2];

          omega = Omega[ths->d];
        }
      }
    }
  }
}

void X(adjoint_direct)(const X(plan) *ths)
{
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;

  memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(C));

  if (ths->d == 1)
  {
    /* specialize for univariate case, rationale: faster */
#ifdef _OPENMP
      INT k_L;
      #pragma omp parallel for default(shared) private(k_L)
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        INT j;
        for (j = 0; j < ths->M_total; j++)
        {
          R omega = K2PI * ((R)(k_L - (ths->N_total/2))) * ths->x[j];
          f_hat[k_L] += f[j] * BASE(II * omega);
        }
      }
#else
      INT j;
      for (j = 0; j < ths->M_total; j++)
      {
        INT k_L;
        for (k_L = 0; k_L < ths->N_total; k_L++)
        {
          R omega = K2PI * ((R)(k_L - ths->N_total / 2)) * ths->x[j];
          f_hat[k_L] += f[j] * BASE(II * omega);
        }
      }
#endif
  }
  else
  {
    /* multivariate case */
    INT j, k_L;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j, k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      INT k[ths->d], k_temp, t;

      k_temp = k_L;

      for (t = ths->d - 1; t >= 0; t--)
      {
        k[t] = k_temp % ths->N[t] - ths->N[t]/2;
        k_temp /= ths->N[t];
      }

      for (j = 0; j < ths->M_total; j++)
      {
        R omega = K(0.0);
        for (t = 0; t < ths->d; t++)
          omega += k[t] * K2PI * ths->x[j * ths->d + t];
        f_hat[k_L] += f[j] * BASE(II * omega);
      }
    }
#else
    for (j = 0; j < ths->M_total; j++)
    {
      R x[ths->d], omega, Omega[ths->d+1];
      INT t, t2, k[ths->d];
      Omega[0] = K(0.0);
      for (t = 0; t < ths->d; t++)
      {
        k[t] = -ths->N[t]/2;
        x[t] = K2PI * ths->x[j * ths->d + t];
        Omega[t+1] = ((R)k[t]) * x[t] + Omega[t];
      }
      omega = Omega[ths->d];
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        f_hat[k_L] += f[j] * BASE(II * omega);

        for (t = ths->d-1; (t >= 1) && (k[t] == ths->N[t]/2-1); t--)
          k[t]-= ths->N[t]-1;

        k[t]++;

        for (t2 = t; t2 < ths->d; t2++)
          Omega[t2+1] = ((R)k[t2]) * x[t2] + Omega[t2];

        omega = Omega[ths->d];
      }
    }
#endif
  }
}

/** fast computation of non-equispaced fourier transforms
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
static inline void uo(const X(plan) *ths, const INT j, INT *up, INT *op,
  const INT act_dim)
{
  const R xj = ths->x[j * ths->d + act_dim];
  INT c = LRINT(FLOOR(xj * (R)(ths->n[act_dim])));

  (*up) = c - (ths->m);
  (*op) = c + 1 + (ths->m);
}

static inline void uo2(INT *u, INT *o, const R x, const INT n, const INT m)
{
  INT c = LRINT(FLOOR(x * (R)(n)));

  *u = (c - m + n) % n;
  *o = (c + 1 + m + n) % n;
}

#define MACRO_D_compute_A \
{ \
  g_hat[k_plain[ths->d]] = f_hat[ks_plain[ths->d]] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_compute_T \
{ \
  f_hat[ks_plain[ths->d]] = g_hat[k_plain[ths->d]] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_init_result_A memset(g_hat, 0, (size_t)(ths->n_total) * sizeof(C));

#define MACRO_D_init_result_T memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(C));

#define MACRO_with_PRE_PHI_HUT * ths->c_phi_inv[t2][ks[t2]];

#define MACRO_without_PRE_PHI_HUT / (PHI_HUT(ths->n[t2],ks[t2]-(ths->N[t2]/2),t2));

#define MACRO_init_k_ks \
{ \
  for (t = ths->d-1; 0 <= t; t--) \
  { \
    kp[t] = k[t] = 0; \
    ks[t] = ths->N[t]/2; \
  } \
  t++; \
}

#define MACRO_update_c_phi_inv_k(which_one) \
{ \
  for (t2 = t; t2 < ths->d; t2++) \
  { \
    c_phi_inv_k[t2+1] = c_phi_inv_k[t2] MACRO_ ##which_one; \
    ks_plain[t2+1] = ks_plain[t2]*ths->N[t2] + ks[t2]; \
    k_plain[t2+1] = k_plain[t2]*ths->n[t2] + k[t2]; \
  } \
}

#define MACRO_count_k_ks \
{ \
  for (t = ths->d-1; (t > 0) && (kp[t] == ths->N[t]-1); t--) \
  { \
    kp[t] = k[t] = 0; \
    ks[t]= ths->N[t]/2; \
  } \
\
  kp[t]++; k[t]++; ks[t]++; \
  if(kp[t] == ths->N[t]/2) \
  { \
    k[t] = ths->n[t] - ths->N[t]/2; \
    ks[t] = 0; \
  } \
} \

/* sub routines for the fast transforms  matrix vector multiplication with D, D^T */
#define MACRO_D(which_one) \
static inline void D_serial_ ## which_one (X(plan) *ths) \
{ \
  C *f_hat, *g_hat; /* local copy */ \
  R c_phi_inv_k[ths->d+1]; /* postfix product of PHI_HUT */ \
  INT t, t2; /* index dimensions */ \
  INT k_L; /* plain index */ \
  INT kp[ths->d]; /* multi index (simple) */ \
  INT k[ths->d]; /* multi index in g_hat */ \
  INT ks[ths->d]; /* multi index in f_hat, c_phi_inv*/ \
  INT k_plain[ths->d+1]; /* postfix plain index */ \
  INT ks_plain[ths->d+1]; /* postfix plain index */ \
 \
  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat; \
  MACRO_D_init_result_ ## which_one; \
\
  c_phi_inv_k[0] = K(1.0); \
  k_plain[0] = 0; \
  ks_plain[0] = 0; \
\
  MACRO_init_k_ks; \
\
  if (ths->flags & PRE_PHI_HUT) \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(with_PRE_PHI_HUT); \
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
  else \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(without_PRE_PHI_HUT); \
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
}

#ifdef _OPENMP
static inline void D_openmp_A(X(plan) *ths)
{
  C *f_hat, *g_hat;                     /**< local copy                     */
  INT k_L;                              /**< plain index                    */

  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat;
  memset(g_hat, 0, ths->n_total * sizeof(C));

  if (ths->flags & PRE_PHI_HUT)
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

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
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

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
        c_phi_inv_k_val /= (PHI_HUT(ths->n[t],ks[t]-(ths->N[t]/2),t));
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*ths->n[t] + k[t];
      }

      g_hat[k_plain_val] = f_hat[ks_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* else(PRE_PHI_HUT) */
}
#endif

#ifndef _OPENMP
MACRO_D(A)
#endif

static inline void D_A(X(plan) *ths)
{
#ifdef _OPENMP
  D_openmp_A(ths);
#else
  D_serial_A(ths);
#endif
}

#ifdef _OPENMP
static void D_openmp_T(X(plan) *ths)
{
  C *f_hat, *g_hat;                     /**< local copy                     */
  INT k_L;                              /**< plain index                    */

  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat;
  memset(f_hat, 0, ths->N_total * sizeof(C));

  if (ths->flags & PRE_PHI_HUT)
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

      for (t = ths->d - 1; t >= 0; t--)
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
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

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
        c_phi_inv_k_val /= (PHI_HUT(ths->n[t],ks[t]-(ths->N[t]/2),t));
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*ths->n[t] + k[t];
      }

      f_hat[ks_plain_val] = g_hat[k_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* else(PRE_PHI_HUT) */
}
#endif

#ifndef _OPENMP
MACRO_D(T)
#endif

static void D_T(X(plan) *ths)
{
#ifdef _OPENMP
  D_openmp_T(ths);
#else
  D_serial_T(ths);
#endif
}

/* sub routines for the fast transforms matrix vector multiplication with B, B^T */
#define MACRO_B_init_result_A memset(ths->f, 0, (size_t)(ths->M_total) * sizeof(C));
#define MACRO_B_init_result_T memset(ths->g, 0, (size_t)(ths->n_total) * sizeof(C));

#define MACRO_B_PRE_FULL_PSI_compute_A \
{ \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]]; \
}

#define MACRO_B_PRE_FULL_PSI_compute_T \
{ \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj); \
}

#define MACRO_B_compute_A \
{ \
  ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]]; \
}

#define MACRO_B_compute_T \
{ \
  ths->g[ll_plain[ths->d]] += phi_prod[ths->d] * ths->f[j]; \
}

#define MACRO_with_FG_PSI fg_psi[t2][lj[t2]]

#define MACRO_with_PRE_PSI ths->psi[(j*ths->d+t2) * (2*ths->m+2)+lj[t2]]

#define MACRO_without_PRE_PSI_improved psij_const[t2 * (2*ths->m+2) + lj[t2]]

#define MACRO_without_PRE_PSI  PHI(ths->n[t2], ths->x[j*ths->d+t2] \
  - ((R) (lj[t2]+u[t2]))/((R)ths->n[t2]), t2)

#define MACRO_init_uo_l_lj_t \
INT l_all[ths->d*(2*ths->m+2)]; \
{ \
  for (t = ths->d-1; t >= 0; t--) \
  { \
    uo(ths,j,&u[t],&o[t],t); \
    INT lj_t; \
    for (lj_t = 0; lj_t < 2*ths->m+2; lj_t++) \
      l_all[t*(2*ths->m+2) + lj_t] = (u[t] + lj_t + ths->n[t]) % ths->n[t]; \
    lj[t] = 0; \
  } \
  t++; \
}

#define MACRO_update_phi_prod_ll_plain(which_one) { \
  for (t2 = t; t2 < ths->d; t2++) \
    { \
      phi_prod[t2+1] = phi_prod[t2] * MACRO_ ## which_one; \
      ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
    } \
}

#define MACRO_count_uo_l_lj_t \
{ \
  for (t = ths->d-1; (t > 0) && (lj[t] == o[t]-u[t]); t--) \
  { \
    lj[t] = 0; \
  } \
 \
  lj[t]++; \
}

#define MACRO_COMPUTE_with_PRE_PSI MACRO_with_PRE_PSI
#define MACRO_COMPUTE_with_PRE_FG_PSI MACRO_with_FG_PSI
#define MACRO_COMPUTE_with_FG_PSI MACRO_with_FG_PSI
#define MACRO_COMPUTE_with_PRE_LIN_PSI MACRO_with_FG_PSI
#define MACRO_COMPUTE_without_PRE_PSI MACRO_without_PRE_PSI_improved
#define MACRO_COMPUTE_without_PRE_PSI_improved MACRO_without_PRE_PSI_improved

#define MACRO_B_COMPUTE_ONE_NODE(whichone_AT,whichone_FLAGS) \
      if (ths->d == 4) \
      { \
        INT l0, l1, l2, l3; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                MACRO_B_compute_ ## whichone_AT; \
              } \
            } \
          } \
        } \
      } /* if(d==4) */ \
      else if (ths->d == 5) \
      { \
        INT l0, l1, l2, l3, l4; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
                for (l4 = 0; l4 < 2*ths->m+2; l4++) \
                { \
                  lj[4] = l4; \
                  t2 = 4; \
                  phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone_FLAGS; \
                  ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                  MACRO_B_compute_ ## whichone_AT; \
                } \
              } \
            } \
          } \
        } \
      } /* if(d==5) */ \
      else { \
        for (l_L = 0; l_L < lprod; l_L++) \
        { \
          MACRO_update_phi_prod_ll_plain(whichone_FLAGS); \
 \
          MACRO_B_compute_ ## whichone_AT; \
 \
          MACRO_count_uo_l_lj_t; \
        } /* for(l_L) */ \
      }

#define MACRO_B(which_one) \
static inline void B_serial_ ## which_one (X(plan) *ths) \
{ \
  INT lprod; /* 'regular bandwidth' of matrix B  */ \
  INT u[ths->d], o[ths->d]; /* multi band with respect to x_j */ \
  INT t, t2; /* index dimensions */ \
  INT k; /* index nodes */ \
  INT l_L, ix; /* index one row of B */ \
  INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */ \
  INT ll_plain[ths->d+1]; /* postfix plain index in g */ \
  R phi_prod[ths->d+1]; /* postfix product of PHI */ \
  R y[ths->d]; \
  R fg_psi[ths->d][2*ths->m+2]; \
  R fg_exp_l[ths->d][2*ths->m+2]; \
  INT l_fg,lj_fg; \
  R tmpEXP1, tmpEXP2, tmpEXP2sq, tmp1, tmp2, tmp3; \
  R ip_w; \
  INT ip_u; \
  INT ip_s = ths->K/(ths->m+2); \
 \
  MACRO_B_init_result_ ## which_one; \
 \
  if (ths->flags & PRE_FULL_PSI) \
  { \
    INT j; \
    C *f, *g; /* local copy */ \
    C *fj; /* local copy */ \
    f = (C*)ths->f; g = (C*)ths->g; \
 \
    for (ix = 0, j = 0, fj = f; j < ths->M_total; j++, fj++) \
    { \
      for (l_L = 0; l_L < ths->psi_index_f[j]; l_L++, ix++) \
      { \
        MACRO_B_PRE_FULL_PSI_compute_ ## which_one; \
      } \
    } \
    return; \
  } \
\
  phi_prod[0] = K(1.0); \
  ll_plain[0] = 0; \
\
  for (t = 0, lprod = 1; t < ths->d; t++) \
    lprod *= (2 * ths->m + 2); \
\
  if (ths->flags & PRE_PSI) \
  { \
    sort(ths); \
 \
    for (k = 0; k < ths->M_total; k++) \
    { \
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      MACRO_B_COMPUTE_ONE_NODE(which_one,with_PRE_PSI); \
    } /* for(j) */ \
    return; \
  } /* if(PRE_PSI) */ \
 \
  if (ths->flags & PRE_FG_PSI) \
  { \
    sort(ths); \
 \
    for(t2 = 0; t2 < ths->d; t2++) \
    { \
      tmpEXP2 = EXP(K(-1.0) / ths->b[t2]); \
      tmpEXP2sq = tmpEXP2*tmpEXP2; \
      tmp2 = K(1.0); \
      tmp3 = K(1.0); \
      fg_exp_l[t2][0] = K(1.0); \
      for (lj_fg = 1; lj_fg <= (2 * ths->m + 2); lj_fg++) \
      { \
        tmp3 = tmp2*tmpEXP2; \
        tmp2 *= tmpEXP2sq; \
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1] * tmp3; \
      } \
    } \
    for (k = 0; k < ths->M_total; k++) \
    { \
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        fg_psi[t2][0] = ths->psi[2*(j*ths->d+t2)]; \
        tmpEXP1 = ths->psi[2*(j*ths->d+t2)+1]; \
        tmp1 = K(1.0); \
        for (l_fg = u[t2]+1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
        } \
      } \
 \
      MACRO_B_COMPUTE_ONE_NODE(which_one,with_FG_PSI); \
    } /* for(j) */ \
    return; \
  } /* if(PRE_FG_PSI) */ \
 \
  if (ths->flags & FG_PSI) \
  { \
    sort(ths); \
 \
    for (t2 = 0; t2 < ths->d; t2++) \
    { \
      tmpEXP2 = EXP(K(-1.0)/ths->b[t2]); \
      tmpEXP2sq = tmpEXP2*tmpEXP2; \
      tmp2 = K(1.0); \
      tmp3 = K(1.0); \
      fg_exp_l[t2][0] = K(1.0); \
      for (lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++) \
      { \
        tmp3 = tmp2*tmpEXP2; \
        tmp2 *= tmpEXP2sq; \
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3; \
      } \
    } \
    for (k = 0; k < ths->M_total; k++) \
    { \
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        fg_psi[t2][0] = (PHI(ths->n[t2], (ths->x[j*ths->d+t2] - ((R)u[t2])/((R)(ths->n[t2]))), t2));\
 \
        tmpEXP1 = EXP(K(2.0) * ((R)(ths->n[t2]) * ths->x[j * ths->d + t2] - (R)(u[t2])) \
          /ths->b[t2]); \
        tmp1 = K(1.0); \
        for (l_fg = u[t2] + 1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
        } \
      } \
 \
      MACRO_B_COMPUTE_ONE_NODE(which_one,with_FG_PSI); \
    } /* for(j) */ \
    return; \
  } /* if(FG_PSI) */ \
 \
  if (ths->flags & PRE_LIN_PSI) \
  { \
    sort(ths); \
 \
    for (k = 0; k<ths->M_total; k++) \
    { \
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        y[t2] = (((R)(ths->n[t2]) * ths->x[j * ths->d + t2] - (R)(u[t2])) \
          * ((R)(ths->K))) / (R)(ths->m + 2); \
        ip_u  = LRINT(FLOOR(y[t2])); \
        ip_w  = y[t2]-ip_u; \
        for (l_fg = u[t2], lj_fg = 0; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s)] \
            * (1-ip_w) + ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s+1)] \
            * (ip_w); \
        } \
      } \
 \
      MACRO_B_COMPUTE_ONE_NODE(which_one,with_FG_PSI); \
    } /* for(j) */ \
    return; \
  } /* if(PRE_LIN_PSI) */ \
 \
  sort(ths); \
 \
  /* no precomputed psi at all */ \
  for (k = 0; k < ths->M_total; k++) \
  { \
    INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k; \
 \
    R psij_const[ths->d * (2*ths->m+2)]; \
 \
    MACRO_init_uo_l_lj_t; \
 \
    for (t2 = 0; t2 < ths->d; t2++) \
    { \
      INT lj_t; \
      for (lj_t = 0; lj_t < 2*ths->m+2; lj_t++) \
        psij_const[t2 * (2*ths->m+2) + lj_t] = PHI(ths->n[t2], ths->x[j*ths->d+t2] \
                - ((R) (lj_t+u[t2]))/((R)ths->n[t2]), t2); \
    } \
 \
    MACRO_B_COMPUTE_ONE_NODE(which_one,without_PRE_PSI_improved); \
  } /* for(j) */ \
} /* nfft_B */ \

#ifndef _OPENMP
MACRO_B(A)
#endif

#ifdef _OPENMP
#define MACRO_B_openmp_A_COMPUTE_BEFORE_LOOP_with_PRE_PSI
#define MACRO_B_openmp_A_COMPUTE_UPDATE_with_PRE_PSI \
  MACRO_update_phi_prod_ll_plain(with_PRE_PSI);

#define MACRO_B_openmp_A_COMPUTE_INIT_FG_PSI \
    for (t2 = 0; t2 < ths->d; t2++) \
    { \
      INT lj_fg; \
      R tmpEXP2 = EXP(K(-1.0)/ths->b[t2]); \
      R tmpEXP2sq = tmpEXP2*tmpEXP2; \
      R tmp2 = K(1.0); \
      R tmp3 = K(1.0); \
      fg_exp_l[t2][0] = K(1.0); \
      for(lj_fg = 1; lj_fg <= (2*ths->m+2); lj_fg++) \
      { \
        tmp3 = tmp2*tmpEXP2; \
        tmp2 *= tmpEXP2sq; \
        fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3; \
      } \
    }
#define MACRO_B_openmp_A_COMPUTE_BEFORE_LOOP_with_PRE_FG_PSI \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        fg_psi[t2][0] = ths->psi[2*(j*ths->d+t2)]; \
        tmpEXP1 = ths->psi[2*(j*ths->d+t2)+1]; \
        tmp1 = K(1.0); \
        for (l_fg = u[t2]+1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
        } \
      }
#define MACRO_B_openmp_A_COMPUTE_UPDATE_with_PRE_FG_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_B_openmp_A_COMPUTE_BEFORE_LOOP_with_FG_PSI \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        fg_psi[t2][0] = (PHI(ths->n[t2],(ths->x[j*ths->d+t2]-((R)u[t2])/((R)ths->n[t2])),t2)); \
 \
        tmpEXP1 = EXP(K(2.0)*(ths->n[t2]*ths->x[j*ths->d+t2] - u[t2]) \
          /ths->b[t2]); \
        tmp1 = K(1.0); \
        for (l_fg = u[t2] + 1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
        } \
      }
#define MACRO_B_openmp_A_COMPUTE_UPDATE_with_FG_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_B_openmp_A_COMPUTE_BEFORE_LOOP_with_PRE_LIN_PSI \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        y[t2] = ((ths->n[t2]*ths->x[j*ths->d+t2]-(R)u[t2]) \
          * ((R)ths->K))/(ths->m+2); \
        ip_u  = LRINT(FLOOR(y[t2])); \
        ip_w  = y[t2]-ip_u; \
        for (l_fg = u[t2], lj_fg = 0; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s)] \
            * (1-ip_w) + ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s+1)] \
            * (ip_w); \
        } \
      }
#define MACRO_B_openmp_A_COMPUTE_UPDATE_with_PRE_LIN_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_B_openmp_A_COMPUTE_BEFORE_LOOP_without_PRE_PSI \
    for (t2 = 0; t2 < ths->d; t2++) \
    { \
      INT lj_t; \
      for (lj_t = 0; lj_t < 2*ths->m+2; lj_t++) \
        psij_const[t2 * (2*ths->m+2) + lj_t] = PHI(ths->n[t2], ths->x[j*ths->d+t2] \
                - ((R) (lj_t+u[t2]))/((R)ths->n[t2]), t2); \
    }
#define MACRO_B_openmp_A_COMPUTE_UPDATE_without_PRE_PSI \
  MACRO_update_phi_prod_ll_plain(without_PRE_PSI_improved);

#define MACRO_B_openmp_A_COMPUTE(whichone) \
{ \
      INT u[ths->d], o[ths->d]; /* multi band with respect to x_j */ \
      INT l_L; /* index one row of B */ \
      INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */ \
      INT ll_plain[ths->d+1]; /* postfix plain index in g */ \
      R phi_prod[ths->d+1]; /* postfix product of PHI */ \
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k; \
 \
      phi_prod[0] = K(1.0); \
      ll_plain[0] = 0; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      MACRO_B_openmp_A_COMPUTE_BEFORE_LOOP_ ##whichone \
 \
      if (ths->d == 4) \
      { \
        INT l0, l1, l2, l3; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]]; \
              } \
            } \
          } \
        } \
      } /* if(d==4) */ \
      else if (ths->d == 5) \
      { \
        INT l0, l1, l2, l3, l4; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
                for (l4 = 0; l4 < 2*ths->m+2; l4++) \
                { \
                  lj[4] = l4; \
                  t2 = 4; \
                  phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                  ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                  ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]]; \
                } \
              } \
            } \
          } \
        } \
      } /* if(d==5) */ \
      else { \
        for (l_L = 0; l_L < lprod; l_L++) \
        { \
          MACRO_B_openmp_A_COMPUTE_UPDATE_ ##whichone \
 \
          ths->f[j] += phi_prod[ths->d] * ths->g[ll_plain[ths->d]]; \
 \
          MACRO_count_uo_l_lj_t; \
        } /* for(l_L) */ \
      } \
}

static inline void B_openmp_A (X(plan) *ths)
{
  INT lprod; /* 'regular bandwidth' of matrix B  */
  INT k;

  memset(ths->f, 0, ths->M_total * sizeof(C));

  for (k = 0, lprod = 1; k < ths->d; k++)
    lprod *= (2*ths->m+2);

  if (ths->flags & PRE_FULL_PSI)
  {
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < ths->M_total; k++)
    {
      INT l;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < lprod; l++)
        ths->f[j] += ths->psi[j*lprod+l] * ths->g[ths->psi_index_g[j*lprod+l]];
    }
    return;
  }

  if (ths->flags & PRE_PSI)
  {
    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < ths->M_total; k++)
    {
      INT t, t2; /* index dimensions */
      MACRO_B_openmp_A_COMPUTE(with_PRE_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_PSI) */

  if (ths->flags & PRE_FG_PSI)
  {
    INT t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];

    MACRO_B_openmp_A_COMPUTE_INIT_FG_PSI

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      R fg_psi[ths->d][2*ths->m+2];
      R tmpEXP1, tmp1;
      INT l_fg,lj_fg;

      MACRO_B_openmp_A_COMPUTE(with_PRE_FG_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_FG_PSI) */

  if (ths->flags & FG_PSI)
  {
    INT t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];

    sort(ths);

    MACRO_B_openmp_A_COMPUTE_INIT_FG_PSI

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      R fg_psi[ths->d][2*ths->m+2];
      R tmpEXP1, tmp1;
      INT l_fg,lj_fg;

      MACRO_B_openmp_A_COMPUTE(with_FG_PSI);
    } /* for(j) */
    return;
  } /* if(FG_PSI) */

  if (ths->flags & PRE_LIN_PSI)
  {
    sort(ths);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k<ths->M_total; k++)
    {
      INT t, t2; /* index dimensions */
      R y[ths->d];
      R fg_psi[ths->d][2*ths->m+2];
      INT l_fg,lj_fg;
      R ip_w;
      INT ip_u;
      INT ip_s = ths->K/(ths->m+2);

      MACRO_B_openmp_A_COMPUTE(with_PRE_LIN_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < ths->M_total; k++)
  {
    INT t, t2; /* index dimensions */
    R psij_const[ths->d * (2*ths->m+2)];

    MACRO_B_openmp_A_COMPUTE(without_PRE_PSI);
  } /* for(j) */
}
#endif

static void B_A(X(plan) *ths)
{
#ifdef _OPENMP
  B_openmp_A(ths);
#else
  B_serial_A(ths);
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
static inline INT index_x_binary_search(const INT *ar_x, const INT len, const INT key)
{
  INT left = 0, right = len - 1;

  if (len == 1)
    return 0;

  while (left < right - 1)
  {
    INT i = (left + right) / 2;
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
static void nfft_adjoint_B_omp_blockwise_init(INT *my_u0, INT *my_o0,
    INT *min_u_a, INT *max_u_a, INT *min_u_b, INT *max_u_b, const INT d,
    const INT *n, const INT m)
{
  const INT n0 = n[0];
  INT k;
  INT nthreads = omp_get_num_threads();
  INT nthreads_used = MIN(nthreads, n0);
  INT size_per_thread = n0 / nthreads_used;
  INT size_left = n0 - size_per_thread * nthreads_used;
  INT size_g[nthreads_used];
  INT offset_g[nthreads_used];
  INT my_id = omp_get_thread_num();
  INT n_prod_rest = 1;

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
    const INT m22 = 2 * m + 2;

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
#ifdef OMP_ASSERT
    assert(*min_u_a <= *max_u_a);
    assert(*min_u_b <= *max_u_b);
    assert(*min_u_b == -1 || *max_u_a < *min_u_b);
#endif
  }
}
#endif

/**
 * Calculates adjoint NFFT for flag PRE_FULL_PSI.
 * Parallel calculation (OpenMP) with and without atomic operations.
 *
 * \arg lprod stride (2*m+2)^d
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_B_compute_full_psi(C *g, const INT *psi_index_g,
    const R *psi, const C *f, const INT M, const INT d, const INT *n,
    const INT m, const unsigned flags, const INT *index_x)
{
  INT k;
  INT lprod;
#ifdef _OPENMP
  INT lprod_m1;
#endif
#ifndef _OPENMP
  UNUSED(n);
#endif
  {
    INT t;
    for(t = 0, lprod = 1; t < d; t++)
        lprod *= 2 * m + 2;
  }
#ifdef _OPENMP
  lprod_m1 = lprod / (2 * m + 2);
#endif

#ifdef _OPENMP
  if (flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel private(k)
    {
      INT my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b;
      const INT *ar_x = index_x;
      INT n_prod_rest = 1;

      for (k = 1; k < d; k++)
        n_prod_rest *= n[k];

      nfft_adjoint_B_omp_blockwise_init(&my_u0, &my_o0, &min_u_a, &max_u_a, &min_u_b, &max_u_b, d, n, m);

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
          INT l0, lrest;
          INT u_prod = ar_x[2*k];
          INT j = ar_x[2*k+1];

          if (u_prod < min_u_a || u_prod > max_u_a)
            break;

          for (l0 = 0; l0 < 2 * m + 2; l0++)
          {
            const INT start_index = psi_index_g[j * lprod + l0 * lprod_m1];

            if (start_index < my_u0 * n_prod_rest || start_index > (my_o0+1) * n_prod_rest - 1)
              continue;

            for (lrest = 0; lrest < lprod_m1; lrest++)
            {
              const INT l = l0 * lprod_m1 + lrest;
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
          INT l0, lrest;
          INT u_prod = ar_x[2*k];
          INT j = ar_x[2*k+1];

          if (u_prod < min_u_b || u_prod > max_u_b)
            break;

          for (l0 = 0; l0 < 2 * m + 2; l0++)
          {
            const INT start_index = psi_index_g[j * lprod + l0 * lprod_m1];

            if (start_index < my_u0 * n_prod_rest || start_index > (my_o0+1) * n_prod_rest - 1)
              continue;
            for (lrest = 0; lrest < lprod_m1; lrest++)
            {
              const INT l = l0 * lprod_m1 + lrest;
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

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT l;
    INT j = (flags & NFFT_SORT_NODES) ? index_x[2*k+1] : k;

    for (l = 0; l < lprod; l++)
    {
#ifdef _OPENMP
      C val = psi[j * lprod + l] * f[j];
      C *gref = g + psi_index_g[j * lprod + l];
      R *gref_real = (R*) gref;

      #pragma omp atomic
      gref_real[0] += CREAL(val);

      #pragma omp atomic
      gref_real[1] += CIMAG(val);
#else
      g[psi_index_g[j * lprod + l]] += psi[j * lprod + l] * f[j];
#endif
    }
  }
}

#ifndef _OPENMP
MACRO_B(T)
#endif


#ifdef _OPENMP

#ifdef OMP_ASSERT
#define MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_A \
{ \
          assert(ar_x[2*k] >= min_u_a || k == M-1); \
          if (k > 0) \
            assert(ar_x[2*k-2] < min_u_a); \
}
#else
#define MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_A
#endif

#ifdef OMP_ASSERT
#define MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_B \
{ \
          assert(ar_x[2*k] >= min_u_b || k == M-1); \
          if (k > 0) \
            assert(ar_x[2*k-2] < min_u_b); \
}
#else
#define MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_B
#endif

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_PRE_PSI
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_PRE_PSI \
  MACRO_update_phi_prod_ll_plain(with_PRE_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_PRE_FG_PSI \
      R fg_psi[ths->d][2*ths->m+2]; \
      R tmpEXP1, tmp1; \
      INT l_fg,lj_fg; \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        fg_psi[t2][0] = ths->psi[2*(j*ths->d+t2)]; \
        tmpEXP1 = ths->psi[2*(j*ths->d+t2)+1]; \
        tmp1 = K(1.0); \
        for (l_fg = u[t2]+1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
        } \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_PRE_FG_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_FG_PSI \
      R fg_psi[ths->d][2*ths->m+2]; \
      R tmpEXP1, tmp1; \
      INT l_fg,lj_fg; \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        fg_psi[t2][0] = (PHI(ths->n[t2],(ths->x[j*ths->d+t2]-((R)u[t2])/((R)ths->n[t2])),t2)); \
 \
        tmpEXP1 = EXP(K(2.0)*((R)ths->n[t2]*ths->x[j*ths->d+t2] - (R)u[t2]) \
          /ths->b[t2]); \
        tmp1 = K(1.0); \
        for (l_fg = u[t2] + 1, lj_fg = 1; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
        } \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_FG_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_PRE_LIN_PSI \
      R y[ths->d]; \
      R fg_psi[ths->d][2*ths->m+2]; \
      INT l_fg,lj_fg; \
      R ip_w; \
      INT ip_u; \
      INT ip_s = ths->K/(ths->m+2); \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        y[t2] = ((((R)ths->n[t2])*ths->x[j*ths->d+t2]-(R)u[t2]) \
          * ((R)ths->K))/((R)ths->m+2); \
        ip_u  = LRINT(FLOOR(y[t2])); \
        ip_w  = y[t2]-ip_u; \
        for (l_fg = u[t2], lj_fg = 0; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s)] \
            * (1-ip_w) + ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s+1)] \
            * (ip_w); \
        } \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_PRE_LIN_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_without_PRE_PSI \
      R psij_const[ths->d * (2*ths->m+2)]; \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        INT lj_t; \
        for (lj_t = 0; lj_t < 2*ths->m+2; lj_t++) \
          psij_const[t2 * (2*ths->m+2) + lj_t] = PHI(ths->n[t2], ths->x[j*ths->d+t2] \
                  - ((R) (lj_t+u[t2]))/((R)ths->n[t2]), t2); \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_without_PRE_PSI \
  MACRO_update_phi_prod_ll_plain(without_PRE_PSI_improved);

#define MACRO_adjoint_nd_B_OMP_BLOCKWISE_COMPUTE(whichone) \
{ \
      INT u[ths->d], o[ths->d]; /* multi band with respect to x_j */ \
      INT t, t2; /* index dimensions */ \
      INT l_L; /* index one row of B */ \
      INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */ \
      INT ll_plain[ths->d+1]; /* postfix plain index in g */ \
      R phi_prod[ths->d+1]; /* postfix product of PHI */ \
 \
      phi_prod[0] = K(1.0); \
      ll_plain[0] = 0; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_ ##whichone \
 \
      if (ths->d == 4) \
      { \
        INT l0, l1, l2, l3; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          if (l_all[lj[0]] < my_u0 || l_all[lj[0]] > my_o0) \
            continue; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                ths->g[ll_plain[ths->d]] += phi_prod[ths->d] * ths->f[j]; \
              } \
            } \
          } \
        } \
      } /* if(d==4) */ \
      else if (ths->d == 5) \
      { \
        INT l0, l1, l2, l3, l4; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          if (l_all[lj[0]] < my_u0 || l_all[lj[0]] > my_o0) \
            continue; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
                for (l4 = 0; l4 < 2*ths->m+2; l4++) \
                { \
                  lj[4] = l4; \
                  t2 = 4; \
                  phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                  ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                  ths->g[ll_plain[ths->d]] += phi_prod[ths->d] * ths->f[j]; \
                } \
              } \
            } \
          } \
        } \
      } /* if(d==5) */ \
      else { \
        l_L = 0; \
        while (l_L < lprod) \
        { \
          if (t == 0 && (l_all[lj[0]] < my_u0 || l_all[lj[0]] > my_o0)) \
          { \
            lj[0]++; \
            l_L += lprodrest; \
            continue; \
          } \
          MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_ ##whichone \
          ths->g[ll_plain[ths->d]] += phi_prod[ths->d] * ths->f[j]; \
          MACRO_count_uo_l_lj_t; \
          l_L++; \
        } /* for(l_L) */ \
      } \
}

#define MACRO_adjoint_nd_B_OMP_BLOCKWISE(whichone) \
{ \
    if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT) \
    { \
      INT lprodrest = 1; \
      for (k = 1; k < ths->d; k++) \
        lprodrest *= (2*ths->m+2); \
      _Pragma("omp parallel private(k)") \
      { \
        INT my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b; \
        INT *ar_x = ths->index_x; \
 \
        nfft_adjoint_B_omp_blockwise_init(&my_u0, &my_o0, &min_u_a, &max_u_a, \
            &min_u_b, &max_u_b, ths->d, ths->n, ths->m); \
 \
        if (min_u_a != -1) \
        { \
          k = index_x_binary_search(ar_x, ths->M_total, min_u_a); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_A \
 \
          while (k < ths->M_total) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_a || u_prod > max_u_a) \
              break; \
 \
            MACRO_adjoint_nd_B_OMP_BLOCKWISE_COMPUTE(whichone) \
 \
            k++; \
          } \
        } \
 \
        if (min_u_b != -1) \
        { \
          INT k = index_x_binary_search(ar_x, ths->M_total, min_u_b); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_B \
 \
          while (k < ths->M_total) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_b || u_prod > max_u_b) \
              break; \
 \
            MACRO_adjoint_nd_B_OMP_BLOCKWISE_COMPUTE(whichone) \
 \
            k++; \
          } \
        } \
      } /* omp parallel */ \
      return; \
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */ \
}

#define MACRO_adjoint_nd_B_OMP_COMPUTE(whichone) \
{ \
      INT u[ths->d], o[ths->d]; /* multi band with respect to x_j */ \
      INT l_L; /* index one row of B */ \
      INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */ \
      INT ll_plain[ths->d+1]; /* postfix plain index in g */ \
      R phi_prod[ths->d+1]; /* postfix product of PHI */ \
 \
      phi_prod[0] = K(1.0); \
      ll_plain[0] = 0; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_ ## whichone \
 \
      if (ths->d == 4) \
      { \
        INT l0, l1, l2, l3; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                C *lhs = ths->g + ll_plain[ths->d]; \
                R *lhs_real = (R*)lhs; \
                C val = phi_prod[ths->d] * ths->f[j]; \
 \
                _Pragma("omp atomic") \
                lhs_real[0] += CREAL(val); \
 \
                _Pragma("omp atomic") \
                lhs_real[1] += CIMAG(val); \
              } \
            } \
          } \
        } \
      } /* if(d==4) */ \
      else if (ths->d == 5) \
      { \
        INT l0, l1, l2, l3, l4; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
                for (l4 = 0; l4 < 2*ths->m+2; l4++) \
                { \
                  lj[4] = l4; \
                  t2 = 4; \
                  phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                  ll_plain[t2+1] = ll_plain[t2] * ths->n[t2] + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                  C *lhs = ths->g + ll_plain[ths->d]; \
                  R *lhs_real = (R*)lhs; \
                  C val = phi_prod[ths->d] * ths->f[j]; \
 \
                  _Pragma("omp atomic") \
                  lhs_real[0] += CREAL(val); \
 \
                  _Pragma("omp atomic") \
                  lhs_real[1] += CIMAG(val); \
                } \
              } \
            } \
          } \
        } \
      } /* if(d==5) */ \
      else { \
        for (l_L = 0; l_L < lprod; l_L++) \
        { \
          C *lhs; \
          R *lhs_real; \
          C val; \
 \
          MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_ ## whichone \
 \
          lhs = ths->g + ll_plain[ths->d]; \
          lhs_real = (R*)lhs; \
          val = phi_prod[ths->d] * ths->f[j]; \
 \
          _Pragma("omp atomic") \
          lhs_real[0] += CREAL(val); \
 \
          _Pragma("omp atomic") \
          lhs_real[1] += CIMAG(val); \
 \
          MACRO_count_uo_l_lj_t; \
        } /* for(l_L) */ \
      } \
}

static inline void B_openmp_T(X(plan) *ths)
{
  INT lprod; /* 'regular bandwidth' of matrix B  */
  INT k;

  memset(ths->g, 0, (size_t)(ths->n_total) * sizeof(C));

  for (k = 0, lprod = 1; k < ths->d; k++)
    lprod *= (2*ths->m+2);

  if (ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(ths->g, ths->psi_index_g, ths->psi, ths->f,
        ths->M_total, ths->d, ths->n, ths->m, ths->flags, ths->index_x);
    return;
  }

  if (ths->flags & PRE_PSI)
  {
    MACRO_adjoint_nd_B_OMP_BLOCKWISE(with_PRE_PSI);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < ths->M_total; k++)
    {
      INT t, t2; /* index dimensions */ \
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_PRE_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_PSI) */

  if (ths->flags & PRE_FG_PSI)
  {
    INT t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];
    for(t2 = 0; t2 < ths->d; t2++)
    {
      INT lj_fg;
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

    MACRO_adjoint_nd_B_OMP_BLOCKWISE(with_PRE_FG_PSI);

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_PRE_FG_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_FG_PSI) */

  if (ths->flags & FG_PSI)
  {
    INT t, t2; /* index dimensions */
    R fg_exp_l[ths->d][2*ths->m+2];

    sort(ths);

    for (t2 = 0; t2 < ths->d; t2++)
    {
      INT lj_fg;
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

    MACRO_adjoint_nd_B_OMP_BLOCKWISE(with_FG_PSI);

    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_FG_PSI);
    } /* for(j) */
    return;
  } /* if(FG_PSI) */

  if (ths->flags & PRE_LIN_PSI)
  {
    sort(ths);

    MACRO_adjoint_nd_B_OMP_BLOCKWISE(with_PRE_LIN_PSI);

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k<ths->M_total; k++)
    {
      INT t, t2; /* index dimensions */
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_PRE_LIN_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

  MACRO_adjoint_nd_B_OMP_BLOCKWISE(without_PRE_PSI);

  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < ths->M_total; k++)
  {
    INT t, t2; /* index dimensions */
    INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
    MACRO_adjoint_nd_B_OMP_COMPUTE(without_PRE_PSI);
  } /* for(j) */
}
#endif

static void B_T(X(plan) *ths)
{
#ifdef _OPENMP
  B_openmp_T(ths);
#else
  B_serial_T(ths);
#endif
}

/* ## specialized version for d=1  ########################################### */

static void nfft_1d_init_fg_exp_l(R *fg_exp_l, const INT m, const R b)
{
  const INT tmp2 = 2*m+2;
  INT l;
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
  const R *xj, const INT n, const INT m)
{
  INT u, o, l;
  const C *gj;
  const R *psij;
  psij = psij_const;

  uo2(&u, &o, *xj, n, m);

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

#ifndef _OPENMP
static void nfft_adjoint_1d_compute_serial(const C *fj, C *g,
    const R *psij_const, const R *xj, const INT n, const INT m)
{
  INT u,o,l;
  C *gj;
  const R *psij;
  psij = psij_const;

  uo2(&u,&o,*xj, n, m);

  if (u < o)
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
#endif

#ifdef _OPENMP
/* adjoint NFFT one-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_1d_compute_omp_atomic(const C f, C *g,
    const R *psij_const, const R *xj, const INT n, const INT m)
{
  INT u,o,l;
  C *gj;
  INT index_temp[2*m+2];

  uo2(&u,&o,*xj, n, m);

  for (l=0; l<=2*m+1; l++)
    index_temp[l] = (l+u)%n;

  for (l = 0, gj = g+u; l <= 2*m+1; l++)
  {
    INT i = index_temp[l];
    C *lhs = g+i;
    R *lhs_real = (R*)lhs;
    C val = psij_const[l] * f;
    #pragma omp atomic
    lhs_real[0] += CREAL(val);

    #pragma omp atomic
    lhs_real[1] += CIMAG(val);
  }
}
#endif

#ifdef _OPENMP
/**
 * Adjoint NFFT for one-dimensional case updating only a specified range of
 * vector g.
 *
 * \arg f input coefficient f[j]
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
static void nfft_adjoint_1d_compute_omp_blockwise(const C f, C *g,
    const R *psij_const, const R *xj, const INT n, const INT m,
    const INT my_u0, const INT my_o0)
{
  INT ar_u,ar_o,l;

  uo2(&ar_u,&ar_o,*xj, n, m);

  if (ar_u < ar_o)
  {
    INT u = MAX(my_u0,ar_u);
    INT o = MIN(my_o0,ar_o);
    INT offset_psij = u-ar_u;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o-u <= 2*m+1);
    assert(offset_psij+o-u <= 2*m+1);
#endif

    for (l = 0; l <= o-u; l++)
      g[u+l] += psij_const[offset_psij+l] * f;
  }
  else
  {
    INT u = MAX(my_u0,ar_u);
    INT o = my_o0;
    INT offset_psij = u-ar_u;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o-u <= 2*m+1);
    assert(offset_psij+o-u <= 2*m+1);
#endif

    for (l = 0; l <= o-u; l++)
      g[u+l] += psij_const[offset_psij+l] * f;

    u = my_u0;
    o = MIN(my_o0,ar_o);
    offset_psij += my_u0-ar_u+n;

#ifdef OMP_ASSERT
    if (u <= o)
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
      g[u+l] += psij_const[offset_psij+l] * f;
  }
}
#endif

static void nfft_trafo_1d_B(X(plan) *ths)
{
  const INT n = ths->n[0], M = ths->M_total, m = ths->m, m2p2 = 2*m+2;
  const C *g = (C*)ths->g;

  if (ths->flags & PRE_FULL_PSI)
  {
    INT k;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT l;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < m2p2; l++)
        ths->f[j] += ths->psi[j*m2p2+l] * g[ths->psi_index_g[j*m2p2+l]];
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if (ths->flags & PRE_PSI)
  {
    INT k;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      nfft_trafo_1d_compute(&ths->f[j], g, ths->psi + j * (2 * m + 2),
        &ths->x[j], n, m);
    }
    return;
  } /* if(PRE_PSI) */

  if (ths->flags & PRE_FG_PSI)
  {
    INT k;
    R fg_exp_l[m2p2];

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      const R fg_psij0 = ths->psi[2 * j], fg_psij1 = ths->psi[2 * j + 1];
      R fg_psij2 = K(1.0);
      R psij_const[m2p2];
      INT l;

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

  if (ths->flags & FG_PSI)
  {
    INT k;
    R fg_exp_l[m2p2];

    sort(ths);

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      INT u, o, l;
      R fg_psij0, fg_psij1, fg_psij2;
      R psij_const[m2p2];

      uo(ths, (INT)j, &u, &o, (INT)0);
      fg_psij0 = (PHI(ths->n[0], ths->x[j] - ((R)(u))/(R)(n), 0));
      fg_psij1 = EXP(K(2.0) * ((R)(n) * ths->x[j] - (R)(u)) / ths->b[0]);
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

  if (ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K, ip_s = K / (m + 2);
    INT k;

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R ip_y, ip_w;
      INT ip_u;
      R psij_const[m2p2];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      uo(ths, (INT)j, &u, &o, (INT)0);

      ip_y = FABS((R)(n) * ths->x[j] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);

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
    INT k;

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[m2p2];
      INT u, o, l;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      uo(ths, (INT)j, &u, &o, (INT)0);

      for (l = 0; l < m2p2; l++)
        psij_const[l] = (PHI(ths->n[0], ths->x[j] - ((R)((u+l))) / (R)(n), 0));

      nfft_trafo_1d_compute(&ths->f[j], g, psij_const, &ths->x[j], n, m);
    }
  }
}


#define MACRO_adjoint_1d_B_OMP_BLOCKWISE_COMPUTE_PRE_PSI \
{ \
            nfft_adjoint_1d_compute_omp_blockwise(ths->f[j], g, \
                ths->psi + j * (2 * m + 2), ths->x + j, n, m, my_u0, my_o0); \
}

#define MACRO_adjoint_1d_B_OMP_BLOCKWISE_COMPUTE_PRE_FG_PSI \
{ \
            R psij_const[2 * m + 2]; \
            INT l; \
            R fg_psij0 = ths->psi[2 * j]; \
            R fg_psij1 = ths->psi[2 * j + 1]; \
            R fg_psij2 = K(1.0); \
 \
            psij_const[0] = fg_psij0; \
            for (l = 1; l <= 2 * m + 1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l]; \
            } \
 \
            nfft_adjoint_1d_compute_omp_blockwise(ths->f[j], g, psij_const, \
                ths->x + j, n, m, my_u0, my_o0); \
}

#define MACRO_adjoint_1d_B_OMP_BLOCKWISE_COMPUTE_FG_PSI \
{ \
            R psij_const[2 * m + 2]; \
            R fg_psij0, fg_psij1, fg_psij2; \
            INT u, o, l; \
 \
            uo(ths, j, &u, &o, (INT)0); \
            fg_psij0 = (PHI(ths->n[0],ths->x[j]-((R)u)/((R)n),0)); \
            fg_psij1 = EXP(K(2.0) * (((R)n) * (ths->x[j]) - (R)u) / ths->b[0]); \
            fg_psij2 = K(1.0); \
            psij_const[0] = fg_psij0; \
            for (l = 1; l <= 2 * m + 1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l]; \
            } \
 \
            nfft_adjoint_1d_compute_omp_blockwise(ths->f[j], g, psij_const, \
                ths->x + j, n, m, my_u0, my_o0); \
}

#define MACRO_adjoint_1d_B_OMP_BLOCKWISE_COMPUTE_PRE_LIN_PSI \
{ \
            R psij_const[2 * m + 2]; \
            INT ip_u; \
            R ip_y, ip_w; \
            INT u, o, l; \
 \
            uo(ths, j, &u, &o, (INT)0); \
 \
            ip_y = FABS(((R)n) * ths->x[j] - (R)u) * ((R)ip_s); \
            ip_u = LRINT(FLOOR(ip_y)); \
            ip_w = ip_y - ip_u; \
            for (l = 0; l < 2 * m + 2; l++) \
              psij_const[l] \
                  = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w) \
                      + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w); \
 \
            nfft_adjoint_1d_compute_omp_blockwise(ths->f[j], g, psij_const, \
                ths->x + j, n, m, my_u0, my_o0); \
}

#define MACRO_adjoint_1d_B_OMP_BLOCKWISE_COMPUTE_NO_PSI \
{ \
            R psij_const[2 * m + 2]; \
            INT u, o, l; \
 \
            uo(ths, j, &u, &o, (INT)0); \
 \
            for (l = 0; l <= 2 * m + 1; l++) \
              psij_const[l] = (PHI(ths->n[0],ths->x[j]-((R)((u+l)))/((R)n),0)); \
 \
            nfft_adjoint_1d_compute_omp_blockwise(ths->f[j], g, psij_const, \
                ths->x + j, n, m, my_u0, my_o0); \
}

#define MACRO_adjoint_1d_B_OMP_BLOCKWISE(whichone) \
{ \
    if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT) \
    { \
      _Pragma("omp parallel private(k)") \
      { \
        INT my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b; \
        INT *ar_x = ths->index_x; \
 \
        nfft_adjoint_B_omp_blockwise_init(&my_u0, &my_o0, &min_u_a, &max_u_a, \
        		                          &min_u_b, &max_u_b, 1, &n, m); \
 \
        if (min_u_a != -1) \
        { \
          k = index_x_binary_search(ar_x, M, min_u_a); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_A \
 \
          while (k < M) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_a || u_prod > max_u_a) \
              break; \
 \
            MACRO_adjoint_1d_B_OMP_BLOCKWISE_COMPUTE_ ##whichone \
 \
            k++; \
          } \
        } \
 \
        if (min_u_b != -1) \
        { \
          k = index_x_binary_search(ar_x, M, min_u_b); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_B \
 \
          while (k < M) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_b || u_prod > max_u_b) \
              break; \
 \
            MACRO_adjoint_1d_B_OMP_BLOCKWISE_COMPUTE_ ##whichone \
 \
            k++; \
          } \
        } \
      } /* omp parallel */ \
      return; \
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */ \
}

static void nfft_adjoint_1d_B(X(plan) *ths)
{
  const INT n = ths->n[0], M = ths->M_total, m = ths->m;
  INT k;
  C *g = (C*)ths->g;

  memset(g, 0, (size_t)(ths->n_total) * sizeof(C));

  if (ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, ths->f, M,
        (INT)1, ths->n, m, ths->flags, ths->index_x);
    return;
  } /* if(PRE_FULL_PSI) */

  if (ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    MACRO_adjoint_1d_B_OMP_BLOCKWISE(PRE_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp_atomic(ths->f[j], g, ths->psi + j * (2 * m + 2), ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(ths->f + j, g, ths->psi + j * (2 * m + 2), ths->x + j, n, m);
#endif
    }

    return;
  } /* if(PRE_PSI) */

  if (ths->flags & PRE_FG_PSI)
  {
    R fg_exp_l[2 * m + 2];

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

#ifdef _OPENMP
    MACRO_adjoint_1d_B_OMP_BLOCKWISE(PRE_FG_PSI)
#endif


#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2 * m + 2];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      INT l;
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
      nfft_adjoint_1d_compute_omp_atomic(ths->f[j], g, psij_const, ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if (ths->flags & FG_PSI)
  {
    R fg_exp_l[2 * m + 2];

    nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

    sort(ths);

#ifdef _OPENMP
    MACRO_adjoint_1d_B_OMP_BLOCKWISE(FG_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      R psij_const[2 * m + 2];
      R fg_psij0, fg_psij1, fg_psij2;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      uo(ths, j, &u, &o, (INT)0);
      fg_psij0 = (PHI(ths->n[0], ths->x[j] - ((R)u) / (R)(n),0));
      fg_psij1 = EXP(K(2.0) * ((R)(n) * (ths->x[j]) - (R)(u)) / ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for (l = 1; l <= 2 * m + 1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0 * fg_psij2 * fg_exp_l[l];
      }

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp_atomic(ths->f[j], g, psij_const, ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if (ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K;
    const INT ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    MACRO_adjoint_1d_B_OMP_BLOCKWISE(PRE_LIN_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      INT ip_u;
      R ip_y, ip_w;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[2 * m + 2];

      uo(ths, j, &u, &o, (INT)0);

      ip_y = FABS((R)(n) * ths->x[j] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for (l = 0; l < 2 * m + 2; l++)
        psij_const[l]
            = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w)
                + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w);

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp_atomic(ths->f[j], g, psij_const, ths->x + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

#ifdef _OPENMP
  MACRO_adjoint_1d_B_OMP_BLOCKWISE(NO_PSI)
#endif

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT u,o,l;
    R psij_const[2 * m + 2];
    INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    uo(ths, j, &u, &o, (INT)0);

    for (l = 0; l <= 2 * m + 1; l++)
      psij_const[l] = (PHI(ths->n[0], ths->x[j] - ((R)((u+l))) / (R)(n),0));

#ifdef _OPENMP
    nfft_adjoint_1d_compute_omp_atomic(ths->f[j], g, psij_const, ths->x + j, n, m);
#else
    nfft_adjoint_1d_compute_serial(ths->f + j, g, psij_const, ths->x + j, n, m);
#endif
  }
}

void X(trafo_1d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->n[0] <= 2*ths->m+2))
  {
    X(trafo_direct)(ths);
    return;
  }
  
  const INT N = ths->N[0], N2 = N/2, n = ths->n[0];
  C *f_hat1 = (C*)ths->f_hat, *f_hat2 = (C*)&ths->f_hat[N2];

  ths->g_hat = ths->g1;
  ths->g = ths->g2;

  {
    C *g_hat1 = (C*)&ths->g_hat[n-N/2], *g_hat2 = (C*)ths->g_hat;
    R *c_phi_inv1, *c_phi_inv2;

    TIC(0)
#ifdef _OPENMP
    {
      INT k;
      #pragma omp parallel for default(shared) private(k)
      for (k = 0; k < ths->n_total; k++)
        ths->g_hat[k] = 0.0;
    }
#else
    memset(ths->g_hat, 0, (size_t)(ths->n_total) * sizeof(C));
#endif
    if(ths->flags & PRE_PHI_HUT)
    {
      INT k;
      c_phi_inv1 = ths->c_phi_inv[0];
      c_phi_inv2 = &ths->c_phi_inv[0][N2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k)
#endif
      for (k = 0; k < N2; k++)
      {
        g_hat1[k] = f_hat1[k] * c_phi_inv1[k];
        g_hat2[k] = f_hat2[k] * c_phi_inv2[k];
      }
    }
    else
    {
      INT k;
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k)
#endif
      for (k = 0; k < N2; k++)
      {
        g_hat1[k] = f_hat1[k] / (PHI_HUT(ths->n[0],k-N2,0));
        g_hat2[k] = f_hat2[k] / (PHI_HUT(ths->n[0],k,0));
      }
    }
    TOC(0)

    TIC_FFTW(1)
    FFTW(execute)(ths->my_fftw_plan1);
    TOC_FFTW(1);

    TIC(2);
    nfft_trafo_1d_B(ths);
    TOC(2);
  }
}

void X(adjoint_1d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->n[0] <= 2*ths->m+2))
  {
    X(adjoint_direct)(ths);
    return;
  }
  
  INT n,N;
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
  FFTW(execute)(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->flags & PRE_PHI_HUT)
  {
    INT k;
    c_phi_inv1=ths->c_phi_inv[0];
    c_phi_inv2=&ths->c_phi_inv[0][N/2];

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < N/2; k++)
    {
      f_hat1[k] = g_hat1[k] * c_phi_inv1[k];
      f_hat2[k] = g_hat2[k] * c_phi_inv2[k];
    }
  }
  else
  {
    INT k;

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < N/2; k++)
    {
      f_hat1[k] = g_hat1[k] / (PHI_HUT(ths->n[0],k-N/2,0));
      f_hat2[k] = g_hat2[k] / (PHI_HUT(ths->n[0],k,0));
    }
  }
  TOC(0)
}


/* ################################################ SPECIFIC VERSIONS FOR d=2 */

static void nfft_2d_init_fg_exp_l(R *fg_exp_l, const INT m, const R b)
{
  INT l;
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

static void nfft_trafo_2d_compute(C *fj, const C *g, const R *psij_const0,
    const R *psij_const1, const R *xj0, const R *xj1, const INT n0,
    const INT n1, const INT m)
{
  INT u0,o0,l0,u1,o1,l1;
  const C *gj;
  const R *psij0,*psij1;

  psij0=psij_const0;
  psij1=psij_const1;

  uo2(&u0,&o0,*xj0, n0, m);
  uo2(&u1,&o1,*xj1, n1, m);

  *fj=0;

  if (u0 < o0)
      if(u1 < o1)
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
static void nfft_adjoint_2d_compute_omp_atomic(const C f, C *g,
            const R *psij_const0, const R *psij_const1, const R *xj0,
            const R *xj1, const INT n0, const INT n1, const INT m)
{
  INT u0,o0,l0,u1,o1,l1;

  INT index_temp0[2*m+2];
  INT index_temp1[2*m+2];

  uo2(&u0,&o0,*xj0, n0, m);
  uo2(&u1,&o1,*xj1, n1, m);

  for (l0=0; l0<=2*m+1; l0++)
    index_temp0[l0] = (u0+l0)%n0;

  for (l1=0; l1<=2*m+1; l1++)
    index_temp1[l1] = (u1+l1)%n1;

  for(l0=0; l0<=2*m+1; l0++)
  {
    for(l1=0; l1<=2*m+1; l1++)
    {
      INT i = index_temp0[l0] * n1 + index_temp1[l1];
      C *lhs = g+i;
      R *lhs_real = (R*)lhs;
      C val = psij_const0[l0] * psij_const1[l1] * f;

      #pragma omp atomic
      lhs_real[0] += CREAL(val);

      #pragma omp atomic
      lhs_real[1] += CIMAG(val);
    }
  }
}
#endif

#ifdef _OPENMP
/** 
 * Adjoint NFFT for two-dimensional case updating only a specified range of
 * vector g.
 *
 * \arg f input coefficient f[j]
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
static void nfft_adjoint_2d_compute_omp_blockwise(const C f, C *g,
            const R *psij_const0, const R *psij_const1, const R *xj0,
            const R *xj1, const INT n0, const INT n1, const INT m,
            const INT my_u0, const INT my_o0)
{
  INT ar_u0,ar_o0,l0,u1,o1,l1;
  INT index_temp1[2*m+2];

  uo2(&ar_u0,&ar_o0,*xj0, n0, m);
  uo2(&u1,&o1,*xj1, n1, m);

  for (l1 = 0; l1 <= 2*m+1; l1++)
    index_temp1[l1] = (u1+l1)%n1;

  if(ar_u0 < ar_o0)
  {
    INT u0 = MAX(my_u0,ar_u0);
    INT o0 = MIN(my_o0,ar_o0);
    INT offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      INT i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
        g[i0 + index_temp1[l1]] += val0 * psij_const1[l1] * f;
    }
  }
  else
  {
    INT u0 = MAX(my_u0,ar_u0);
    INT o0 = my_o0;
    INT offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      INT i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
        g[i0 + index_temp1[l1]] += val0 * psij_const1[l1] * f;
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
      INT i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
        g[i0 + index_temp1[l1]] += val0 * psij_const1[l1] * f;
    }
  }
}
#endif

#ifndef _OPENMP
static void nfft_adjoint_2d_compute_serial(const C *fj, C *g,
            const R *psij_const0, const R *psij_const1, const R *xj0,
            const R *xj1, const INT n0, const INT n1, const INT m)
{
  INT u0,o0,l0,u1,o1,l1;
  C *gj;
  const R *psij0,*psij1;

  psij0=psij_const0;
  psij1=psij_const1;

  uo2(&u0,&o0,*xj0, n0, m);
  uo2(&u1,&o1,*xj1, n1, m);

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
#endif

static void nfft_trafo_2d_B(X(plan) *ths)
{
  const C *g = (C*)ths->g;
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT M = ths->M_total;
  const INT m = ths->m;

  INT k;

  if(ths->flags & PRE_FULL_PSI)
  {
    const INT lprod = (2*m+2) * (2*m+2);
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT l;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < lprod; l++)
        ths->f[j] += ths->psi[j*lprod+l] * g[ths->psi_index_g[j*lprod+l]];
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      nfft_trafo_2d_compute(ths->f+j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m);
    }

      return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+2)];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      INT l;
      R fg_psij0 = ths->psi[2*j*2];
      R fg_psij1 = ths->psi[2*j*2+1];
      R fg_psij2 = K(1.0);

      psij_const[0] = fg_psij0;
      for (l = 1; l <= 2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      fg_psij0 = ths->psi[2*(j*2+1)];
      fg_psij1 = ths->psi[2*(j*2+1)+1];
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for (l = 1; l <= 2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      nfft_trafo_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->flags & FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R fg_psij0, fg_psij1, fg_psij2;
      R psij_const[2*(2*m+2)];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      uo(ths, j, &u, &o, (INT)0);
      fg_psij0 = (PHI(ths->n[0], ths->x[2*j] - ((R)u) / (R)(n0),0));
      fg_psij1 = EXP(K(2.0) * ((R)(n0) * (ths->x[2*j]) - (R)(u)) / ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for (l = 1; l <= 2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      uo(ths,j,&u,&o, (INT)1);
      fg_psij0 = (PHI(ths->n[1], ths->x[2*j+1] - ((R)u) / (R)(n1),1));
      fg_psij1 = EXP(K(2.0) * ((R)(n1) * (ths->x[2*j+1]) - (R)(u)) / ths->b[1]);
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

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K, ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R ip_y, ip_w;
      INT ip_u;
      R psij_const[2*(2*m+2)];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      uo(ths,j,&u,&o,(INT)0);
      ip_y = FABS((R)(n0) * ths->x[2*j] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)LRINT(FLOOR(ip_y));
      ip_w = ip_y - (R)(ip_u);
      for (l = 0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths,j,&u,&o,(INT)1);
      ip_y = FABS((R)(n1) * ths->x[2*j+1] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for (l = 0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_trafo_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
    }
      return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */

  sort(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    R psij_const[2*(2*m+2)];
    INT u, o, l;
    INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    uo(ths,j,&u,&o,(INT)0);
    for (l = 0; l <= 2*m+1; l++)
      psij_const[l]=(PHI(ths->n[0], ths->x[2*j] - ((R)((u+l))) / (R)(n0),0));

    uo(ths,j,&u,&o,(INT)1);
    for (l = 0; l <= 2*m+1; l++)
      psij_const[2*m+2+l] = (PHI(ths->n[1], ths->x[2*j+1] - ((R)((u+l)))/(R)(n1),1));

    nfft_trafo_2d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
  }
}

#define MACRO_adjoint_2d_B_OMP_BLOCKWISE_COMPUTE_PRE_PSI \
            nfft_adjoint_2d_compute_omp_blockwise(ths->f[j], g, \
                ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), \
                ths->x+2*j, ths->x+2*j+1, n0, n1, m, my_u0, my_o0);

#define MACRO_adjoint_2d_B_OMP_BLOCKWISE_COMPUTE_PRE_FG_PSI \
{ \
            R psij_const[2*(2*m+2)]; \
            INT l; \
            R fg_psij0 = ths->psi[2*j*2]; \
            R fg_psij1 = ths->psi[2*j*2+1]; \
            R fg_psij2 = K(1.0); \
 \
            psij_const[0] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l]; \
            } \
 \
            fg_psij0 = ths->psi[2*(j*2+1)]; \
            fg_psij1 = ths->psi[2*(j*2+1)+1]; \
            fg_psij2 = K(1.0); \
            psij_const[2*m+2] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l]; \
            } \
 \
            nfft_adjoint_2d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, \
                n0, n1, m, my_u0, my_o0); \
}

#define MACRO_adjoint_2d_B_OMP_BLOCKWISE_COMPUTE_FG_PSI \
{ \
            R psij_const[2*(2*m+2)]; \
            R fg_psij0, fg_psij1, fg_psij2; \
            INT u, o, l; \
 \
            uo(ths,j,&u,&o,(INT)0); \
            fg_psij0 = (PHI(ths->n[0],ths->x[2*j]-((R)u)/((R)n0),0)); \
            fg_psij1 = EXP(K(2.0)*(((R)n0)*(ths->x[2*j]) - (R)u)/ths->b[0]); \
            fg_psij2 = K(1.0); \
            psij_const[0] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l]; \
            } \
 \
            uo(ths,j,&u,&o,(INT)1); \
            fg_psij0 = (PHI(ths->n[1],ths->x[2*j+1]-((R)u)/((R)n1),1)); \
            fg_psij1 = EXP(K(2.0)*(((R)n1)*(ths->x[2*j+1]) - (R)u)/ths->b[1]); \
            fg_psij2 = K(1.0); \
            psij_const[2*m+2] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l]; \
            } \
 \
            nfft_adjoint_2d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, \
                n0, n1, m, my_u0, my_o0); \
}

#define MACRO_adjoint_2d_B_OMP_BLOCKWISE_COMPUTE_PRE_LIN_PSI \
{ \
            R psij_const[2*(2*m+2)]; \
            INT u, o, l; \
            INT ip_u; \
            R ip_y, ip_w; \
 \
            uo(ths,j,&u,&o,(INT)0); \
            ip_y = FABS(((R)n0)*(ths->x[2*j]) - (R)u)*((R)ip_s); \
            ip_u = LRINT(FLOOR(ip_y)); \
            ip_w = ip_y-ip_u; \
            for(l=0; l < 2*m+2; l++) \
              psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + \
                ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w); \
 \
            uo(ths,j,&u,&o,(INT)1); \
            ip_y = FABS(((R)n1)*(ths->x[2*j+1]) - (R)u)*((R)ip_s); \
            ip_u = LRINT(FLOOR(ip_y)); \
            ip_w = ip_y-ip_u; \
            for(l=0; l < 2*m+2; l++) \
              psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + \
                ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w); \
 \
            nfft_adjoint_2d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, \
                n0, n1, m, my_u0, my_o0); \
}

#define MACRO_adjoint_2d_B_OMP_BLOCKWISE_COMPUTE_NO_PSI \
{ \
            R psij_const[2*(2*m+2)]; \
            INT u, o, l; \
 \
            uo(ths,j,&u,&o,(INT)0); \
            for(l=0;l<=2*m+1;l++) \
              psij_const[l]=(PHI(ths->n[0],ths->x[2*j]-((R)((u+l)))/((R)n0),0)); \
 \
            uo(ths,j,&u,&o,(INT)1); \
            for(l=0;l<=2*m+1;l++) \
              psij_const[2*m+2+l]=(PHI(ths->n[1],ths->x[2*j+1]-((R)((u+l)))/((R)n1),1)); \
 \
            nfft_adjoint_2d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, \
                n0, n1, m, my_u0, my_o0); \
}

#define MACRO_adjoint_2d_B_OMP_BLOCKWISE(whichone) \
{ \
    if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT) \
    { \
      _Pragma("omp parallel private(k)") \
      { \
        INT my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b; \
        INT *ar_x = ths->index_x; \
 \
        nfft_adjoint_B_omp_blockwise_init(&my_u0, &my_o0, &min_u_a, &max_u_a, \
            &min_u_b, &max_u_b, 2, ths->n, m); \
 \
        if (min_u_a != -1) \
        { \
          k = index_x_binary_search(ar_x, M, min_u_a); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_A \
 \
          while (k < M) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_a || u_prod > max_u_a) \
              break; \
 \
            MACRO_adjoint_2d_B_OMP_BLOCKWISE_COMPUTE_ ##whichone \
 \
            k++; \
          } \
        } \
 \
        if (min_u_b != -1) \
        { \
          INT k = index_x_binary_search(ar_x, M, min_u_b); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_B \
 \
          while (k < M) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_b || u_prod > max_u_b) \
              break; \
 \
            MACRO_adjoint_2d_B_OMP_BLOCKWISE_COMPUTE_ ##whichone \
 \
            k++; \
          } \
        } \
      } /* omp parallel */ \
      return; \
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */ \
}


static void nfft_adjoint_2d_B(X(plan) *ths)
{
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT M = ths->M_total;
  const INT m = ths->m;
  C* g = (C*) ths->g;
  INT k;

  memset(g, 0, (size_t)(ths->n_total) * sizeof(C));

  if(ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, ths->f, M,
        (INT)2, ths->n, m, ths->flags, ths->index_x);
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    MACRO_adjoint_2d_B_OMP_BLOCKWISE(PRE_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp_atomic(ths->f[j], g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(ths->f+j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

#ifdef _OPENMP
    MACRO_adjoint_2d_B_OMP_BLOCKWISE(PRE_FG_PSI)
#endif


#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+2)];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      INT l;
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
      nfft_adjoint_2d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->flags & FG_PSI)
  {
    R fg_exp_l[2*(2*m+2)];

    nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

    sort(ths);

#ifdef _OPENMP
    MACRO_adjoint_2d_B_OMP_BLOCKWISE(FG_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R fg_psij0, fg_psij1, fg_psij2;
      R psij_const[2*(2*m+2)];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      uo(ths,j,&u,&o,(INT)0);
      fg_psij0 = (PHI(ths->n[0], ths->x[2*j] - ((R)u)/(R)(n0),0));
      fg_psij1 = EXP(K(2.0) * ((R)(n0) * (ths->x[2*j]) - (R)(u)) / ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      uo(ths,j,&u,&o,(INT)1);
      fg_psij0 = (PHI(ths->n[1], ths->x[2*j+1] - ((R)u) / (R)(n1),1));
      fg_psij1 = EXP(K(2.0) * ((R)(n1) * (ths->x[2*j+1]) - (R)(u)) / ths->b[1]);
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K;
    const INT ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    MACRO_adjoint_2d_B_OMP_BLOCKWISE(PRE_LIN_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      INT ip_u;
      R ip_y, ip_w;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[2*(2*m+2)];

      uo(ths,j,&u,&o,(INT)0);
      ip_y = FABS((R)(n0) * (ths->x[2*j]) - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths,j,&u,&o,(INT)1);
      ip_y = FABS((R)(n1) * (ths->x[2*j+1]) - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
  }
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

#ifdef _OPENMP
  MACRO_adjoint_2d_B_OMP_BLOCKWISE(NO_PSI)
#endif

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT u,o,l;
    R psij_const[2*(2*m+2)];
    INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    uo(ths,j,&u,&o,(INT)0);
    for(l=0;l<=2*m+1;l++)
      psij_const[l]=(PHI(ths->n[0], ths->x[2*j] - ((R)((u+l))) / (R)(n0),0));

    uo(ths,j,&u,&o,(INT)1);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*m+2+l]=(PHI(ths->n[1], ths->x[2*j+1] - ((R)((u+l))) / (R)(n1),1));

#ifdef _OPENMP
    nfft_adjoint_2d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#else
    nfft_adjoint_2d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, ths->x+2*j, ths->x+2*j+1, n0, n1, m);
#endif
  }
}


void X(trafo_2d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2))
  {
    X(trafo_direct)(ths);
    return;
  }
  
  INT k0,k1,n0,n1,N0,N1;
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
  memset(ths->g_hat, 0, (size_t)(ths->n_total) * sizeof(C));
#endif
  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,c_phi_inv11,c_phi_inv12,g_hat11,f_hat11,g_hat21,f_hat21,g_hat12,f_hat12,g_hat22,f_hat22,ck11,ck12)
#endif
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
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,ck11,ck12)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));
      g_hat[(n0-N0/2+k0)*n1+n1-N1/2+k1] = f_hat[k0*N1+k1]             * ck01 * ck11;
      g_hat[k0*n1+n1-N1/2+k1]           = f_hat[(N0/2+k0)*N1+k1]      * ck02 * ck11;
      g_hat[(n0-N0/2+k0)*n1+k1]         = f_hat[k0*N1+N1/2+k1]        * ck01 * ck12;
      g_hat[k0*n1+k1]                   = f_hat[(N0/2+k0)*N1+N1/2+k1] * ck02 * ck12;
    }
      }

  TOC(0)

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan1);
  TOC_FFTW(1);

  TIC(2);
  nfft_trafo_2d_B(ths);
  TOC(2);
}

void X(adjoint_2d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2))
  {
    X(adjoint_direct)(ths);
    return;
  }
  
  INT k0,k1,n0,n1,N0,N1;
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
  FFTW(execute)(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,c_phi_inv11,c_phi_inv12,g_hat11,f_hat11,g_hat21,f_hat21,g_hat12,f_hat12,g_hat22,f_hat22,ck11,ck12)
#endif
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
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,ck11,ck12)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));
      f_hat[k0*N1+k1]             = g_hat[(n0-N0/2+k0)*n1+n1-N1/2+k1] * ck01 * ck11;
      f_hat[(N0/2+k0)*N1+k1]      = g_hat[k0*n1+n1-N1/2+k1]           * ck02 * ck11;
      f_hat[k0*N1+N1/2+k1]        = g_hat[(n0-N0/2+k0)*n1+k1]         * ck01 * ck12;
      f_hat[(N0/2+k0)*N1+N1/2+k1] = g_hat[k0*n1+k1]                   * ck02 * ck12;
    }
      }
  TOC(0)
}

/* ################################################ SPECIFIC VERSIONS FOR d=3 */

static void nfft_3d_init_fg_exp_l(R *fg_exp_l, const INT m, const R b)
{
  INT l;
  R fg_exp_b0, fg_exp_b1, fg_exp_b2, fg_exp_b0_sq;

  fg_exp_b0 = EXP(-K(1.0) / b);
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

static void nfft_trafo_3d_compute(C *fj, const C *g, const R *psij_const0,
    const R *psij_const1, const R *psij_const2, const R *xj0, const R *xj1,
    const R *xj2, const INT n0, const INT n1, const INT n2, const INT m)
{
  INT u0, o0, l0, u1, o1, l1, u2, o2, l2;
  const C *gj;
  const R *psij0, *psij1, *psij2;

  psij0 = psij_const0;
  psij1 = psij_const1;
  psij2 = psij_const2;

  uo2(&u0, &o0, *xj0, n0, m);
  uo2(&u1, &o1, *xj1, n1, m);
  uo2(&u2, &o2, *xj2, n2, m);

  *fj = 0;

  if (u0 < o0)
    if (u1 < o1)
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      else
        /* asserts (u2>o2)*/
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
    else /* asserts (u1>o1)*/
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      else/* asserts (u2>o2) */
      {
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      }
  else /* asserts (u0>o0) */
    if (u1 < o1)
      if (u2 < o2)
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      } else/* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      }
    else /* asserts (u1>o1) */
      if (u2 < o2)
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      } else/* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
            gj = g + (l0 * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*fj) += (*psij0) * (*psij1) * (*psij2++) * (*gj++);
          }
        }
      }
}

#ifdef _OPENMP
/** 
 * Adjoint NFFT for three-dimensional case updating only a specified range of
 * vector g.
 *
 * \arg f input coefficient f[j]
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
static void nfft_adjoint_3d_compute_omp_blockwise(const C f, C *g,
    const R *psij_const0, const R *psij_const1, const R *psij_const2,
    const R *xj0, const R *xj1, const R *xj2,
    const INT n0, const INT n1, const INT n2, const INT m,
    const INT my_u0, const INT my_o0)
{
  INT ar_u0,ar_o0,l0,u1,o1,l1,u2,o2,l2;

  INT index_temp1[2*m+2];
  INT index_temp2[2*m+2];

  uo2(&ar_u0,&ar_o0,*xj0, n0, m);
  uo2(&u1,&o1,*xj1, n1, m);
  uo2(&u2,&o2,*xj2, n2, m);

  for (l1=0; l1<=2*m+1; l1++)
    index_temp1[l1] = (u1+l1)%n1;

  for (l2=0; l2<=2*m+1; l2++)
    index_temp2[l2] = (u2+l2)%n2;

  if(ar_u0<ar_o0)
  {
    INT u0 = MAX(my_u0,ar_u0);
    INT o0 = MIN(my_o0,ar_o0);
    INT offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      const INT i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
      {
        const INT i1 = (i0 + index_temp1[l1]) * n2;
        const C val1 = psij_const1[l1];

        for(l2=0; l2<=2*m+1; l2++)
          g[i1 + index_temp2[l2]] += val0 * val1 * psij_const2[l2] * f;
      }
    }
  }
  else
  {
    INT u0 = MAX(my_u0,ar_u0);
    INT o0 = my_o0;
    INT offset_psij = u0-ar_u0;
#ifdef OMP_ASSERT
    assert(offset_psij >= 0);
    assert(o0-u0 <= 2*m+1);
    assert(offset_psij+o0-u0 <= 2*m+1);
#endif

    for (l0 = 0; l0 <= o0-u0; l0++)
    {
      INT i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
      {
        const INT i1 = (i0 + index_temp1[l1]) * n2;
        const C val1 = psij_const1[l1];

        for(l2=0; l2<=2*m+1; l2++)
          g[i1 + index_temp2[l2]] += val0 * val1 * psij_const2[l2] * f;
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
      INT i0 = (u0+l0) * n1;
      const C val0 = psij_const0[offset_psij+l0];

      for(l1=0; l1<=2*m+1; l1++)
      {
        const INT i1 = (i0 + index_temp1[l1]) * n2;
        const C val1 = psij_const1[l1];

        for(l2=0; l2<=2*m+1; l2++)
          g[i1 + index_temp2[l2]] += val0 * val1 * psij_const2[l2] * f;
      }
    }
  }
}
#endif

#ifdef _OPENMP
/* adjoint NFFT three-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_3d_compute_omp_atomic(const C f, C *g,
    const R *psij_const0, const R *psij_const1, const R *psij_const2,
    const R *xj0, const R *xj1, const R *xj2,
    const INT n0, const INT n1, const INT n2, const INT m)
{
  INT u0,o0,l0,u1,o1,l1,u2,o2,l2;

  INT index_temp0[2*m+2];
  INT index_temp1[2*m+2];
  INT index_temp2[2*m+2];

  uo2(&u0,&o0,*xj0, n0, m);
  uo2(&u1,&o1,*xj1, n1, m);
  uo2(&u2,&o2,*xj2, n2, m);

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
        INT i = (index_temp0[l0] * n1 + index_temp1[l1]) * n2 + index_temp2[l2];
        C *lhs = g+i;
        R *lhs_real = (R*)lhs;
        C val = psij_const0[l0] * psij_const1[l1] * psij_const2[l2] * f;

#pragma omp atomic
        lhs_real[0] += CREAL(val);

#pragma omp atomic
        lhs_real[1] += CIMAG(val);
      }
    }
  }
}
#endif

#ifndef _OPENMP
static void nfft_adjoint_3d_compute_serial(const C *fj, C *g,
    const R *psij_const0, const R *psij_const1, const R *psij_const2, const R *xj0,
    const R *xj1, const R *xj2, const INT n0, const INT n1, const INT n2,
    const INT m)
{
  INT u0, o0, l0, u1, o1, l1, u2, o2, l2;
  C *gj;
  const R *psij0, *psij1, *psij2;

  psij0 = psij_const0;
  psij1 = psij_const1;
  psij2 = psij_const2;

  uo2(&u0, &o0, *xj0, n0, m);
  uo2(&u1, &o1, *xj1, n1, m);
  uo2(&u2, &o2, *xj2, n2, m);

  if (u0 < o0)
    if (u1 < o1)
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      else
        /* asserts (u2>o2)*/
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
    else /* asserts (u1>o1)*/
      if (u2 < o2)
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      else/* asserts (u2>o2) */
      {
        for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      }
  else /* asserts (u0>o0) */
    if (u1 < o1)
      if (u2 < o2)
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      } else/* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 <= 2 * m + 1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      }
    else /* asserts (u1>o1) */
      if (u2 < o2)
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 <= 2 * m + 1; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      } else/* asserts (u2>o2) */
      {
        for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + ((u0 + l0) * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + ((u0 + l0) * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }

        for (l0 = 0; l0 <= o0; l0++, psij0++)
        {
          psij1 = psij_const1;
          for (l1 = 0; l1 < 2 * m + 1 - o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + (u1 + l1)) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + (l0 * n1 + (u1 + l1)) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
          for (l1 = 0; l1 <= o1; l1++, psij1++)
          {
            psij2 = psij_const2;
            gj = g + (l0 * n1 + l1) * n2 + u2;
            for (l2 = 0; l2 < 2 * m + 1 - o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
            gj = g + (l0 * n1 + l1) * n2;
            for (l2 = 0; l2 <= o2; l2++)
              (*gj++) += (*psij0) * (*psij1) * (*psij2++) * (*fj);
          }
        }
      }
}
#endif

static void nfft_trafo_3d_B(X(plan) *ths)
{
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT n2 = ths->n[2];
  const INT M = ths->M_total;
  const INT m = ths->m;

  const C* g = (C*) ths->g;

  INT k;

  if(ths->flags & PRE_FULL_PSI)
  {
    const INT lprod = (2*m+2) * (2*m+2) * (2*m+2);
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT l;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      ths->f[j] = K(0.0);
      for (l = 0; l < lprod; l++)
        ths->f[j] += ths->psi[j*lprod+l] * g[ths->psi_index_g[j*lprod+l]];
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      nfft_trafo_3d_compute(ths->f+j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      INT l;
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

  if(ths->flags & FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      INT u, o, l;
      R psij_const[3*(2*m+2)];
      R fg_psij0, fg_psij1, fg_psij2;

      uo(ths,j,&u,&o,(INT)0);
      fg_psij0 = (PHI(ths->n[0], ths->x[3*j] - ((R)u) / (R)(n0),0));
      fg_psij1 = EXP(K(2.0) * ((R)(n0) * (ths->x[3*j]) - (R)(u)) / ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      uo(ths,j,&u,&o,(INT)1);
      fg_psij0 = (PHI(ths->n[1], ths->x[3*j+1] - ((R)u) / (R)(n1),1));
      fg_psij1 = EXP(K(2.0) * ((R)(n1) * (ths->x[3*j+1]) - (R)(u)) / ths->b[1]);
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      uo(ths,j,&u,&o,(INT)2);
      fg_psij0 = (PHI(ths->n[2], ths->x[3*j+2] - ((R)u) / (R)(n2),2));
      fg_psij1 = EXP(K(2.0) * ((R)(n2) * (ths->x[3*j+2]) - (R)(u)) / ths->b[2]);
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

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K, ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R ip_y, ip_w;
      INT ip_u;
      R psij_const[3*(2*m+2)];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

      uo(ths,j,&u,&o,(INT)0);
      ip_y = FABS((R)(n0) * ths->x[3*j+0] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths,j,&u,&o,(INT)1);
      ip_y = FABS((R)(n1) * ths->x[3*j+1] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths,j,&u,&o,(INT)2);
      ip_y = FABS((R)(n2) * ths->x[3*j+2] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_trafo_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */

  sort(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    R psij_const[3*(2*m+2)];
    INT u, o, l;
    INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    uo(ths,j,&u,&o,(INT)0);
    for(l=0;l<=2*m+1;l++)
      psij_const[l]=(PHI(ths->n[0], ths->x[3*j] - ((R)((u+l))) / (R)(n0),0));

    uo(ths,j,&u,&o,(INT)1);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*m+2+l]=(PHI(ths->n[1], ths->x[3*j+1] - ((R)((u+l))) / (R)(n1),1));

    uo(ths,j,&u,&o,(INT)2);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*(2*m+2)+l]=(PHI(ths->n[2], ths->x[3*j+2] - ((R)((u+l))) / (R)(n2),2));

    nfft_trafo_3d_compute(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
  }
}

#define MACRO_adjoint_3d_B_OMP_BLOCKWISE_COMPUTE_PRE_PSI \
            nfft_adjoint_3d_compute_omp_blockwise(ths->f[j], g, \
                ths->psi+j*3*(2*m+2), \
                ths->psi+(j*3+1)*(2*m+2), \
                ths->psi+(j*3+2)*(2*m+2), \
                ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, \
                n0, n1, n2, m, my_u0, my_o0);

#define MACRO_adjoint_3d_B_OMP_BLOCKWISE_COMPUTE_PRE_FG_PSI \
{ \
            INT l; \
            R psij_const[3*(2*m+2)]; \
            R fg_psij0 = ths->psi[2*j*3]; \
            R fg_psij1 = ths->psi[2*j*3+1]; \
            R fg_psij2 = K(1.0); \
 \
            psij_const[0] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l]; \
            } \
 \
            fg_psij0 = ths->psi[2*(j*3+1)]; \
            fg_psij1 = ths->psi[2*(j*3+1)+1]; \
            fg_psij2 = K(1.0); \
            psij_const[2*m+2] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l]; \
            } \
 \
            fg_psij0 = ths->psi[2*(j*3+2)]; \
            fg_psij1 = ths->psi[2*(j*3+2)+1]; \
            fg_psij2 = K(1.0); \
            psij_const[2*(2*m+2)] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l]; \
            } \
 \
            nfft_adjoint_3d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, \
                ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, \
                n0, n1, n2, m, my_u0, my_o0); \
}

#define MACRO_adjoint_3d_B_OMP_BLOCKWISE_COMPUTE_FG_PSI \
{ \
            INT u, o, l; \
            R psij_const[3*(2*m+2)]; \
            R fg_psij0, fg_psij1, fg_psij2; \
 \
            uo(ths,j,&u,&o,(INT)0); \
            fg_psij0 = (PHI(ths->n[0],ths->x[3*j]-((R)u)/((R)n0),0)); \
            fg_psij1 = EXP(K(2.0)*(((R)n0)*(ths->x[3*j]) - (R)u)/ths->b[0]); \
            fg_psij2 = K(1.0); \
            psij_const[0] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l]; \
            } \
 \
            uo(ths,j,&u,&o,(INT)1); \
            fg_psij0 = (PHI(ths->n[1],ths->x[3*j+1]-((R)u)/((R)n1),1)); \
            fg_psij1 = EXP(K(2.0)*(((R)n1)*(ths->x[3*j+1]) - (R)u)/ths->b[1]); \
            fg_psij2 = K(1.0); \
            psij_const[2*m+2] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l]; \
            } \
 \
            uo(ths,j,&u,&o,(INT)2); \
            fg_psij0 = (PHI(ths->n[2],ths->x[3*j+2]-((R)u)/((R)n2),2)); \
            fg_psij1 = EXP(K(2.0)*(((R)n2)*(ths->x[3*j+2]) - (R)u)/ths->b[2]); \
            fg_psij2 = K(1.0); \
            psij_const[2*(2*m+2)] = fg_psij0; \
            for(l=1; l<=2*m+1; l++) \
            { \
              fg_psij2 *= fg_psij1; \
              psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l]; \
            } \
 \
            nfft_adjoint_3d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, \
                ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, \
                n0, n1, n2, m, my_u0, my_o0); \
}

#define MACRO_adjoint_3d_B_OMP_BLOCKWISE_COMPUTE_PRE_LIN_PSI \
{ \
            INT u, o, l; \
            R psij_const[3*(2*m+2)]; \
            INT ip_u; \
            R ip_y, ip_w; \
 \
            uo(ths,j,&u,&o,(INT)0); \
            ip_y = FABS(((R)n0)*ths->x[3*j+0] - (R)u)*((R)ip_s); \
            ip_u = LRINT(FLOOR(ip_y)); \
            ip_w = ip_y-ip_u; \
            for(l=0; l < 2*m+2; l++) \
              psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + \
                ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w); \
 \
            uo(ths,j,&u,&o,(INT)1); \
            ip_y = FABS(((R)n1)*ths->x[3*j+1] - (R)u)*((R)ip_s); \
            ip_u = LRINT(FLOOR(ip_y)); \
            ip_w = ip_y-ip_u; \
            for(l=0; l < 2*m+2; l++) \
              psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + \
                ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w); \
 \
            uo(ths,j,&u,&o,(INT)2); \
            ip_y = FABS(((R)n2)*ths->x[3*j+2] - (R)u)*((R)ip_s); \
            ip_u = LRINT(FLOOR(ip_y)); \
            ip_w = ip_y-ip_u; \
            for(l=0; l < 2*m+2; l++) \
              psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + \
                ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w); \
 \
            nfft_adjoint_3d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, \
                ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, \
                n0, n1, n2, m, my_u0, my_o0); \
}

#define MACRO_adjoint_3d_B_OMP_BLOCKWISE_COMPUTE_NO_PSI \
{ \
            INT u, o, l; \
            R psij_const[3*(2*m+2)]; \
 \
            uo(ths,j,&u,&o,(INT)0); \
            for(l=0;l<=2*m+1;l++) \
              psij_const[l]=(PHI(ths->n[0],ths->x[3*j]-((R)((u+l)))/((R) n0),0)); \
 \
            uo(ths,j,&u,&o,(INT)1); \
            for(l=0;l<=2*m+1;l++) \
              psij_const[2*m+2+l]=(PHI(ths->n[1],ths->x[3*j+1]-((R)((u+l)))/((R) n1),1)); \
 \
            uo(ths,j,&u,&o,(INT)2); \
            for(l=0;l<=2*m+1;l++) \
              psij_const[2*(2*m+2)+l]=(PHI(ths->n[2],ths->x[3*j+2]-((R)((u+l)))/((R) n2),2)); \
 \
            nfft_adjoint_3d_compute_omp_blockwise(ths->f[j], g, \
                psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, \
                ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, \
                n0, n1, n2, m, my_u0, my_o0); \
}

#define MACRO_adjoint_3d_B_OMP_BLOCKWISE(whichone) \
{ \
    if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT) \
    { \
      _Pragma("omp parallel private(k)") \
      { \
        INT my_u0, my_o0, min_u_a, max_u_a, min_u_b, max_u_b; \
        INT *ar_x = ths->index_x; \
 \
        nfft_adjoint_B_omp_blockwise_init(&my_u0, &my_o0, &min_u_a, &max_u_a, \
            &min_u_b, &max_u_b, 3, ths->n, m); \
 \
        if (min_u_a != -1) \
        { \
          k = index_x_binary_search(ar_x, M, min_u_a); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_A \
 \
          while (k < M) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_a || u_prod > max_u_a) \
              break; \
 \
            MACRO_adjoint_3d_B_OMP_BLOCKWISE_COMPUTE_ ##whichone \
 \
            k++; \
          } \
        } \
 \
        if (min_u_b != -1) \
        { \
          INT k = index_x_binary_search(ar_x, M, min_u_b); \
 \
          MACRO_adjoint_nd_B_OMP_BLOCKWISE_ASSERT_B \
 \
          while (k < M) \
          { \
            INT u_prod = ar_x[2*k]; \
            INT j = ar_x[2*k+1]; \
 \
            if (u_prod < min_u_b || u_prod > max_u_b) \
              break; \
 \
            MACRO_adjoint_3d_B_OMP_BLOCKWISE_COMPUTE_ ##whichone \
 \
            k++; \
          } \
        } \
      } /* omp parallel */ \
      return; \
    } /* if(NFFT_OMP_BLOCKWISE_ADJOINT) */ \
}

static void nfft_adjoint_3d_B(X(plan) *ths)
{
  INT k;
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT n2 = ths->n[2];
  const INT M = ths->M_total;
  const INT m = ths->m;

  C* g = (C*) ths->g;

  memset(g, 0, (size_t)(ths->n_total) * sizeof(C));

  if(ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, ths->f, M,
        (INT)3, ths->n, m, ths->flags, ths->index_x);
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    MACRO_adjoint_3d_B_OMP_BLOCKWISE(PRE_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp_atomic(ths->f[j], g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(ths->f+j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

#ifdef _OPENMP
    MACRO_adjoint_3d_B_OMP_BLOCKWISE(PRE_FG_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[3*(2*m+2)];
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      INT l;
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
      nfft_adjoint_3d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->flags & FG_PSI)
  {
    R fg_exp_l[3*(2*m+2)];

    nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
    nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

    sort(ths);

#ifdef _OPENMP
    MACRO_adjoint_3d_B_OMP_BLOCKWISE(FG_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[3*(2*m+2)];
      R fg_psij0, fg_psij1, fg_psij2;

      uo(ths,j,&u,&o,(INT)0);
      fg_psij0 = (PHI(ths->n[0], ths->x[3*j] - ((R)u) / (R)(n0),0));
      fg_psij1 = EXP(K(2.0) * ((R)(n0) * (ths->x[3*j]) - (R)(u))/ths->b[0]);
      fg_psij2 = K(1.0);
      psij_const[0] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
      }

      uo(ths,j,&u,&o,(INT)1);
      fg_psij0 = (PHI(ths->n[1], ths->x[3*j+1] - ((R)u) / (R)(n1),1));
      fg_psij1 = EXP(K(2.0) * ((R)(n1) * (ths->x[3*j+1]) - (R)(u))/ths->b[1]);
      fg_psij2 = K(1.0);
      psij_const[2*m+2] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
      }

      uo(ths,j,&u,&o,(INT)2);
      fg_psij0 = (PHI(ths->n[2], ths->x[3*j+2] - ((R)u) / (R)(n2),2));
      fg_psij1 = EXP(K(2.0) * ((R)(n2) * (ths->x[3*j+2]) - (R)(u))/ths->b[2]);
      fg_psij2 = K(1.0);
      psij_const[2*(2*m+2)] = fg_psij0;
      for(l=1; l<=2*m+1; l++)
      {
        fg_psij2 *= fg_psij1;
        psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
      }

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K;
    const INT ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    MACRO_adjoint_3d_B_OMP_BLOCKWISE(PRE_LIN_PSI)
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      INT ip_u;
      R ip_y, ip_w;
      INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;
      R psij_const[3*(2*m+2)];

      uo(ths,j,&u,&o,(INT)0);
      ip_y = FABS((R)(n0) * ths->x[3*j+0] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths,j,&u,&o,(INT)1);
      ip_y = FABS((R)(n1) * ths->x[3*j+1] - (R)(u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths,j,&u,&o,(INT)2);
      ip_y = FABS((R)(n2) * ths->x[3*j+2] - (R)(u))*((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

#ifdef _OPENMP
  MACRO_adjoint_3d_B_OMP_BLOCKWISE(NO_PSI)
#endif

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT u,o,l;
    R psij_const[3*(2*m+2)];
    INT j = (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2*k+1] : k;

    uo(ths,j,&u,&o,(INT)0);
    for(l=0;l<=2*m+1;l++)
      psij_const[l]=(PHI(ths->n[0], ths->x[3*j] - ((R)((u+l))) / (R)(n0),0));

    uo(ths,j,&u,&o,(INT)1);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*m+2+l]=(PHI(ths->n[1], ths->x[3*j+1] - ((R)((u+l))) / (R)(n1),1));

    uo(ths,j,&u,&o,(INT)2);
    for(l=0;l<=2*m+1;l++)
      psij_const[2*(2*m+2)+l]=(PHI(ths->n[2], ths->x[3*j+2] - ((R)((u+l))) / (R)(n2),2));

#ifdef _OPENMP
    nfft_adjoint_3d_compute_omp_atomic(ths->f[j], g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#else
    nfft_adjoint_3d_compute_serial(ths->f+j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, ths->x+3*j, ths->x+3*j+1, ths->x+3*j+2, n0, n1, n2, m);
#endif
  }
}


void X(trafo_3d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->N[2] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2) || (ths->n[2] <= 2*ths->m+2))
  {
    X(trafo_direct)(ths);
    return;
  }
  
  INT k0,k1,k2,n0,n1,n2,N0,N1,N2;
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
  memset(ths->g_hat, 0, (size_t)(ths->n_total) * sizeof(C));
#endif

  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,c_phi_inv11,c_phi_inv12,ck11,ck12,c_phi_inv21,c_phi_inv22,g_hat111,f_hat111,g_hat211,f_hat211,g_hat121,f_hat121,g_hat221,f_hat221,g_hat112,f_hat112,g_hat212,f_hat212,g_hat122,f_hat122,g_hat222,f_hat222,ck21,ck22)
#endif
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
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,ck11,ck12,ck21,ck22)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));

      for(k2=0;k2<N2/2;k2++)
        {
    ck21=K(1.0)/(PHI_HUT(ths->n[2],k2-N2/2,2));
    ck22=K(1.0)/(PHI_HUT(ths->n[2],k2,2));

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
  FFTW(execute)(ths->my_fftw_plan1);
  TOC_FFTW(1);

  TIC(2);
  nfft_trafo_3d_B(ths);
  TOC(2);
}

void X(adjoint_3d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->N[2] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2) || (ths->n[2] <= 2*ths->m+2))
  {
    X(adjoint_direct)(ths);
    return;
  }
  
  INT k0,k1,k2,n0,n1,n2,N0,N1,N2;
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
  FFTW(execute)(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,c_phi_inv11,c_phi_inv12,ck11,ck12,c_phi_inv21,c_phi_inv22,g_hat111,f_hat111,g_hat211,f_hat211,g_hat121,f_hat121,g_hat221,f_hat221,g_hat112,f_hat112,g_hat212,f_hat212,g_hat122,f_hat122,g_hat222,f_hat222,ck21,ck22)
#endif
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
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,ck11,ck12,ck21,ck22)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));

      for(k2=0;k2<N2/2;k2++)
        {
    ck21=K(1.0)/(PHI_HUT(ths->n[2],k2-N2/2,2));
    ck22=K(1.0)/(PHI_HUT(ths->n[2],k2,2));

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
void X(trafo)(X(plan) *ths)
{
  /* use direct transform if degree N is too low */
  for (int j = 0; j < ths->d; j++)
  {
    if((ths->N[j] <= ths->m) || (ths->n[j] <= 2*ths->m+2))
    {
      X(trafo_direct)(ths);
      return;
    }
  }
  
  switch(ths->d)
  {
    case 1: X(trafo_1d)(ths); break;
    case 2: X(trafo_2d)(ths); break;
    case 3: X(trafo_3d)(ths); break;
    default:
    {
      /* use ths->my_fftw_plan1 */
      ths->g_hat = ths->g1;
      ths->g = ths->g2;

      /** form \f$ \hat g_k = \frac{\hat f_k}{c_k\left(\phi\right)} \text{ for }
       *  k \in I_N \f$
       */
      TIC(0)
      D_A(ths);
      TOC(0)

      /** compute by d-variate discrete Fourier transform
       *  \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
       *  \text{ for } l \in I_n \f$
       */
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_plan1);
      TOC_FFTW(1)

      /** set \f$ f_j =\sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
       *  \text{ for } j=0,\dots,M_total-1 \f$
       */
      TIC(2)
      B_A(ths);
      TOC(2)
    }
  }
} /* nfft_trafo */

void X(adjoint)(X(plan) *ths)
{
  /* use direct transform if degree N is too low */
  for (int j = 0; j < ths->d; j++)
  {
    if((ths->N[j] <= ths->m) || (ths->n[j] <= 2*ths->m+2))
    {
      X(adjoint_direct)(ths);
      return;
    }
  }
  
  switch(ths->d)
  {
    case 1: X(adjoint_1d)(ths); break;
    case 2: X(adjoint_2d)(ths); break;
    case 3: X(adjoint_3d)(ths); break;
    default:
    {
      /* use ths->my_fftw_plan2 */
      ths->g_hat=ths->g1;
      ths->g=ths->g2;
      
      /** set \f$ g_l = \sum_{j=0}^{M_total-1} f_j \psi\left(x_j-\frac{l}{n}\right)
       *  \text{ for } l \in I_n,m(x_j) \f$
       */
      TIC(2)
      B_T(ths);
      TOC(2)

      /** compute by d-variate discrete Fourier transform
       *  \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{+2\pi {\rm i} \frac{kl}{n}}
       *  \text{ for }  k \in I_N\f$
       */
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_plan2);
      TOC_FFTW(1)

      /** form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
       *  k \in I_N \f$
       */
      TIC(0)
      D_T(ths);
      TOC(0)
    }
  }
} /* nfft_adjoint */


/** initialisation of direct transform
 */
static void precompute_phi_hut(X(plan) *ths)
{
  INT ks[ths->d]; /* index over all frequencies */
  INT t; /* index over all dimensions */

  ths->c_phi_inv = (R**) Y(malloc)((size_t)(ths->d) * sizeof(R*));

  for (t = 0; t < ths->d; t++)
  {
    ths->c_phi_inv[t] = (R*)Y(malloc)((size_t)(ths->N[t]) * sizeof(R));

    for (ks[t] = 0; ks[t] < ths->N[t]; ks[t]++)
    {
      ths->c_phi_inv[t][ks[t]]= K(1.0) / (PHI_HUT(ths->n[t], ks[t] - ths->N[t] / 2,t));
    }
  }
} /* nfft_phi_hut */

/** create a lookup table, but NOT for each node
 *  good idea K=2^xx
 *  assumes an EVEN window function
 */
void X(precompute_lin_psi)(X(plan) *ths)
{
  INT t;                                /**< index over all dimensions       */
  INT j;                                /**< index over all nodes            */
  R step;                          /**< step size in [0,(m+2)/n]        */

  for (t=0; t<ths->d; t++)
    {
      step = ((R)(ths->m+2)) / ((R)(ths->K * ths->n[t]));
      for(j = 0;j <= ths->K; j++)
  {
    ths->psi[(ths->K+1)*t + j] = PHI(ths->n[t], (R)(j) * step,t);
  } /* for(j) */
    } /* for(t) */
}

void X(precompute_fg_psi)(X(plan) *ths)
{
  INT t;                                /**< index over all dimensions       */
  INT u, o;                             /**< depends on x_j                  */

  sort(ths);

  for (t=0; t<ths->d; t++)
  {
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j,u,o)
#endif
    for (j = 0; j < ths->M_total; j++)
      {
  uo(ths,j,&u,&o,t);

        ths->psi[2*(j*ths->d+t)]=
            (PHI(ths->n[t] ,(ths->x[j*ths->d+t] - ((R)u) / (R)(ths->n[t])),t));

        ths->psi[2*(j*ths->d+t)+1]=
            EXP(K(2.0) * ((R)(ths->n[t]) * ths->x[j*ths->d+t] - (R)(u)) / ths->b[t]);
      } /* for(j) */
  }
  /* for(t) */
} /* nfft_precompute_fg_psi */

void X(precompute_psi)(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT l; /* index u<=l<=o */
  INT lj; /* index 0<=lj<u+o+1 */
  INT u, o; /* depends on x_j */

  sort(ths);

  for (t=0; t<ths->d; t++)
  {
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j,l,lj,u,o)
#endif
    for (j = 0; j < ths->M_total; j++)
    {
      uo(ths,j,&u,&o,t);

      for(l = u, lj = 0; l <= o; l++, lj++)
        ths->psi[(j * ths->d + t) * (2 * ths->m + 2) + lj] =
            (PHI(ths->n[t], (ths->x[j*ths->d+t] - ((R)l) / (R)(ths->n[t])), t));
    } /* for(j) */
  }
  /* for(t) */
} /* nfft_precompute_psi */

#ifdef _OPENMP
static void nfft_precompute_full_psi_omp(X(plan) *ths)
{
  INT j;                                /**< index over all nodes            */
  INT lprod;                            /**< 'bandwidth' of matrix B         */

  {
    INT t;
    for(t=0,lprod = 1; t<ths->d; t++)
        lprod *= 2*ths->m+2;
  }

  #pragma omp parallel for default(shared) private(j)
  for(j=0; j<ths->M_total; j++)
    {
      INT t,t2;                             /**< index over all dimensions       */
      INT l_L;                              /**< plain index 0<=l_L<lprod        */
      INT lj[ths->d];                       /**< multi index 0<=lj<u+o+1         */
      INT ll_plain[ths->d+1];               /**< postfix plain index             */

      INT u[ths->d], o[ths->d];             /**< depends on x_j                  */

      R phi_prod[ths->d+1];
      INT ix = j*lprod;

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

void X(precompute_full_psi)(X(plan) *ths)
{
#ifdef _OPENMP
  sort(ths);

  nfft_precompute_full_psi_omp(ths);
#else
  INT t, t2; /* index over all dimensions */
  INT j; /* index over all nodes */
  INT l_L; /* plain index 0 <= l_L < lprod */
  INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */
  INT ll_plain[ths->d+1]; /* postfix plain index */
  INT lprod; /* 'bandwidth' of matrix B */
  INT u[ths->d], o[ths->d]; /* depends on x_j */

  R phi_prod[ths->d+1];

  INT ix, ix_old;

  sort(ths);

  phi_prod[0] = K(1.0);
  ll_plain[0] = 0;

  for (t = 0, lprod = 1; t < ths->d; t++)
    lprod *= 2 * ths->m + 2;

  for (j = 0, ix = 0, ix_old = 0; j < ths->M_total; j++)
  {
    MACRO_init_uo_l_lj_t;

    for (l_L = 0; l_L < lprod; l_L++, ix++)
    {
      MACRO_update_phi_prod_ll_plain(without_PRE_PSI);

      ths->psi_index_g[ix] = ll_plain[ths->d];
      ths->psi[ix] = phi_prod[ths->d];

      MACRO_count_uo_l_lj_t;
    } /* for(l_L) */

    ths->psi_index_f[j] = ix - ix_old;
    ix_old = ix;
  } /* for(j) */
#endif
}

void X(precompute_one_psi)(X(plan) *ths)
{
  if(ths->flags & PRE_LIN_PSI)
    X(precompute_lin_psi)(ths);
  if(ths->flags & PRE_FG_PSI)
    X(precompute_fg_psi)(ths);
  if(ths->flags & PRE_PSI)
    X(precompute_psi)(ths);
  if(ths->flags & PRE_FULL_PSI)
    X(precompute_full_psi)(ths);
}

static void init_help(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT lprod; /* 'bandwidth' of matrix B */

  if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    ths->flags |= NFFT_SORT_NODES;

  ths->N_total = intprod(ths->N, 0, ths->d);
  ths->n_total = intprod(ths->n, 0, ths->d);

  ths->sigma = (R*) Y(malloc)((size_t)(ths->d) * sizeof(R));

  for(t = 0;t < ths->d; t++)
    ths->sigma[t] = ((R)ths->n[t]) / (R)(ths->N[t]);

  WINDOW_HELP_INIT;

  if(ths->flags & MALLOC_X)
    ths->x = (R*)Y(malloc)((size_t)(ths->d * ths->M_total) * sizeof(R));

  if(ths->flags & MALLOC_F_HAT)
    ths->f_hat = (C*)Y(malloc)((size_t)(ths->N_total) * sizeof(C));

  if(ths->flags & MALLOC_F)
    ths->f = (C*)Y(malloc)((size_t)(ths->M_total) * sizeof(C));

  if(ths->flags & PRE_PHI_HUT)
    precompute_phi_hut(ths);

  if (ths->flags & PRE_LIN_PSI)
  {
      if (ths->K == 0)
      {
        ths->K = Y(m2K)(ths->m);
      }
      ths->psi = (R*) Y(malloc)((size_t)((ths->K+1) * ths->d) * sizeof(R));
  }

  if(ths->flags & PRE_FG_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * 2) * sizeof(R));

  if(ths->flags & PRE_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * (2 * ths->m + 2)) * sizeof(R));

  if(ths->flags & PRE_FULL_PSI)
  {
      for (t = 0, lprod = 1; t < ths->d; t++)
        lprod *= 2 * ths->m + 2;

      ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(R));

      ths->psi_index_f = (INT*) Y(malloc)((size_t)(ths->M_total) * sizeof(INT));
      ths->psi_index_g = (INT*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(INT));
  }

  if(ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    INT nthreads = Y(get_num_threads)();
#endif

    ths->g1 = (C*)Y(malloc)((size_t)(ths->n_total) * sizeof(C));

    if(ths->flags & FFT_OUT_OF_PLACE)
      ths->g2 = (C*) Y(malloc)((size_t)(ths->n_total) * sizeof(C));
    else
      ths->g2 = ths->g1;

#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    FFTW(plan_with_nthreads)(nthreads);
#endif
    {
      int *_n = Y(malloc)((size_t)(ths->d) * sizeof(int));

      for (t = 0; t < ths->d; t++)
        _n[t] = (int)(ths->n[t]);

      ths->my_fftw_plan1 = FFTW(plan_dft)((int)ths->d, _n, ths->g1, ths->g2, FFTW_FORWARD, ths->fftw_flags);
      ths->my_fftw_plan2 = FFTW(plan_dft)((int)ths->d, _n, ths->g2, ths->g1, FFTW_BACKWARD, ths->fftw_flags);
      Y(free)(_n);
    }
#ifdef _OPENMP
}
#endif
  }

  if(ths->flags & NFFT_SORT_NODES)
    ths->index_x = (INT*) Y(malloc)(sizeof(INT) * 2U * (size_t)(ths->M_total));
  else
    ths->index_x = NULL;

  ths->mv_trafo = (void (*) (void* ))X(trafo);
  ths->mv_adjoint = (void (*) (void* ))X(adjoint);
}

void X(init)(X(plan) *ths, int d, int *N, int M_total)
{
  INT t; /* index over all dimensions */

  ths->d = (INT)d;

  ths->N = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->M_total = (INT)M_total;

  ths->n = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = 2 * (Y(next_power_of_2)(ths->N[t]));

  ths->m = WINDOW_HELP_ESTIMATE_m;

  if (d > 1)
  {
#ifdef _OPENMP
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | NFFT_SORT_NODES |
                 NFFT_OMP_BLOCKWISE_ADJOINT;
#else
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | NFFT_SORT_NODES;
#endif
  }
  else
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE;

  ths->fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;

  ths->K = 0;
  init_help(ths);
}

void X(init_guru)(X(plan) *ths, int d, int *N, int M_total, int *n, int m,
  unsigned flags, unsigned fftw_flags)
{
  INT t; /* index over all dimensions */

  ths->d = (INT)d;
  ths->M_total = (INT)M_total;
  ths->N = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->n = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = (INT)n[t];

  ths->m = (INT)m;

  ths->flags = flags;
  ths->fftw_flags = fftw_flags;

  ths->K = 0;
  init_help(ths);
}

void X(init_lin)(X(plan) *ths, int d, int *N, int M_total, int *n, int m, int K,
  unsigned flags, unsigned fftw_flags)
{
  INT t; /* index over all dimensions */

  ths->d = (INT)d;
  ths->M_total = (INT)M_total;
  ths->N = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->n = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = (INT)n[t];

  ths->m = (INT)m;

  ths->flags = flags;
  ths->fftw_flags = fftw_flags;

  ths->K = K;
  init_help(ths);
}

void X(init_1d)(X(plan) *ths, int N1, int M_total)
{
  int N[1];

  N[0] = N1;

  X(init)(ths, 1, N, M_total);
}

void X(init_2d)(X(plan) *ths, int N1, int N2, int M_total)
{
  int N[2];

  N[0] = N1;
  N[1] = N2;
  X(init)(ths, 2, N, M_total);
}

void X(init_3d)(X(plan) *ths, int N1, int N2, int N3, int M_total)
{
  int N[3];

  N[0] = N1;
  N[1] = N2;
  N[2] = N3;
  X(init)(ths, 3, N, M_total);
}

const char* X(check)(X(plan) *ths)
{
  INT j;

  if (!ths->f)
      return "Member f not initialized.";

  if (!ths->x)
      return "Member x not initialized.";

  if (!ths->f_hat)
      return "Member f_hat not initialized.";

  if ((ths->flags & PRE_LIN_PSI) && ths->K < ths->M_total)
    return "Number of nodes too small to use PRE_LIN_PSI.";

  for (j = 0; j < ths->M_total * ths->d; j++)
  {
    if ((ths->x[j]<-K(0.5)) || (ths->x[j]>= K(0.5)))
    {
      return "ths->x out of range [-0.5,0.5)";
    }
  }

  for (j = 0; j < ths->d; j++)
  {
    if (ths->sigma[j] <= 1)
      return "Oversampling factor too small";
    
    /* Automatically calls trafo_direct if 
    if(ths->N[j] <= ths->m)
      return "Polynomial degree N is <= cut-off m";
    */

    if(ths->N[j]%2 == 1)
      return "polynomial degree N has to be even";
  }
  return 0;
}

void X(finalize)(X(plan) *ths)
{
  INT t; /* index over dimensions */

  if(ths->flags & NFFT_SORT_NODES)
    Y(free)(ths->index_x);

  if(ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_plan2);
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_plan1);

    if(ths->flags & FFT_OUT_OF_PLACE)
      Y(free)(ths->g2);

    Y(free)(ths->g1);
  }

  if(ths->flags & PRE_FULL_PSI)
  {
    Y(free)(ths->psi_index_g);
    Y(free)(ths->psi_index_f);
    Y(free)(ths->psi);
  }

  if(ths->flags & PRE_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_FG_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_LIN_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_PHI_HUT)
  {
    for (t = 0; t < ths->d; t++)
        Y(free)(ths->c_phi_inv[t]);
    Y(free)(ths->c_phi_inv);
  }

  if(ths->flags & MALLOC_F)
    Y(free)(ths->f);

  if(ths->flags & MALLOC_F_HAT)
    Y(free)(ths->f_hat);

  if(ths->flags & MALLOC_X)
    Y(free)(ths->x);

  WINDOW_HELP_FINALIZE;

  Y(free)(ths->sigma);
  Y(free)(ths->n);
  Y(free)(ths->N);
}
