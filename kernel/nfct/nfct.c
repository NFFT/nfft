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

/* Nonequispaced fast cosine transform */

/* Author: Steffen Klatt 2004-2006, Jens Keiner 2010 */

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
#define X(name) NFCT(name)

/** Compute aggregated product of integer array. */
static inline INT intprod(const INT *vec, const INT d)
{
  INT t, p;

  p = 1;
  for (t = 0; t < d; t++)
    p *= vec[t];

  return p;
}

/* handy shortcuts */
#define NFCT_DEFAULT_FLAGS PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | \
  MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE

#define FFTW_DEFAULT_FLAGS FFTW_ESTIMATE | FFTW_DESTROY_INPUT

#define NFCT_SUMMANDS (2 * ths->m + 2)
#define NODE(p,r) (ths->x[(p) * ths->d + (r)])

#define NFCT_WINDOW_HELP_INIT WINDOW_HELP_INIT

#define MACRO_with_PRE_PSI ths->psi[(j * ths->d + t) * NFCT_SUMMANDS + lc[t]]
#define MACRO_compute_PSI PHI(2 * (ths->n[t] - 1), NODE(j,t) - ((R)(lc[t] + lb[t])) / (K(2.0)*((R)(ths->n[t])-K(1.0))), t)

/** direct computation of non equispaced cosine transforms
 *  nfct_trafo_direct,  nfct_adjoint_direct
 *  require O(M N^d) arithemtical operations
 *
 * direct computation of the nfct_trafo_direct, formula (1.1)
 * nfct_trafo_direct:
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * cos(2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint_direct, formula (1.2)
 * nfct_adjoint_direct:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * cos(2 (pi) k x[j])
 */

void X(trafo_direct)(const X(plan) *ths)
{
  R *f_hat = (R*)ths->f_hat, *f = (R*)ths->f;

  memset(f, 0, ths->M_total * sizeof(R));

  if (ths->d == 1)
  {
    /* specialize for univariate case, rationale: faster */
    INT j;
    #pragma omp parallel for default(shared) private(j)
    for (j = 0; j < ths->M_total; j++)
    {
      INT k_L;
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        R omega = K2PI * k_L * ths->x[j];
        f[j] += f_hat[k_L] * COS(omega);
      }
    }
  }
  else
  {
    /* multivariate case */
    INT j;
    #pragma omp parallel for default(shared) private(j)
    for (j = 0; j < ths->M_total; j++)
    {
      R x[ths->d], omega, Omega[ths->d + 1];
      INT t, t2, k_L, k[ths->d];
      Omega[0] = K(1.0);
      for (t = 0; t < ths->d; t++)
      {
        k[t] = 0;
        x[t] = K2PI * ths->x[j * ths->d + t];
        Omega[t+1] = COS(((R)k[t]) * x[t]) * Omega[t];
      }
      omega = Omega[ths->d];

      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        f[j] += f_hat[k_L] * omega;
        {
          for (t = ths->d-1; (t >= 1) && (k[t] == (ths->N[t] - 1)); t--)
            k[t] = 0;

          k[t]++;

          for (t2 = t; t2 < ths->d; t2++)
            Omega[t2+1] = COS(((R)k[t2]) * x[t2]) * Omega[t2];

          omega = Omega[ths->d];
        }
      }
    }
  }
}

void X(adjoint_direct)(const X(plan) *ths)
{
  R *f_hat = (R*)ths->f_hat, *f = (R*)ths->f;

  memset(f_hat, 0, ths->N_total * sizeof(R));

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
          R omega = K2PI * k_L * ths->x[j];
          f_hat[k_L] += f[j] * COS(omega);
        }
      }
#else
      INT j;
      for (j = 0; j < ths->M_total; j++)
      {
        INT k_L;
        for (k_L = 0; k_L < ths->N_total; k_L++)
        {
          R omega = K2PI * k_L * ths->x[j];
          f_hat[k_L] += f[j] * COS(omega);
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
        k[t] = k_temp % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (j = 0; j < ths->M_total; j++)
      {
        R omega = K(1.0);
        for (t = 0; t < ths->d; t++)
          omega *= COS(K2PI * k[t] * ths->x[j * ths->d + t]);
        f_hat[k_L] += f[j] * omega;
      }
    }
#else
    for (j = 0; j < ths->M_total; j++)
    {
      R x[ths->d], omega, Omega[ths->d+1];
      INT t, t2, k[ths->d];
      Omega[0] = K(1.0);
      for (t = 0; t < ths->d; t++)
      {
        k[t] = 0;
        x[t] = K2PI * ths->x[j * ths->d + t];
        Omega[t+1] = COS(((R)k[t]) * x[t]) * Omega[t];
      }
      omega = Omega[ths->d];
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        f_hat[k_L] += f[j] * omega;

        for (t = ths->d-1; (t >= 1) && (k[t] == ths->N[t] - 1); t--)
          k[t] = 0;

        k[t]++;

        for (t2 = t; t2 < ths->d; t2++)
          Omega[t2+1] = COS(((R)k[t2]) * x[t2]) * Omega[t2];

        omega = Omega[ths->d];
      }
    }
#endif
  }
}

/** fast computation of non equispaced cosine transforms
*  require O(N^d log(N) + M) arithemtical operations
*
* fast computation of the nfct_trafo, formula (1.1)
* nfct_trafo:
* for j=0,...,M-1
*  f[j] = sum_{k in I_N^d} f_hat[k] * cos(2 (pi) k x[j])
*
* direct computation of the nfct_adjoint, formula (1.2)
* nfct_adjoint:
* for k in I_N^d
*  f_hat[k] = sum_{j=0}^{M-1} f[j] * cos(2 (pi) k x[j])
*/

#define MACRO_nfct_D_compute_A \
{ \
  g_hat[kg_plain[ths->d]] = f_hat[k_L] * c_phi_inv_k[ths->d]; \
}

#define MACRO_nfct_D_compute_T \
{ \
  f_hat[k_L] = g_hat[kg_plain[ths->d]] * c_phi_inv_k[ths->d]; \
}

#define MACRO_nfct_D_init_result_A memset(g_hat, 0, ths->n_total * sizeof(R));

#define MACRO_nfct_D_init_result_T memset(f_hat, 0, ths->N_total * sizeof(R));

#define MACRO_with_PRE_PHI_HUT ths->c_phi_inv[t][kg[t]]

#define MACRO_compute_PHI_HUT_INV (K(1.0) / (PHI_HUT(2 * (ths->n[t] - 1), kg[t], t)))

#define MACRO_init_k_ks \
{ \
  for (t = 0; t < ths->d; t++) \
  { \
    kg[t] = 0; \
  } \
  i = 0; \
}

#define MACRO_update_c_phi_inv_k(what_kind, which_phi) \
{ \
  for (t = i; t < ths->d; t++) \
  { \
    MACRO_update_c_phi_inv_k_ ## what_kind(which_phi); \
    kg_plain[t+1] = kg_plain[t] * ths->n[t] + kg[t]; \
  } \
}

#define MACRO_update_c_phi_inv_k_A(which_phi) \
{ \
  c_phi_inv_k[t+1] = (kg[t] == 0 ? K(1.0) : K(0.5)) * c_phi_inv_k[t] * MACRO_ ## which_phi; \
}

#define MACRO_update_c_phi_inv_k_T(which_phi) \
{ \
  c_phi_inv_k[t+1] = c_phi_inv_k[t] * MACRO_ ## which_phi; \
}

#define MACRO_count_k_ks \
{ \
  kg[ths->d - 1]++; \
  i = ths->d - 1; \
\
  while ((kg[i] == ths->N[i]) && (i > 0)) \
  { \
    kg[i - 1]++; \
    kg[i] = 0; \
    i--; \
  } \
}

/* sub routines for the fast transforms  matrix vector multiplication with D, D^T */
#define MACRO_nfct_D(which_one) \
static inline void D_ ## which_one (X(plan) *ths) \
{ \
  R *g_hat, *f_hat; /* local copy */ \
  R c_phi_inv_k[ths->d+1]; /* postfix product of PHI_HUT */ \
  INT i, t; \
  INT k_L; /* plain index */ \
  INT kg[ths->d]; /* multi index in g_hat */ \
  INT kg_plain[ths->d+1]; /* postfix plain index */ \
\
  g_hat = ths->g_hat; f_hat = ths->f_hat; \
  MACRO_nfct_D_init_result_ ## which_one; \
\
  c_phi_inv_k[0] = K(1.0); \
  kg_plain[0] = 0; \
\
  MACRO_init_k_ks; \
\
  if (ths->flags & PRE_PHI_HUT) \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(which_one, with_PRE_PHI_HUT); \
      MACRO_nfct_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
  else \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(which_one,compute_PHI_HUT_INV); \
      MACRO_nfct_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
}

MACRO_nfct_D(A)

MACRO_nfct_D(T)

/* sub routines for the fast transforms matrix vector multiplication with B, B^T */
#define MACRO_nfct_B_init_result_A memset(f, 0, ths->M_total * sizeof(R));
#define MACRO_nfct_B_init_result_T memset(g, 0, ths->n_total * sizeof(R));

#define MACRO_nfct_B_PRE_FULL_PSI_compute_A \
{ \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]]; \
}

#define MACRO_nfct_B_PRE_FULL_PSI_compute_T \
{ \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj); \
}

#define MACRO_nfct_B_compute_A \
{ \
  (*fj) += phi_tilde[ths->d] * g[lg_plain[ths->d]]; \
}

#define MACRO_nfct_B_compute_T \
{ \
  g[lg_plain[ths->d]] += phi_tilde[ths->d] * (*fj); \
}

#define MACRO_compute_lg_offset__count_lg(i0) \
{ \
  /* determine index in g-array corresponding to lb[(i0)] */ \
  if (lb[(i0)] < 0) \
    lg_offset[(i0)] = \
      (lb[(i0)] % (2 * (ths->n[(i0)] - 1))) + (2 * (ths->n[(i0)] - 1)); \
  else \
    lg_offset[(i0)] = lb[(i0)] % (2 * (ths->n[(i0)] - 1)); \
    if (lg_offset[(i0)] >= ths->n[(i0)]) \
      lg_offset[(i0)] = -(2 * (ths->n[(i0)] - 1) - lg_offset[(i0)]); \
}

#define MACRO_set__lg__to__lg_offset \
{ \
  if (lg_offset[i] <= 0) \
  { \
    lg[i] = -lg_offset[i]; \
    count_lg[i] = -1; \
  } \
  else \
  { \
    lg[i] = +lg_offset[i]; \
    count_lg[i] = +1; \
  } \
}

#define MACRO_count__lg(dim) \
{ \
  /* turn around if we hit one of the boundaries */ \
  if ((lg[(dim)] == 0) || (lg[(dim)] == ths->n[(dim)]-1)) \
    count_lg[(dim)] *= -1; \
  /* move array index */ \
  lg[(dim)] += count_lg[(dim)]; \
}

#define MACRO_nfct__lower_boundary(j,act_dim) \
{ \
  lb[(act_dim)] = \
    LRINT(NODE((j),(act_dim)) * (2 * (ths->n[(act_dim)] - 1))) - ths->m; \
}

#define MACRO_init_lb_lg_lc \
{ \
  for (i = 0; i < ths->d; i++) \
  { \
    MACRO_nfct__lower_boundary(j, i); \
    MACRO_compute_lg_offset__count_lg(i); \
    MACRO_set__lg__to__lg_offset; \
    /* counter for lg */ \
    lc[i] = 0; \
   } \
   i = 0; \
}

#define MACRO_count__lg_lc \
{ \
  MACRO_count__lg(ths->d-1); \
  lc[ths->d - 1]++; \
  i = ths->d - 1; \
  while ((lc[i] == NFCT_SUMMANDS) && (i > 0)) \
  { \
    lc[i - 1]++; \
    lc[i] = 0; \
    /* ansonsten lg[i-1] verschieben */ \
    MACRO_count__lg(i - 1); \
    /* lg[i] = anfangswert */ \
    MACRO_set__lg__to__lg_offset; \
    i--; \
  } \
}

#define  MACRO_update_phi_tilde_lg_plain(which_one, which_psi) \
{ \
  for (t = i; t < ths->d; t++) \
  { \
    MACRO__phi_tilde__ ## which_one(which_psi); \
    lg_plain[t+1]  = lg_plain[t]  * ths->n[t] + lg[t]; \
  } \
}

#define MACRO__phi_tilde__A(which_psi) \
{ \
  phi_tilde[t+1] = phi_tilde[t] * MACRO_ ## which_psi; \
}

#define MACRO__phi_tilde__T(which_psi) \
{ \
  if(lg[t] == 0 || lg[t] == ths->n[t] - 1) \
  { \
    phi_tilde[t+1] = phi_tilde[t] * MACRO_ ## which_psi; \
  } \
  else \
  { \
    phi_tilde[t+1] = K(0.5) * phi_tilde[t] * MACRO_ ## which_psi; \
  } \
}

#define MACRO_nfct_B(which_one) \
static inline void B_ ## which_one (nfct_plan *ths) \
{ /* MACRO_nfct_B */ \
  INT lb[ths->d]; /* multi band with respect to x_j */ \
  INT j, t, i; /* index nodes, help vars */ \
  INT lprod, l_L, ix; /* index one row of B */ \
  INT lc[ths->d]; /* multi index 0<=lc<2m+2 */ \
  INT lg[ths->d]; /* real index of g in array */ \
  INT lg_offset[ths->d]; /* offset in g according to u */ \
  INT count_lg[ths->d]; /* count summands (2m+2) */ \
  INT lg_plain[ths->d+1]; /* index of g in multi_array */ \
  R *f, *g; /* local copy */ \
  R phi_tilde[ths->d+1]; /* holds values for psi */ \
  R *fj; /* pointer to final result */ \
\
  f = ths->f; g = ths->g; \
\
  MACRO_nfct_B_init_result_ ## which_one \
\
  /* both flags are set */ \
  if ((ths->flags & PRE_PSI)&&(ths->flags & PRE_FULL_PSI)) \
  { \
    for (ix = 0, j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1) \
      for (l_L = 0; l_L < ths->psi_index_f[j]; l_L++, ix++) \
      { \
        MACRO_nfct_B_PRE_FULL_PSI_compute_ ## which_one; \
      } \
  } \
  else \
  { \
    phi_tilde[0] = K(1.0); \
    lg_plain[0]  = 0; \
\
    for (t = 0, lprod = 1; t < ths->d; t++) \
      lprod *= NFCT_SUMMANDS; \
\
    /* PRE_PSI flag is set */ \
    if (ths->flags & PRE_PSI) \
      for (j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1) \
        { \
          MACRO_init_lb_lg_lc; \
          for (l_L = 0; l_L < lprod; l_L++) \
          { \
            MACRO_update_phi_tilde_lg_plain(which_one, with_PRE_PSI); \
            MACRO_nfct_B_compute_ ## which_one; \
            MACRO_count__lg_lc; \
          } /* for(l_L) */ \
        } /* for(j) */ \
\
    /* no PSI flag is set */ \
    else \
      for (j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1) \
      { \
        MACRO_init_lb_lg_lc; \
        for (l_L = 0; l_L < lprod; l_L++) \
        { \
          MACRO_update_phi_tilde_lg_plain(which_one,compute_PSI); \
          MACRO_nfct_B_compute_ ## which_one; \
          MACRO_count__lg_lc; \
        } /* for (l_L) */ \
      } /* for (j) */ \
  } /* else (PRE_PSI && FULL_PRE_PSI) */ \
} /* nfct_B */

MACRO_nfct_B(A)
MACRO_nfct_B(T)

/* more memory, but faster */
#define MACRO_nfct_full_psi(which_one) \
static inline void full_psi__ ## which_one(nfct_plan *ths) \
{ \
  INT t, i; /* index over all dimensions */ \
  INT j; /* node index */ \
  INT l_L; /* plain index 0 <= l_L < lprod */ \
  INT lc[ths->d]; /* multi index 0<=lj<u+o+1 */ \
  INT lg_plain[ths->d+1]; /* postfix plain index */ \
  INT count_lg[ths->d]; \
  INT lg_offset[ths->d]; \
  INT lg[ths->d]; \
  INT lprod; /* 'bandwidth' of matrix B */ \
  INT lb[ths->d]; /* depends on x_j */ \
\
  R phi_tilde[ths->d+1]; \
  R eps = ths->nfct_full_psi_eps; \
\
  INT *index_g, *index_f; \
  R *new_psi; \
  INT ix, ix_old, size_psi; \
\
  phi_tilde[0] = K(1.0); \
  lg_plain[0]  =   0; \
 \
  if (ths->flags & PRE_PSI) \
  { \
    size_psi = ths->M_total; \
    index_f = (INT*)Y(malloc)(ths->M_total  * sizeof(INT)); \
    index_g = (INT*)Y(malloc)(size_psi * sizeof(INT)); \
    new_psi = (R*)Y(malloc)(size_psi * sizeof(R)); \
\
    for (t = 0,lprod = 1; t < ths->d; t++) \
    { \
      lprod *= NFCT_SUMMANDS; \
      eps *= PHI(2 * (ths->n[t] - 1), K(0.0), t); \
    } \
\
    for (ix = 0, ix_old = 0, j = 0; j < ths->M_total; j++) \
    { \
      MACRO_init_lb_lg_lc; \
\
      for (l_L = 0; l_L < lprod; l_L++) \
      { \
        MACRO_update_phi_tilde_lg_plain(which_one, with_PRE_PSI); \
\
        if (phi_tilde[ths->d] > eps) \
        { \
          index_g[ix] = lg_plain[ths->d]; \
          new_psi[ix] = phi_tilde[ths->d]; \
\
          ix++; \
          if (ix == size_psi) \
          { \
            size_psi += ths->M_total; \
            index_g = (INT*)realloc(index_g, size_psi * sizeof(INT)); \
            new_psi = (R*)realloc(new_psi, size_psi * sizeof(R)); \
          } \
        } \
\
        MACRO_count__lg_lc; \
\
      } /* for (l_L) */ \
\
      index_f[j] = ix - ix_old; \
      ix_old = ix; \
\
    } /* for(j) */ \
\
    Y(free)(ths->psi); \
    size_psi = ix; \
    ths->size_psi = size_psi; \
\
    index_g = (INT*)realloc(index_g, size_psi * sizeof(INT)); \
    new_psi = (R*)realloc(new_psi, size_psi * sizeof(R)); \
\
    ths->psi = new_psi; \
    ths->psi_index_g = index_g; \
    ths->psi_index_f = index_f; \
\
  } /* if(PRE_PSI) */ \
}

MACRO_nfct_full_psi(A)
MACRO_nfct_full_psi(T)

/* user routines */

void X(trafo)(X(plan) *ths)
{
  /* use ths->my_fftw_r2r_plan */
  ths->g_hat = ths->g1;
  ths->g = ths->g2;

  /* form \f$ \hat g_k = \frac{\hat f_k}{c_k\left(\phi\right)} \text{ for }
   * k \in I_N \f$ */
  TIC(0)
  D_A(ths);
  TOC(0)

  /* Compute by d-variate discrete Fourier transform
   * \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   * \text{ for } l \in I_n \f$ */
  TIC(1)
  Z(execute)(ths->my_fftw_r2r_plan);
  TOC(1)

  if (ths->flags & PRE_FULL_PSI)
    full_psi__A(ths);

  /* Set \f$ f_j = \sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
   * \text{ for } j=0,\hdots,M-1 \f$ */
  TIC(2)
  B_A(ths);
  TOC(2)

  if (ths->flags & PRE_FULL_PSI)
  {
    Y(free)(ths->psi_index_g);
    Y(free)(ths->psi_index_f);
  }
} /* nfct_trafo */

void X(adjoint)(X(plan) *ths)
{
  /* use ths->my_fftw_plan */
  ths->g_hat = ths->g2;
  ths->g = ths->g1;

  if (ths->flags & PRE_FULL_PSI)
    full_psi__T(ths);

  /* Set \f$ g_l = \sum_{j=0}^{M-1} f_j \psi\left(x_j-\frac{l}{n}\right)
   * \text{ for } l \in I_n,m(x_j) \f$ */
  TIC(2)
  B_T(ths);
  TOC(2)

  if (ths->flags & PRE_FULL_PSI)
  {
    Y(free)(ths->psi_index_g);
    Y(free)(ths->psi_index_f);
  }

  /* Compute by d-variate discrete cosine transform
   * \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   * \text{ for }  k \in I_N\f$ */
  TIC(1)
  Z(execute)(ths->my_fftw_r2r_plan);
  TOC(1)

  /* Form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
   * k \in I_N \f$ */
  TIC(0)
  D_T(ths);
  TOC(0)

} /* nfct_adjoint */

static inline void precompute_phi_hut(X(plan) *ths)
{
  INT kg[ths->d]; /* index over all frequencies */
  INT t; /* index over all dimensions */

  ths->c_phi_inv = (R**)Y(malloc)(ths->d * sizeof(R*));

  for (t = 0; t < ths->d; t++)
  {
    ths->c_phi_inv[t] = (R*)Y(malloc)(ths->N[t] * sizeof(R));

    for (kg[t] = 0; kg[t] < ths->N[t]; kg[t]++)
    {
      ths->c_phi_inv[t][kg[t]] = MACRO_compute_PHI_HUT_INV;
    }
  }
} /* nfct_phi_hut */

void X(precompute_psi)(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT j; /* index over all nodes */
  INT lc[ths->d]; /* index 0<=lj<u+o+1 */
  INT lb[ths->d]; /* depends on x_j */

  for (t = 0; t < ths->d; t++)
  {
    for (j = 0; j < ths->M_total; j++)
    {
      MACRO_nfct__lower_boundary(j, t);
      for(lc[t] = 0; lc[t] < NFCT_SUMMANDS; lc[t]++)
	      ths->psi[(j * ths->d + t) * NFCT_SUMMANDS + lc[t]] = MACRO_compute_PSI;
    } /* for (j) */
  } /* for (t) */
} /* nfct_precompute_psi */

static inline void init_help(X(plan) *ths)
{
  INT t; /* index over all dimensions */

  ths->N_total = intprod(ths->N, ths->d);
  ths->n_total = intprod(ths->n, ths->d);
  ths->sigma = (R*)Y(malloc)(ths->d * sizeof(R));

  for (t = 0; t < ths->d; t++)
    ths->sigma[t] = ((R)(ths->n[t] - 1)) / ths->N[t];

  /* Assign r2r transform kinds for each dimension */
  ths->r2r_kind = (Z(r2r_kind)*)Y(malloc)(ths->d * sizeof (Z(r2r_kind)));
  for (t = 0; t < ths->d; t++)
    ths->r2r_kind[t] = FFTW_REDFT00;

  NFCT_WINDOW_HELP_INIT;

  if (ths->flags & MALLOC_X)
    ths->x = (R*)Y(malloc)(ths->d * ths->M_total * sizeof(R));

  if (ths->flags & MALLOC_F_HAT)
    ths->f_hat = (R*)Y(malloc)(ths->N_total * sizeof(R));

  if (ths->flags & MALLOC_F)
    ths->f = (R*)Y(malloc)(ths->M_total * sizeof(R));

  if (ths->flags & PRE_PHI_HUT)
    precompute_phi_hut(ths);

  /* NO FFTW_MALLOC HERE */
  if (ths->flags & PRE_PSI)
  {
    ths->psi =
      (R*)Y(malloc)(ths->M_total * ths->d * NFCT_SUMMANDS * sizeof(R));

    /* Set default for full_psi_eps */
    ths->nfct_full_psi_eps = POW(K(10.0), K(-10.0));
  }

  if (ths->flags & FFTW_INIT)
  {
    ths->g1 =
      (R*)Y(malloc)(intprod(ths->n, ths->d) * sizeof(R));

    if (ths->flags & FFT_OUT_OF_PLACE)
      ths->g2 =
	      (R*) Y(malloc)(intprod(ths->n, ths->d) * sizeof(R));
    else
      ths->g2 = ths->g1;

    {
      int *n = Y(malloc)(ths->d * sizeof(int));
      for (t = 0; t < ths->d; t++)
        n[t] = (int)(ths->n[t]);
      ths->my_fftw_r2r_plan =
        Z(plan_r2r)(ths->d, n, ths->g1, ths->g2, ths->r2r_kind,
          ths->fftw_flags);
    }
  }

  ths->mv_trafo = (void (*) (void* ))X(trafo);
  ths->mv_adjoint = (void (*) (void* ))X(adjoint);
}

void X(init)(X(plan) *ths, int d, int *N, int M_total)
{
  int t;

  ths->d = d;
  ths->M_total = M_total;
  ths->N = (INT*) Y(malloc)(ths->d * sizeof(INT));

  for (t = 0;t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->n = (INT*) Y(malloc)(ths->d * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = 2 * (Y(next_power_of_2)(ths->N[t]) - 1);

  ths->m = WINDOW_HELP_ESTIMATE_m;

  ths->flags = NFCT_DEFAULT_FLAGS;
  ths->fftw_flags = FFTW_DEFAULT_FLAGS;

  init_help(ths);
}

/* Was macht diese Funktion. Wird sie gebraucht? Bei NFST ist sie auch in
 * nfft3.h deklariert.
void nfct_init_m(nfct_plan *ths, int d, int *N, int M_total, int m)
{
  int t, n[d];

  for(t = 0; t < d; t++)
    n[t] = 2 * (Y(next_power_of_2)(N[t]) - 1);

  nfct_init_guru(ths, d, N, M_total, n, m, NFCT_DEFAULT_FLAGS, FFTW_DEFAULT_FLAGS);
}
*/

void X(init_guru)(X(plan) *ths, int d, int *N, int M_total, int *n, int m,
  unsigned flags, unsigned fftw_flags)
{
  INT t; /* index over all dimensions */

  ths->d = d;
  ths->M_total = M_total;

  ths->N = (INT*)Y(malloc)(ths->d * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->n = (INT*)Y(malloc)(ths->d * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = (INT)n[t];

  ths->m = m;

  ths->flags = flags;
  ths->fftw_flags = fftw_flags;

  init_help(ths);
}

void X(init_1d)(X(plan) *ths, int N0, int M_total)
{
  int N[1];

  N[0] = N0;
  X(init)(ths, 1, N, M_total);
}

void X(init_2d)(X(plan) *ths, int N0, int N1, int M_total)
{
  int N[2];

  N[0] = N0;
  N[1] = N1;
  X(init)(ths, 2, N, M_total);
}

void X(init_3d)(X(plan) *ths, int N0, int N1, int N2, int M_total)
{
  int N[3];

  N[0] = N0;
  N[1] = N1;
  N[2] = N2;
  X(init)(ths, 3, N, M_total);
}

const char* X(check)(nfct_plan *ths)
{
  INT j;

  for(j=0;j<ths->M_total*ths->d;j++)
    if((ths->x[j] < K(0.0)) || (ths->x[j] >= K(0.5)))
    {
      fprintf(stderr, "\nj = %d, x[j] = " __FE__ "\n", j, ths->x[j]);
      return "ths->x out of range [0.0,0.5)";
    }

  for(j=0;j<ths->d;j++)
  {
    if(ths->sigma[j]<=1)
      return "nfft_check: oversampling factor too small";

    if(ths->N[j] - 1 <= ths->m)
      return "Polynomial degree N is smaller than cut-off m";

    if(ths->N[j]%2==1)
      return "polynomial degree N has to be even";
  }
  return 0;
}

void X(finalize)(X(plan) *ths)
{
  INT t; /* dimension index*/

  if (ths->flags & FFTW_INIT)
  {
    Z(destroy_plan)(ths->my_fftw_r2r_plan);

    if (ths->flags & FFT_OUT_OF_PLACE)
      Y(free)(ths->g2);

    Y(free)(ths->g1);
  }

  /* NO FFTW_FREE HERE */
  if (ths->flags & PRE_PSI)
  {
    Y(free)(ths->psi);
  }

  if (ths->flags & PRE_PHI_HUT)
  {
    for (t = 0; t < ths->d; t++)
      Y(free)(ths->c_phi_inv[t]);
    Y(free)(ths->c_phi_inv);
  }

  if (ths->flags & MALLOC_F)
    Y(free)(ths->f);

  if(ths->flags & MALLOC_F_HAT)
    Y(free)(ths->f_hat);

  if (ths->flags & MALLOC_X)
    Y(free)(ths->x);

  WINDOW_HELP_FINALIZE;

  Y(free)(ths->N);
  Y(free)(ths->n);
  Y(free)(ths->sigma);

  Y(free)(ths->r2r_kind);
} /* nfct_finalize */
