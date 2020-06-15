/*
 * Copyright (c) 2002, 2020 Jens Keiner, Stefan Kunis, Daniel Potts
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
#define X(name) NFST(name)

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
#define BASE(x) SIN(x)
#define NN(x) (x + 1)
#define OFFSET 1
#define FOURIER_TRAFO FFTW_RODFT00
#define FFTW_DEFAULT_FLAGS FFTW_ESTIMATE | FFTW_DESTROY_INPUT

#define NODE(p,r) (ths->x[(p) * ths->d + (r)])

#define MACRO_with_FG_PSI fg_psi[t][lj[t]]
#define MACRO_with_PRE_PSI ths->psi[(j * ths->d + t) * (2 * ths->m + 2) + lj[t]]
#define MACRO_without_PRE_PSI PHI((2 * NN(ths->n[t])), ((ths->x[(j) * ths->d + t]) \
  - ((R)(lj[t] + u[t])) / (K(2.0) * ((R)NN(ths->n[t])))), t)
#define MACRO_compute_PSI PHI((2 * NN(ths->n[t])), (NODE(j,t) - ((R)(lj[t] + u[t])) / (K(2.0) * ((R)NN(ths->n[t])))), t)

/**
 * Direct computation of non equispaced sine transforms
 *  nfst_trafo_direct,  nfst_adjoint_direct
 *  require O(M N^d) arithemtical operations
 *
 * direct computation of the nfst_trafo_direct, formula (1.1)
 * nfst_trafo_direct:
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * sin(2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint_direct, formula (1.2)
 * nfst_adjoint_direct:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * sin(2 (pi) k x[j])
 */
void X(trafo_direct)(const X(plan) *ths)
{
  R *f_hat = (R*)ths->f_hat, *f = (R*)ths->f;

  memset(f, 0, (size_t)(ths->M_total) * sizeof(R));

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
        R omega = K2PI * ((R)(k_L + OFFSET)) * ths->x[j];
        f[j] += f_hat[k_L] * BASE(omega);
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
      Omega[0] = K(1.0);
      for (t = 0; t < ths->d; t++)
      {
        k[t] = OFFSET;
        x[t] = K2PI * ths->x[j * ths->d + t];
        Omega[t+1] = BASE(((R)(k[t])) * x[t]) * Omega[t];
      }
      omega = Omega[ths->d];

      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        f[j] += f_hat[k_L] * omega;
        {
          for (t = ths->d - 1; (t >= 1) && (k[t] == (ths->N[t] - 1)); t--)
            k[t] = OFFSET;

          k[t]++;

          for (t2 = t; t2 < ths->d; t2++)
            Omega[t2+1] = BASE(((R)(k[t2])) * x[t2]) * Omega[t2];

          omega = Omega[ths->d];
        }
      }
    }
  }
}

void X(adjoint_direct)(const X(plan) *ths)
{
  R *f_hat = (R*)ths->f_hat, *f = (R*)ths->f;

  memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(R));

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
          R omega = K2PI * ((R)(k_L + OFFSET)) * ths->x[j];
          f_hat[k_L] += f[j] * BASE(omega);
        }
      }
#else
      INT j;
      for (j = 0; j < ths->M_total; j++)
      {
        INT k_L;
        for (k_L = 0; k_L < ths->N_total; k_L++)
        {
          R omega = K2PI * ((R)(k_L + OFFSET)) * ths->x[j];
          f_hat[k_L] += f[j] * BASE(omega);
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
          omega *= BASE(K2PI * (k[t] + OFFSET) * ths->x[j * ths->d + t]);
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
        k[t] = OFFSET;
        x[t] = K2PI * ths->x[j * ths->d + t];
        Omega[t+1] = BASE(((R)(k[t])) * x[t]) * Omega[t];
      }
      omega = Omega[ths->d];
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        f_hat[k_L] += f[j] * omega;

        for (t = ths->d-1; (t >= 1) && (k[t] == ths->N[t] - 1); t--)
          k[t] = OFFSET;

        k[t]++;

        for (t2 = t; t2 < ths->d; t2++)
          Omega[t2+1] = BASE(((R)(k[t2])) * x[t2]) * Omega[t2];

        omega = Omega[ths->d];
      }
    }
#endif
  }
}

/** fast computation of non equispaced sine transforms
 *  require O(N^d log(N) + M) arithemtical operations
 *
 * fast computation of the nfst_trafo, formula (1.1)
 * nfst_trafo:
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * sin(2 (pi) k x[j])
 *
 * direct computation of the nfst_adjoint, formula (1.2)
 * nfst_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * sin(2 (pi) k x[j])
 */

/** macros and small sub routines for the fast transforms
 */

/** computes 2m+2 indices for the matrix B
 */
static inline void uo(const X(plan) *ths, const INT j, INT *up, INT *op,
  const INT act_dim)
{
  const R xj = ths->x[j * ths->d + act_dim];
  INT c = LRINT(xj * (2 * NN(ths->n[(act_dim)])));

  (*up) = c - (ths->m);
  (*op) = c + 1 + (ths->m);
}

#define MACRO_D_compute_A \
{ \
  g_hat[kg_plain[ths->d]] = f_hat[k_L] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_compute_T \
{ \
  f_hat[k_L] = g_hat[kg_plain[ths->d]] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_init_result_A memset(g_hat, 0, (size_t)(ths->n_total) * sizeof(R));

#define MACRO_D_init_result_T memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(R));

#define MACRO_with_PRE_PHI_HUT ths->c_phi_inv[t][kg[t]]

#define MACRO_compute_PHI_HUT_INV (K(1.0) / (PHI_HUT((2 * NN(ths->n[t])), kg[t] + OFFSET, t)))

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
  c_phi_inv_k[t+1] = K(0.5) * c_phi_inv_k[t] * MACRO_ ## which_phi; \
}

#define MACRO_update_c_phi_inv_k_T(which_phi) \
{ \
  c_phi_inv_k[t+1] = K(0.5) * c_phi_inv_k[t] * MACRO_ ## which_phi; \
}

#define MACRO_count_k_ks \
{ \
  kg[ths->d - 1]++; \
  i = ths->d - 1; \
\
  while ((kg[i] == ths->N[i] - 1) && (i > 0)) \
  { \
    kg[i - 1]++; \
    kg[i] = 0; \
    i--; \
  } \
}

/* sub routines for the fast transforms  matrix vector multiplication with D, D^T */
#define MACRO_D(which_one) \
static inline void D_ ## which_one (X(plan) *ths) \
{ \
  R *g_hat, *f_hat; /* local copy */ \
  R c_phi_inv_k[ths->d+1]; /* postfix product of PHI_HUT */ \
  INT t; /* index dimensions */ \
  INT i; \
  INT k_L; /* plain index */ \
  INT kg[ths->d]; /* multi index in g_hat */ \
  INT kg_plain[ths->d+1]; /* postfix plain index */ \
\
  f_hat = (R*)ths->f_hat; g_hat = (R*)ths->g_hat; \
  MACRO_D_init_result_ ## which_one; \
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
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
  else \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(which_one,compute_PHI_HUT_INV); \
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
}

MACRO_D(A)
MACRO_D(T)

/* sub routines for the fast transforms matrix vector multiplication with B, B^T */
#define MACRO_B_init_result_A memset(f, 0, (size_t)(ths->M_total) * sizeof(R));
#define MACRO_B_init_result_T memset(g, 0, (size_t)(ths->n_total) * sizeof(R));

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
  (*fj) += phi_prod[ths->d] * g[ll_plain[ths->d]]; \
}

#define MACRO_B_compute_T \
{ \
  g[ll_plain[ths->d]] += phi_prod[ths->d] * (*fj); \
}

#define MACRO_init_uo_l_lj_t \
{ \
  for (t2 = 0; t2 < ths->d; t2++) \
  { \
    uo(ths, j, &u[t2], &o[t2], t2); \
    \
    /* determine index in g-array corresponding to u[(t2)] */ \
    if (u[(t2)] < 0) \
      lg_offset[(t2)] = \
        (u[(t2)] % (2 * NN(ths->n[(t2)]))) + (2 * NN(ths->n[(t2)])); \
    else \
      lg_offset[(t2)] = u[(t2)] % (2 * NN(ths->n[(t2)])); \
      if (lg_offset[(t2)] > NN(ths->n[(t2)])) \
        lg_offset[(t2)] = -(2 * NN(ths->n[(t2)]) - lg_offset[(t2)]); \
    \
    if (lg_offset[t2] <= 0) \
    { \
      l[t2] = -lg_offset[t2]; \
      count_lg[t2] = -1; \
    } \
    else \
    { \
      l[t2] = +lg_offset[t2]; \
      count_lg[t2] = +1; \
    } \
 \
    lj[t2] = 0; \
   } \
   t2 = 0; \
}

#define FOO_A ((R)count_lg[t])

#define FOO_T ((R)count_lg[t])

#define MACRO_update_phi_prod_ll_plain(which_one,which_psi) \
{ \
  for (t = t2; t < ths->d; t++) \
  { \
    if ((l[t] != 0) && (l[t] != NN(ths->n[t]))) \
    { \
      phi_prod[t+1] = (FOO_ ## which_one) * phi_prod[t] * (MACRO_ ## which_psi); \
      ll_plain[t+1]  = ll_plain[t] * ths->n[t] + l[t] - 1; \
    } \
    else \
    { \
      phi_prod[t + 1] = K(0.0); \
      ll_plain[t+1]  = ll_plain[t] * ths->n[t]; \
    } \
  } \
}

#define MACRO_count_uo_l_lj_t \
{ \
  /* turn around if we hit one of the boundaries */ \
  if ((l[(ths->d-1)] == 0) || (l[(ths->d-1)] == NN(ths->n[(ths->d-1)]))) \
    count_lg[(ths->d-1)] *= -1; \
 \
  /* move array index */ \
  l[(ths->d-1)] += count_lg[(ths->d-1)]; \
 \
  lj[ths->d - 1]++; \
  t2 = ths->d - 1; \
 \
  while ((lj[t2] == (2 * ths->m + 2)) && (t2 > 0)) \
  { \
    lj[t2 - 1]++; \
    lj[t2] = 0; \
    /* ansonsten lg[i-1] verschieben */ \
 \
    /* turn around if we hit one of the boundaries */ \
    if ((l[(t2 - 1)] == 0) || (l[(t2 - 1)] == NN(ths->n[(t2 - 1)]))) \
      count_lg[(t2 - 1)] *= -1; \
    /* move array index */ \
    l[(t2 - 1)] += count_lg[(t2 - 1)]; \
 \
    /* lg[i] = anfangswert */ \
    if (lg_offset[t2] <= 0) \
    { \
      l[t2] = -lg_offset[t2]; \
      count_lg[t2] = -1; \
    } \
    else \
    { \
      l[t2] = +lg_offset[t2]; \
      count_lg[t2] = +1; \
    } \
 \
    t2--; \
  } \
}

#define MACRO_B(which_one) \
static inline void B_ ## which_one (X(plan) *ths) \
{ \
  INT lprod; /* 'regular bandwidth' of matrix B  */ \
  INT u[ths->d], o[ths->d]; /* multi band with respect to x_j */ \
  INT t, t2; /* index dimensions */ \
  INT j; /* index nodes */ \
  INT l_L, ix; /* index one row of B */ \
  INT l[ths->d]; /* multi index u<=l<=o (real index of g in array) */ \
  INT lj[ths->d]; /* multi index 0<=lc<2m+2 */ \
  INT ll_plain[ths->d+1]; /* postfix plain index in g */ \
  R phi_prod[ths->d+1]; /* postfix product of PHI */ \
  R *f, *g; /* local copy */ \
  R *fj; /* local copy */ \
  R y[ths->d]; \
  R fg_psi[ths->d][2*ths->m+2]; \
  R fg_exp_l[ths->d][2*ths->m+2]; \
  INT l_fg,lj_fg; \
  R tmpEXP1, tmpEXP2, tmpEXP2sq, tmp1, tmp2, tmp3; \
  R ip_w; \
  INT ip_u; \
  INT ip_s = ths->K/(ths->m+2); \
  INT lg_offset[ths->d]; /* offset in g according to u */ \
  INT count_lg[ths->d]; /* count summands (2m+2) */ \
\
  f = (R*)ths->f; g = (R*)ths->g; \
\
  MACRO_B_init_result_ ## which_one \
\
  if (ths->flags & PRE_FULL_PSI) \
  { \
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
    for (j = 0, fj = f; j < ths->M_total; j++, fj++) \
    { \
      MACRO_init_uo_l_lj_t; \
 \
      for (l_L = 0; l_L < lprod; l_L++) \
      { \
        MACRO_update_phi_prod_ll_plain(which_one, with_PRE_PSI); \
 \
        MACRO_B_compute_ ## which_one; \
 \
        MACRO_count_uo_l_lj_t; \
      } /* for(l_L) */ \
    } /* for(j) */ \
    return; \
  } /* if(PRE_PSI) */ \
 \
  if (ths->flags & PRE_FG_PSI) \
  { \
    for (t = 0; t < ths->d; t++) \
    { \
      tmpEXP2 = EXP(K(-1.0) / ths->b[t]); \
      tmpEXP2sq = tmpEXP2 * tmpEXP2; \
      tmp2 = K(1.0); \
      tmp3 = K(1.0); \
      fg_exp_l[t][0] = K(1.0); \
 \
      for (lj_fg = 1; lj_fg <= (2 * ths->m + 2); lj_fg++) \
      { \
        tmp3 = tmp2 * tmpEXP2; \
        tmp2 *= tmpEXP2sq; \
        fg_exp_l[t][lj_fg] = fg_exp_l[t][lj_fg-1] * tmp3; \
      } \
    } \
 \
    for (j = 0, fj = f; j < ths->M_total; j++, fj++) \
    { \
      MACRO_init_uo_l_lj_t; \
 \
      for (t = 0; t < ths->d; t++) \
      { \
        fg_psi[t][0] = ths->psi[2 * (j * ths->d + t)]; \
        tmpEXP1 = ths->psi[2 * (j * ths->d + t) + 1]; \
        tmp1 = K(1.0); \
 \
        for (l_fg = u[t] + 1, lj_fg = 1; l_fg <= o[t]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t][lj_fg] = fg_psi[t][0] * tmp1 * fg_exp_l[t][lj_fg]; \
        } \
      } \
 \
      for (l_L= 0; l_L < lprod; l_L++) \
      { \
        MACRO_update_phi_prod_ll_plain(which_one, with_FG_PSI); \
 \
        MACRO_B_compute_ ## which_one; \
 \
        MACRO_count_uo_l_lj_t; \
      } \
    } \
    return; \
  } \
 \
  if (ths->flags & FG_PSI) \
  { \
    for (t = 0; t < ths->d; t++) \
    { \
      tmpEXP2 = EXP(K(-1.0) / ths->b[t]); \
      tmpEXP2sq = tmpEXP2 * tmpEXP2; \
      tmp2 = K(1.0); \
      tmp3 = K(1.0); \
      fg_exp_l[t][0] = K(1.0); \
      for (lj_fg = 1; lj_fg <= (2 * ths->m + 2); lj_fg++) \
      { \
        tmp3 = tmp2 * tmpEXP2; \
        tmp2 *= tmpEXP2sq; \
        fg_exp_l[t][lj_fg] = fg_exp_l[t][lj_fg-1] * tmp3; \
      } \
    } \
 \
    for (j = 0, fj = f; j < ths->M_total; j++, fj++) \
    { \
      MACRO_init_uo_l_lj_t; \
 \
      for (t = 0; t < ths->d; t++) \
      { \
        fg_psi[t][0] = (PHI((2 * NN(ths->n[t])), (ths->x[j*ths->d+t] - ((R)u[t])/(2 * NN(ths->n[t]))),(t)));\
 \
        tmpEXP1 = EXP(K(2.0) * ((2 * NN(ths->n[t])) * ths->x[j * ths->d + t] - u[t]) / ths->b[t]); \
        tmp1 = K(1.0); \
        for (l_fg = u[t] + 1, lj_fg = 1; l_fg <= o[t]; l_fg++, lj_fg++) \
        { \
          tmp1 *= tmpEXP1; \
          fg_psi[t][lj_fg] = fg_psi[t][0] * tmp1 * fg_exp_l[t][lj_fg]; \
        } \
      } \
  \
      for (l_L = 0; l_L < lprod; l_L++) \
      { \
        MACRO_update_phi_prod_ll_plain(which_one, with_FG_PSI); \
 \
        MACRO_B_compute_ ## which_one; \
 \
        MACRO_count_uo_l_lj_t; \
      } \
    } \
    return; \
  } \
 \
  if (ths->flags & PRE_LIN_PSI) \
  { \
    for (j = 0, fj = f; j < ths->M_total; j++, fj++) \
    { \
      MACRO_init_uo_l_lj_t; \
  \
      for (t = 0; t < ths->d; t++) \
      { \
        y[t] = (((2 * NN(ths->n[t])) * ths->x[j * ths->d + t] - (R)u[t]) \
                * ((R)ths->K))/(ths->m + 2); \
        ip_u  = LRINT(FLOOR(y[t])); \
        ip_w  = y[t]-ip_u; \
        for (l_fg = u[t], lj_fg = 0; l_fg <= o[t]; l_fg++, lj_fg++) \
        { \
          fg_psi[t][lj_fg] = ths->psi[(ths->K+1)*t + ABS(ip_u-lj_fg*ip_s)] \
            * (1-ip_w) + ths->psi[(ths->K+1)*t + ABS(ip_u-lj_fg*ip_s+1)] \
            * (ip_w); \
        } \
      } \
  \
      for (l_L = 0; l_L < lprod; l_L++) \
      { \
        MACRO_update_phi_prod_ll_plain(which_one, with_FG_PSI); \
 \
        MACRO_B_compute_ ## which_one; \
 \
        MACRO_count_uo_l_lj_t; \
      }  /* for(l_L) */  \
    } /* for(j) */  \
    return; \
  } /* if(PRE_LIN_PSI) */ \
  \
  /* no precomputed psi at all */ \
  for (j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1) \
  { \
    MACRO_init_uo_l_lj_t; \
 \
    for (l_L = 0; l_L < lprod; l_L++) \
    { \
      MACRO_update_phi_prod_ll_plain(which_one, without_PRE_PSI); \
 \
      MACRO_B_compute_ ## which_one; \
 \
      MACRO_count_uo_l_lj_t; \
    } /* for (l_L) */ \
  } /* for (j) */ \
} /* B */

MACRO_B(A)
MACRO_B(T)

/**
 * user routines
 */
void X(trafo)(X(plan) *ths)
{
  switch(ths->d)
  {
    default:
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
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_r2r_plan);
      TOC_FFTW(1)

      /*if (ths->flags & PRE_FULL_PSI)
        full_psi__A(ths);*/

      /* Set \f$ f_j = \sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
       * \text{ for } j=0,\dots,M-1 \f$ */
      TIC(2)
      B_A(ths);
      TOC(2)

      /*if (ths->flags & PRE_FULL_PSI)
      {
        Y(free)(ths->psi_index_g);
        Y(free)(ths->psi_index_f);
      }*/
    }
  }
} /* trafo */

void X(adjoint)(X(plan) *ths)
{
  switch(ths->d)
  {
    default:
    {
      /* use ths->my_fftw_plan */
      ths->g_hat = ths->g2;
      ths->g = ths->g1;

      /*if (ths->flags & PRE_FULL_PSI)
        full_psi__T(ths);*/

      /* Set \f$ g_l = \sum_{j=0}^{M-1} f_j \psi\left(x_j-\frac{l}{n}\right)
       * \text{ for } l \in I_n,m(x_j) \f$ */
      TIC(2)
      B_T(ths);
      TOC(2)

      /* Compute by d-variate discrete cosine transform
       * \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
       * \text{ for }  k \in I_N\f$ */
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_r2r_plan);
      TOC_FFTW(1)

      /* Form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
       * k \in I_N \f$ */
      TIC(0)
      D_T(ths);
      TOC(0)
    }
  }
} /* adjoint */

/** initialisation of direct transform
 */
static inline void precompute_phi_hut(X(plan) *ths)
{
  INT ks[ths->d]; /* index over all frequencies */
  INT t; /* index over all dimensions */

  ths->c_phi_inv = (R**) Y(malloc)((size_t)(ths->d) * sizeof(R*));

  for (t = 0; t < ths->d; t++)
  {
    ths->c_phi_inv[t] = (R*)Y(malloc)((size_t)(ths->N[t] - OFFSET) * sizeof(R));

    for (ks[t] = 0; ks[t] < ths->N[t] - OFFSET; ks[t]++)
    {
      ths->c_phi_inv[t][ks[t]] = (K(1.0) / (PHI_HUT((2 * NN(ths->n[t])), ks[t] + OFFSET, t)));
    }
  }
} /* phi_hut */

/** create a lookup table, but NOT for each node
 *  good idea K=2^xx
 *  TODO: estimate K, call from init
 *  assumes an EVEN window function
 */
void X(precompute_lin_psi)(X(plan) *ths)
{
  INT t; /**< index over all dimensions */
  INT j; /**< index over all nodes */
  R step; /**< step size in [0,(m+2)/n] */

  for (t = 0; t < ths->d; t++)
  {
    step = ((R)(ths->m+2)) / (((R)ths->K) * (2 * NN(ths->n[t])));

    for (j = 0; j <= ths->K; j++)
    {
      ths->psi[(ths->K + 1) * t + j] = PHI((2 * NN(ths->n[t])), (j * step), t);
    } /* for(j) */
  } /* for(t) */
}

void X(precompute_fg_psi)(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT u, o; /* depends on x_j */

//  sort(ths);

  for (t = 0; t < ths->d; t++)
  {
    INT j;
//    #pragma omp parallel for default(shared) private(j,u,o)
    for (j = 0; j < ths->M_total; j++)
    {
      uo(ths, j, &u, &o, t);

      ths->psi[2 * (j*ths->d + t)] = (PHI((2 * NN(ths->n[t])),(ths->x[j * ths->d + t] - ((R)u) / (2 * NN(ths->n[t]))),(t)));
      ths->psi[2 * (j*ths->d + t) + 1] = EXP(K(2.0) * ( (2 * NN(ths->n[t])) * ths->x[j * ths->d + t] - u) / ths->b[t]);
      } /* for(j) */
  }
  /* for(t) */
} /* nfft_precompute_fg_psi */

void X(precompute_psi)(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT lj; /* index 0<=lj<u+o+1 */
  INT u, o; /* depends on x_j */

  //sort(ths);

  for (t = 0; t < ths->d; t++)
  {
    INT j;

    for (j = 0; j < ths->M_total; j++)
    {
      uo(ths, j, &u, &o, t);

      for(lj = 0; lj < (2 * ths->m + 2); lj++)
        ths->psi[(j * ths->d + t) * (2 * ths->m + 2) + lj] =
            (PHI((2 * NN(ths->n[t])), ((ths->x[(j) * ths->d + (t)]) - ((R)(lj + u)) / (K(2.0) * ((R)NN(ths->n[t])))), t));
    } /* for (j) */
  } /* for (t) */
} /* precompute_psi */

void X(precompute_full_psi)(X(plan) *ths)
{
//#ifdef _OPENMP
//  sort(ths);
//
//  nfft_precompute_full_psi_omp(ths);
//#else
  INT t, t2; /* index over all dimensions */
  INT j; /* index over all nodes */
  INT l_L; /* plain index 0 <= l_L < lprod */
  INT l[ths->d]; /* multi index u<=l<=o */
  INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */
  INT ll_plain[ths->d+1]; /* postfix plain index */
  INT lprod; /* 'bandwidth' of matrix B */
  INT u[ths->d], o[ths->d]; /* depends on x_j */
  INT count_lg[ths->d];
  INT lg_offset[ths->d];

  R phi_prod[ths->d+1];

  INT ix, ix_old;

  //sort(ths);

  phi_prod[0] = K(1.0);
  ll_plain[0]  = 0;

  for (t = 0, lprod = 1; t < ths->d; t++)
    lprod *= 2 * ths->m + 2;

  for (j = 0, ix = 0, ix_old = 0; j < ths->M_total; j++)
  {
    MACRO_init_uo_l_lj_t;

    for (l_L = 0; l_L < lprod; l_L++, ix++)
    {
      MACRO_update_phi_prod_ll_plain(A, without_PRE_PSI);

      ths->psi_index_g[ix] = ll_plain[ths->d];
      ths->psi[ix] = phi_prod[ths->d];

      MACRO_count_uo_l_lj_t;
    } /* for (l_L) */

    ths->psi_index_f[j] = ix - ix_old;
    ix_old = ix;
  } /* for(j) */
//#endif
}

void X(precompute_one_psi)(X(plan) *ths)
{
  if(ths->flags & PRE_PSI)
    X(precompute_psi)(ths);
  if(ths->flags & PRE_FULL_PSI)
    X(precompute_full_psi)(ths);
  if(ths->flags & PRE_FG_PSI)
    X(precompute_fg_psi)(ths);
  if(ths->flags & PRE_LIN_PSI)
    X(precompute_lin_psi)(ths);
}

static inline void init_help(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT lprod; /* 'bandwidth' of matrix B */

  if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    ths->flags |= NFFT_SORT_NODES;

  ths->N_total = intprod(ths->N, OFFSET, ths->d);
  ths->n_total = intprod(ths->n, 0, ths->d);

  ths->sigma = (R*)Y(malloc)((size_t)(ths->d) * sizeof(R));

  for (t = 0; t < ths->d; t++)
    ths->sigma[t] = ((R)NN(ths->n[t])) / ths->N[t];

  /* Assign r2r transform kinds for each dimension */
  ths->r2r_kind = (FFTW(r2r_kind)*)Y(malloc)((size_t)(ths->d) * sizeof (FFTW(r2r_kind)));
  for (t = 0; t < ths->d; t++)
    ths->r2r_kind[t] = FOURIER_TRAFO;

  WINDOW_HELP_INIT;

  if (ths->flags & MALLOC_X)
    ths->x = (R*)Y(malloc)((size_t)(ths->d * ths->M_total) * sizeof(R));

  if (ths->flags & MALLOC_F_HAT)
    ths->f_hat = (R*)Y(malloc)((size_t)(ths->N_total) * sizeof(R));

  if (ths->flags & MALLOC_F)
    ths->f = (R*)Y(malloc)((size_t)(ths->M_total) * sizeof(R));

  if (ths->flags & PRE_PHI_HUT)
    precompute_phi_hut(ths);

  if(ths->flags & PRE_LIN_PSI)
  {
      ths->K = (1U<< 10) * (ths->m+2);
      ths->psi = (R*) Y(malloc)((size_t)((ths->K + 1) * ths->d) * sizeof(R));
  }

  if(ths->flags & PRE_FG_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * 2) * sizeof(R));

  if (ths->flags & PRE_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * (2 * ths->m + 2 )) * sizeof(R));

  if(ths->flags & PRE_FULL_PSI)
  {
      for (t = 0, lprod = 1; t < ths->d; t++)
        lprod *= 2 * ths->m + 2;

      ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(R));

      ths->psi_index_f = (INT*) Y(malloc)((size_t)(ths->M_total) * sizeof(INT));
      ths->psi_index_g = (INT*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(INT));
  }

  if (ths->flags & FFTW_INIT)
  {
    ths->g1 = (R*)Y(malloc)((size_t)(ths->n_total) * sizeof(R));

    if (ths->flags & FFT_OUT_OF_PLACE)
      ths->g2 = (R*) Y(malloc)((size_t)(ths->n_total) * sizeof(R));
    else
      ths->g2 = ths->g1;

    {
      int *_n = Y(malloc)((size_t)(ths->d) * sizeof(int));

      for (t = 0; t < ths->d; t++)
        _n[t] = (int)(ths->n[t]);

      ths->my_fftw_r2r_plan = FFTW(plan_r2r)((int)ths->d, _n, ths->g1, ths->g2, ths->r2r_kind, ths->fftw_flags);
      Y(free)(_n);
    }
  }

//  if(ths->flags & NFFT_SORT_NODES)
//    ths->index_x = (INT*) Y(malloc)(sizeof(INT)*2*ths->M_total);
//  else
//    ths->index_x = NULL;

  ths->mv_trafo = (void (*) (void* ))X(trafo);
  ths->mv_adjoint = (void (*) (void* ))X(adjoint);
}

void X(init)(X(plan) *ths, int d, int *N, int M_total)
{
  int t; /* index over all dimensions */

  ths->d = (INT)d;

  ths->N = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->M_total = (INT)M_total;

  ths->n = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = 2 * (Y(next_power_of_2)(ths->N[t]) - 1) + OFFSET;

  ths->m = WINDOW_HELP_ESTIMATE_m;

  if (d > 1)
  {
//#ifdef _OPENMP
//    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
//                      FFTW_INIT | NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT;
//#else
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | NFFT_SORT_NODES;
//#endif
  }
  else
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE;

  ths->fftw_flags = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;

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

  for (j = 0; j < ths->M_total * ths->d; j++)
  {
    if ((ths->x[j] < K(0.0)) || (ths->x[j] >= K(0.5)))
    {
      return "ths->x out of range [0.0,0.5)";
    }
  }

  for (j = 0; j < ths->d; j++)
  {
    if (ths->sigma[j] <= 1)
      return "Oversampling factor too small";

    if(ths->N[j] - 1 <= ths->m)
      return "Polynomial degree N is smaller than cut-off m";
  }
  return 0;
}

void X(finalize)(X(plan) *ths)
{
  INT t; /* index over dimensions */

//  if(ths->flags & NFFT_SORT_NODES)
//    Y(free)(ths->index_x);

  if (ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_r2r_plan);

    if (ths->flags & FFT_OUT_OF_PLACE)
      Y(free)(ths->g2);

    Y(free)(ths->g1);
  }

  if(ths->flags & PRE_FULL_PSI)
  {
    Y(free)(ths->psi_index_g);
    Y(free)(ths->psi_index_f);
    Y(free)(ths->psi);
  }

  if (ths->flags & PRE_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_FG_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_LIN_PSI)
    Y(free)(ths->psi);

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
} /* finalize */
