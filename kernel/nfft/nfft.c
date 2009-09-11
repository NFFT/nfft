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

/**
 * Simple and fast computation of the NDFT.
 * authors: D. Potts, S. Kunis 2002-2006
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"

/** direct computation of non equispaced fourier transforms
 *  ndft_trafo, ndft_conjugated, ndft_adjoint, ndft_transposed
 *  require O(M_total N^d) arithemtical operations
 *
 * direct computation of the ndft_trafo and ndft_conjugated, formula (1.1)
 * ndft_trafo:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 * ndft_conjugated:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 *
 * direct computation of the ndft_adjoint and ndft_transposed, formula (1.2)
 * ndft_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(+2(pi) k x[j])
 * ndft_transposed:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(-2(pi) k x[j])
 */
/** macros and small sub routines for the direct transforms
 */
#define MACRO_ndft_init_result_trafo memset(f,0,ths->M_total*                 \
                                            sizeof(double _Complex));
#define MACRO_ndft_init_result_conjugated MACRO_ndft_init_result_trafo
#define MACRO_ndft_init_result_adjoint memset(f_hat,0,ths->N_total*           \
					      sizeof(double _Complex));
#define MACRO_ndft_init_result_transposed MACRO_ndft_init_result_adjoint

#define MACRO_ndft_sign_trafo      +2*PI*ths->x[j*ths->d+t]
#define MACRO_ndft_sign_conjugated -2*PI*ths->x[j*ths->d+t]
#define MACRO_ndft_sign_adjoint    +2*PI*ths->x[j*ths->d+t]
#define MACRO_ndft_sign_transposed -2*PI*ths->x[j*ths->d+t]

#define MACRO_init_k_N_Omega_x(which_one) {                                   \
for(t=0; t<ths->d; t++)                                                       \
  {                                                                           \
    k[t]=-ths->N[t]/2;                                                        \
    x[t]= MACRO_ndft_sign_ ## which_one;                                      \
    Omega[t+1]=k[t]*x[t]+Omega[t];                                            \
  }                                                                           \
omega=Omega[ths->d];                                                          \
}                                                                             \

#define MACRO_count_k_N_Omega {                                               \
for(t = ths->d-1; (t >= 1) && (k[t] == ths->N[t]/2-1); t--)                   \
  k[t]-= ths->N[t]-1;                                                         \
                                                                              \
k[t]++;                                                                       \
                                                                              \
for(t2 = t; t2<ths->d; t2++)                                                  \
  Omega[t2+1]=k[t2]*x[t2]+Omega[t2];                                          \
                                                                              \
omega=Omega[ths->d];                                                          \
}

#define MACRO_ndft_compute_trafo (*fj) += (*f_hat_k)*cexp(-_Complex_I*omega);

#define MACRO_ndft_compute_conjugated MACRO_ndft_compute_trafo

#define MACRO_ndft_compute_adjoint (*f_hat_k) += (*fj)*cexp(+ _Complex_I*omega);

#define MACRO_ndft_compute_transposed MACRO_ndft_compute_adjoint

#define MACRO_ndft(which_one)                                                 \
void ndft_ ## which_one (nfft_plan *ths)                                      \
{                                                                             \
  int j;                                /**< index over all nodes           */\
  int t,t2;                             /**< index for dimensions           */\
  int k_L;                              /**< plain index for summation      */\
  double _Complex *f_hat, *f;           /**< dito                           */\
  double _Complex *f_hat_k;             /**< actual Fourier coefficient     */\
  double _Complex *fj;                  /**< actual sample                  */\
  double x[ths->d];                     /**< actual node x[d*j+t]           */\
  int k[ths->d];                        /**< multi index for summation      */\
  double omega, Omega[ths->d+1];        /**< sign times 2*pi*k*x            */\
                                                                              \
  f_hat=ths->f_hat; f=ths->f;                                                 \
                                                                              \
  MACRO_ndft_init_result_ ## which_one                                        \
                                                                              \
  if(ths->d==1) /* univariate case (due to performance) */                    \
    {                                                                         \
      t=0;                                                                    \
      for(j=0, fj = f; j<ths->M_total; j++, fj++)                             \
        {                                                                     \
	  for(k_L=0, f_hat_k = f_hat; k_L<ths->N_total; k_L++, f_hat_k++)     \
	    {                                                                 \
	      omega=(k_L-ths->N_total/2)* MACRO_ndft_sign_ ## which_one;      \
              MACRO_ndft_compute_ ## which_one;                               \
	    }                                                                 \
        }                                                                     \
    }                                                                         \
  else /* multivariate case */					              \
    {                                                                         \
      Omega[0]=0;                                                             \
      for(j=0, fj=f; j<ths->M_total; j++, fj++)                               \
        {                                                                     \
          MACRO_init_k_N_Omega_x(which_one);                                  \
          for(k_L=0, f_hat_k=f_hat; k_L<ths->N_total; k_L++, f_hat_k++)       \
	    {                                                                 \
              MACRO_ndft_compute_ ## which_one;                               \
	      MACRO_count_k_N_Omega;                                          \
	    } /* for(k_L) */                                                  \
        } /* for(j) */                                                        \
    } /* else */                                                              \
} /* ndft_trafo */


/** user routines
 */
MACRO_ndft(trafo)
MACRO_ndft(adjoint)

/** fast computation of non equispaced fourier transforms
 *  require O(N^d log(N) + M_total) arithemtical operations
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
static void nfft_uo(const nfft_plan *ths,const int j,int *up,int *op,const int act_dim)
{
  double xj=ths->x[j*ths->d+act_dim];
  int c = LRINT(xj * ths->n[act_dim]);

  if(xj < 0)
    {
      (*up) = c-1-(ths->m);
      (*op) = c  +(ths->m);
    }
  else
    {
      (*up) = c  -(ths->m);
      (*op) = c+1+(ths->m);
    }
}

static void nfft_uo2(int *u, int *o, const double x, const int n, const int m)
{
  int c = LRINT(x * n);

  if(x < 0)
    {
      *u=(c-1-m+n)%n;
      *o=(c+  m+n)%n;
    }
  else
    {
      *u=(c  -m+n)%n;
      *o=(c+1+m+n)%n;
    }
}


#define MACRO_nfft_D_compute_A {                                              \
 g_hat[k_plain[ths->d]] = f_hat[ks_plain[ths->d]] * c_phi_inv_k[ths->d];      \
}

#define MACRO_nfft_D_compute_T {                                              \
 f_hat[ks_plain[ths->d]] = g_hat[k_plain[ths->d]] * c_phi_inv_k[ths->d];      \
}

#define MACRO_nfft_D_init_result_A  memset(g_hat,0,ths->n_total*              \
					   sizeof(double _Complex));
#define MACRO_nfft_D_init_result_T memset(f_hat,0,ths->N_total*               \
                                           sizeof(double _Complex));

#define MACRO_with_PRE_PHI_HUT * ths->c_phi_inv[t2][ks[t2]];
#define MACRO_without_PRE_PHI_HUT / (PHI_HUT(ks[t2]-ths->N[t2]/2,t2));

#define MACRO_init_k_ks {                                                     \
  for(t = ths->d-1; t>=0; t--)                                                \
    {                                                                         \
      kp[t]= 0;                                                               \
      k[t] = 0;                                                               \
      ks[t] = ths->N[t]/2;                                                    \
    }                                                                         \
  t++;                                                                        \
}

#define MACRO_update_c_phi_inv_k(which_one) {                                 \
  for(t2=t; t2<ths->d; t2++)                                                  \
    {                                                                         \
      c_phi_inv_k[t2+1]= c_phi_inv_k[t2] MACRO_ ##which_one;                  \
      ks_plain[t2+1]= ks_plain[t2]*ths->N[t2]+ks[t2];                         \
      k_plain[t2+1]= k_plain[t2]*ths->n[t2]+k[t2];                            \
    }                                                                         \
}

#define MACRO_count_k_ks {                                                    \
  for(t=ths->d-1; (t>0)&& (kp[t]==ths->N[t]-1); t--)                          \
    {                                                                         \
      kp[t]= 0;                                                               \
      k[t]= 0;                                                                \
      ks[t]= ths->N[t]/2;                                                     \
    }                                                                         \
                                                                              \
  kp[t]++; k[t]++; ks[t]++;                                                   \
  if(kp[t]==ths->N[t]/2)                                                      \
    {                                                                         \
      k[t]= ths->n[t]-ths->N[t]/2;                                            \
      ks[t]= 0;                                                               \
    }                                                                         \
}                                                                             \


/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$D, D^T\f$
 */
#define MACRO_nfft_D(which_one)                                               \
static inline void nfft_D_ ## which_one (nfft_plan *ths)                      \
{                                                                             \
  int t, t2;                            /**< index dimensions               */\
  int k_L;                              /**< plain index                    */\
  int kp[ths->d];                       /**< multi index (simple)           */\
  int k[ths->d];                        /**< multi index in g_hat           */\
  int ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/\
  double c_phi_inv_k[ths->d+1];         /**< postfix product of PHI_HUT     */\
  int k_plain[ths->d+1];                /**< postfix plain index            */\
  int ks_plain[ths->d+1];               /**< postfix plain index            */\
  double _Complex *f_hat, *g_hat;       /**< local copy                     */\
                                                                              \
  f_hat=ths->f_hat; g_hat=ths->g_hat;                                         \
  MACRO_nfft_D_init_result_ ## which_one;                                     \
                                                                              \
  c_phi_inv_k[0]=1;                                                           \
  k_plain[0]=0;                                                               \
  ks_plain[0]=0;                                                              \
                                                                              \
  if(ths->nfft_flags & PRE_PHI_HUT)                                           \
    {                                                                         \
      MACRO_init_k_ks;                                                        \
                                                                              \
      for(k_L=0; k_L<ths->N_total; k_L++)                                     \
	{                                                                     \
          MACRO_update_c_phi_inv_k(with_PRE_PHI_HUT);                         \
                                                                              \
	  MACRO_nfft_D_compute_ ## which_one;                                 \
	                                                                      \
	  MACRO_count_k_ks;                                                   \
	} /* for(k_L) */                                                      \
    } /* if(PRE_PHI_HUT) */                                                   \
  else                                                                        \
    {                                                                         \
      MACRO_init_k_ks;                                                        \
                                                                              \
      for(k_L=0; k_L<ths->N_total; k_L++)                                     \
	  {                                                                   \
        MACRO_update_c_phi_inv_k(without_PRE_PHI_HUT);                        \
                                                                              \
	    MACRO_nfft_D_compute_ ## which_one;                               \
	                                                                      \
	    MACRO_count_k_ks;                                                 \
	  } /* for(k_L) */                                                    \
    } /* else(PRE_PHI_HUT) */                                                 \
} /* nfft_D */

MACRO_nfft_D(A)
MACRO_nfft_D(T)

/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$B, B^{\rm T}\f$
 */
#define MACRO_nfft_B_init_result_A  memset(f,0,ths->M_total*                  \
                                           sizeof(double _Complex));
#define MACRO_nfft_B_init_result_T memset(g,0,ths->n_total*                   \
                                          sizeof(double _Complex));

#define MACRO_nfft_B_PRE_FULL_PSI_compute_A {                                 \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]];			      \
}

#define MACRO_nfft_B_PRE_FULL_PSI_compute_T {                                 \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj);                            \
}

#define MACRO_nfft_B_compute_A {                                              \
  (*fj) += phi_prod[ths->d] * g[ll_plain[ths->d]];                            \
}

#define MACRO_nfft_B_compute_T {                                              \
  g[ll_plain[ths->d]] += phi_prod[ths->d] * (*fj);                            \
}

#define MACRO_with_FG_PSI fg_psi[t2][lj[t2]]

#define MACRO_with_PRE_PSI     ths->psi[(j*ths->d+t2)*(2*ths->m+2)+lj[t2]]

#define MACRO_without_PRE_PSI  PHI(ths->x[j*ths->d+t2]-                       \
                                   ((double)l[t2])/ths->n[t2], t2)

#define MACRO_init_uo_l_lj_t {                                                \
  for(t = ths->d-1; t>=0; t--)                                                \
    {                                                                         \
      nfft_uo(ths,j,&u[t],&o[t],t);                                           \
      l[t] = u[t];                                                            \
      lj[t] = 0;                                                              \
    } /* for(t) */                                                            \
  t++;                                                                        \
}

#define MACRO_update_phi_prod_ll_plain(which_one) {                           \
  for(t2=t; t2<ths->d; t2++)                                                  \
    {                                                                         \
      phi_prod[t2+1]=phi_prod[t2]* MACRO_ ## which_one;                       \
      ll_plain[t2+1]=ll_plain[t2]*ths->n[t2] +(l[t2]+ths->n[t2])%ths->n[t2];  \
    } /* for(t2) */                                                           \
}

#define MACRO_count_uo_l_lj_t {                                               \
  for(t = ths->d-1; (t>0)&&(l[t]==o[t]); t--)                                 \
    {                                                                         \
      l[t] = u[t];                                                            \
      lj[t] = 0;                                                              \
    } /* for(t) */                                                            \
                                                                              \
  l[t]++;                                                                     \
  lj[t]++;                                                                    \
}

#define MACRO_nfft_B(which_one)                                               \
static inline void nfft_B_ ## which_one (nfft_plan *ths)                      \
{                                                                             \
  int lprod;                            /**< 'regular bandwidth' of matrix B*/\
  int u[ths->d], o[ths->d];             /**< multi band with respect to x_j */\
  int t, t2;                            /**< index dimensions               */\
  int j;                                /**< index nodes                    */\
  int l_L, ix;                          /**< index one row of B             */\
  int l[ths->d];                        /**< multi index u<=l<=o            */\
  int lj[ths->d];                       /**< multi index 0<=lj<u+o+1        */\
  int ll_plain[ths->d+1];               /**< postfix plain index in g       */\
  double phi_prod[ths->d+1];            /**< postfix product of PHI         */\
  double _Complex *f, *g;               /**< local copy                     */\
  double _Complex *fj;                  /**< local copy                     */\
  double y[ths->d];                                                           \
  double fg_psi[ths->d][2*ths->m+2];                                          \
  double fg_exp_l[ths->d][2*ths->m+2];                                        \
  int l_fg,lj_fg;                                                             \
  double tmpEXP1, tmpEXP2, tmpEXP2sq, tmp1, tmp2, tmp3;                       \
  double ip_w;                                                                \
  int ip_u;                                                                   \
  int ip_s=ths->K/(ths->m+1);                                                 \
                                                                              \
  f=ths->f; g=ths->g;                                                         \
                                                                              \
  MACRO_nfft_B_init_result_ ## which_one;                                     \
                                                                              \
  if(ths->nfft_flags & PRE_FULL_PSI)                                          \
    {                                                                         \
      for(ix=0, j=0, fj=f; j<ths->M_total; j++, fj++)                         \
        for(l_L=0; l_L<ths->psi_index_f[j]; l_L++, ix++)                      \
	  MACRO_nfft_B_PRE_FULL_PSI_compute_ ## which_one;                    \
      return;                                                                 \
    }                                                                         \
                                                                              \
  phi_prod[0]=1;                                                              \
  ll_plain[0]=0;                                                              \
                                                                              \
  for(t=0,lprod = 1; t<ths->d; t++)                                           \
    lprod *= (2*ths->m+2);                                                    \
                                                                              \
  if(ths->nfft_flags & PRE_PSI)                                               \
    {                                                                         \
      for(j=0, fj=f; j<ths->M_total; j++, fj++)                               \
	{                                                                     \
          MACRO_init_uo_l_lj_t;                                               \
                                                                              \
	  for(l_L=0; l_L<lprod; l_L++)                                        \
	    {                                                                 \
              MACRO_update_phi_prod_ll_plain(with_PRE_PSI);                   \
                                                                              \
	      MACRO_nfft_B_compute_ ## which_one;                             \
		                                                              \
	      MACRO_count_uo_l_lj_t;                                          \
            } /* for(l_L) */                                                  \
	} /* for(j) */                                                        \
      return;                                                                 \
    } /* if(PRE_PSI) */                                                       \
                                                                              \
  if(ths->nfft_flags & PRE_FG_PSI)                                            \
    {                                                                         \
      for(t2=0; t2<ths->d; t2++)                                              \
        {                                                                     \
          tmpEXP2 = exp(-1.0/ths->b[t2]);                                     \
          tmpEXP2sq = tmpEXP2*tmpEXP2;                                        \
          tmp2 = 1.0;                                                         \
          tmp3 = 1.0;                                                         \
          fg_exp_l[t2][0] = 1.0;                                              \
          for(lj_fg=1; lj_fg <= (2*ths->m+2); lj_fg++)                        \
            {                                                                 \
              tmp3 = tmp2*tmpEXP2;                                            \
              tmp2 *= tmpEXP2sq;                                              \
              fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;               \
            }                                                                 \
        }                                                                     \
      for(j=0, fj=f; j<ths->M_total; j++, fj++)                               \
	{                                                                     \
          MACRO_init_uo_l_lj_t;                                               \
                                                                              \
          for(t2=0; t2<ths->d; t2++)                                          \
            {                                                                 \
              fg_psi[t2][0] = ths->psi[2*(j*ths->d+t2)];                      \
              tmpEXP1 = ths->psi[2*(j*ths->d+t2)+1];                          \
              tmp1 = 1.0;                                                     \
              for(l_fg=u[t2]+1, lj_fg=1; l_fg <= o[t2]; l_fg++, lj_fg++)      \
                {                                                             \
                  tmp1 *= tmpEXP1;                                            \
                  fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
                }                                                             \
            }                                                                 \
                                                                              \
	  for(l_L=0; l_L<lprod; l_L++)                                        \
	    {                                                                 \
              MACRO_update_phi_prod_ll_plain(with_FG_PSI);                    \
                                                                              \
	      MACRO_nfft_B_compute_ ## which_one;                             \
		                                                              \
	      MACRO_count_uo_l_lj_t;                                          \
            } /* for(l_L) */                                                  \
	} /* for(j) */                                                        \
      return;                                                                 \
    } /* if(PRE_FG_PSI) */                                                    \
                                                                              \
  if(ths->nfft_flags & FG_PSI)                                                \
    {                                                                         \
      for(t2=0; t2<ths->d; t2++)                                              \
        {                                                                     \
          tmpEXP2 = exp(-1.0/ths->b[t2]);                                     \
          tmpEXP2sq = tmpEXP2*tmpEXP2;                                        \
          tmp2 = 1.0;                                                         \
          tmp3 = 1.0;                                                         \
          fg_exp_l[t2][0] = 1.0;                                              \
          for(lj_fg=1; lj_fg <= (2*ths->m+2); lj_fg++)                        \
            {                                                                 \
              tmp3 = tmp2*tmpEXP2;                                            \
              tmp2 *= tmpEXP2sq;                                              \
              fg_exp_l[t2][lj_fg] = fg_exp_l[t2][lj_fg-1]*tmp3;               \
            }                                                                 \
        }                                                                     \
      for(j=0, fj=f; j<ths->M_total; j++, fj++)                               \
	{                                                                     \
          MACRO_init_uo_l_lj_t;                                               \
                                                                              \
          for(t2=0; t2<ths->d; t2++)                                          \
            {                                                                 \
              fg_psi[t2][0] =                                                 \
                (PHI((ths->x[j*ths->d+t2]-((double)u[t2])/ths->n[t2]),t2));   \
                                                                              \
              tmpEXP1 = exp(2.0*(ths->n[t2]*ths->x[j*ths->d+t2] - u[t2])      \
                      / ths->b[t2]);                                          \
              tmp1 = 1.0;                                                     \
              for(l_fg=u[t2]+1, lj_fg=1; l_fg <= o[t2]; l_fg++, lj_fg++)      \
                {                                                             \
                  tmp1 *= tmpEXP1;                                            \
                  fg_psi[t2][lj_fg] = fg_psi[t2][0]*tmp1*fg_exp_l[t2][lj_fg]; \
                }                                                             \
            }                                                                 \
                                                                              \
	  for(l_L=0; l_L<lprod; l_L++)                                        \
	    {                                                                 \
              MACRO_update_phi_prod_ll_plain(with_FG_PSI);                    \
                                                                              \
	      MACRO_nfft_B_compute_ ## which_one;                             \
		                                                              \
	      MACRO_count_uo_l_lj_t;                                          \
            } /* for(l_L) */                                                  \
	} /* for(j) */                                                        \
      return;                                                                 \
    } /* if(FG_PSI) */                                                        \
                                                                              \
                                                                              \
  if(ths->nfft_flags & PRE_LIN_PSI)                                           \
    {                                                                         \
      for(j=0, fj=f; j<ths->M_total; j++, fj++)                               \
	{                                                                     \
          MACRO_init_uo_l_lj_t;                                               \
                                                                              \
          for(t2=0; t2<ths->d; t2++)                                          \
            {                                                                 \
              y[t2] = ((ths->n[t2]*ths->x[j*ths->d+t2]-                       \
                          (double)u[t2]) * ((double)ths->K))/(ths->m+1);      \
              ip_u  = LRINT(floor(y[t2]));                                    \
              ip_w  = y[t2]-ip_u;                                             \
              for(l_fg=u[t2], lj_fg=0; l_fg <= o[t2]; l_fg++, lj_fg++)        \
                {                                                             \
                  fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2+                 \
					       abs(ip_u-lj_fg*ip_s)]*         \
                                      (1-ip_w) +                              \
                                      ths->psi[(ths->K+1)*t2+                 \
					       abs(ip_u-lj_fg*ip_s+1)]*       \
                                       (ip_w);                                \
              }                                                               \
            }                                                                 \
                                                                              \
	  for(l_L=0; l_L<lprod; l_L++)                                        \
	    {                                                                 \
              MACRO_update_phi_prod_ll_plain(with_FG_PSI);                    \
                                                                              \
	      MACRO_nfft_B_compute_ ## which_one;                             \
		                                                              \
	      MACRO_count_uo_l_lj_t;                                          \
            } /* for(l_L) */                                                  \
	} /* for(j) */                                                        \
      return;                                                                 \
    } /* if(PRE_LIN_PSI) */                                                   \
                                                                              \
  /* no precomputed psi at all */                                             \
  for(j=0, fj=f; j<ths->M_total; j++, fj++)                                   \
    {                                                                         \
      MACRO_init_uo_l_lj_t;                                                   \
	                                                                      \
      for(l_L=0; l_L<lprod; l_L++)                                            \
     	{                                                                     \
          MACRO_update_phi_prod_ll_plain(without_PRE_PSI);                    \
                                                                              \
          MACRO_nfft_B_compute_ ## which_one;                                 \
		                                                              \
          MACRO_count_uo_l_lj_t;                                              \
	} /* for(l_L) */                                                      \
    } /* for(j) */                                                            \
} /* nfft_B */                                                                \

MACRO_nfft_B(A)
MACRO_nfft_B(T)

/* ############################################################ SPECIFIC VERSIONS FOR d=1 */

static void nfft_1d_init_fg_exp_l(double *fg_exp_l, const int m, const double b)
{
  int l;
  double fg_exp_b0, fg_exp_b1, fg_exp_b2, fg_exp_b0_sq;

  fg_exp_b0 = exp(-1.0/b);
  fg_exp_b0_sq = fg_exp_b0*fg_exp_b0;
  fg_exp_b1 = 1.0;
  fg_exp_b2 = 1.0;
  fg_exp_l[0] = 1.0;
  for(l=1; l <= 2*m+1; l++)
    {
      fg_exp_b2 = fg_exp_b1*fg_exp_b0;
      fg_exp_b1 *= fg_exp_b0_sq;
      fg_exp_l[l] = fg_exp_l[l-1]*fg_exp_b2;
    }
}

static void nfft_trafo_1d_compute(double _Complex *fj, const double _Complex *g,const double *psij_const, const double *xj, const int n, const int m)
{
  int u,o,l;
  const double _Complex *gj;
  const double *psij;
  psij=psij_const;

  nfft_uo2(&u,&o,*xj, n, m);

  if(u<o)
    for(l=1,gj=g+u,(*fj)=(*psij++) * (*gj++); l<=2*m+1; l++)
      (*fj) += (*psij++) * (*gj++);
  else
    {
      for(l=1,gj=g+u,(*fj)=(*psij++) * (*gj++); l<2*m+1-o; l++)
	(*fj) += (*psij++) * (*gj++);
      for(l=0,gj=g; l<=o; l++)
	(*fj) += (*psij++) * (*gj++);
    }
}

static void nfft_adjoint_1d_compute(const double _Complex *fj, double _Complex *g,const double *psij_const, const double *xj, const int n, const int m)
{
  int u,o,l;
  double _Complex *gj;
  const double *psij;
  psij=psij_const;

  nfft_uo2(&u,&o,*xj, n, m);
  
  if(u<o)
    for(l=0,gj=g+u; l<=2*m+1; l++)
      (*gj++) += (*psij++) * (*fj);
  else
    {
      for(l=0,gj=g+u; l<2*m+1-o; l++)
	(*gj++) += (*psij++) * (*fj);
      for(l=0,gj=g; l<=o; l++)
	(*gj++) += (*psij++) * (*fj);
    }
}

static void nfft_trafo_1d_B(nfft_plan *ths)
{
  int n,N,u,o,j,M,l,m, *psi_index_g,K,ip_s,ip_u;
  double _Complex *fj,*g;
  double *psij, *psij_const, *xj, ip_y,ip_w;

  double *fg_exp_l, fg_psij0, fg_psij1, fg_psij2;

  N=ths->N[0];
  n=ths->n[0];
  M=ths->M_total;
  m=ths->m;

  g=ths->g;


  if(ths->nfft_flags & PRE_FULL_PSI)
    {
      psi_index_g=ths->psi_index_g;
      for(j=0, fj=ths->f, psij=ths->psi; j<M; j++, fj++)
        for(l=1, (*fj)=(*psij++) * g[(*psi_index_g++)]; l<=2*m+1; l++)
	  (*fj) += (*psij++) * g[(*psi_index_g++)];
      return;
    } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
    {
      for(j=0,fj=ths->f,xj=ths->x; j<M; j++,fj++,xj++)
	nfft_trafo_1d_compute(fj, g, ths->psi+j*(2*m+2), xj, n, m);
      return;
    } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
    {
      psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc((2*m+2)*sizeof(double));

      nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
	{
	  fg_psij0 = ths->psi[2*j];
	  fg_psij1 = ths->psi[2*j+1];
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_trafo_1d_compute(fj, g, psij_const, xj, n, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
    {
      psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc((2*m+2)*sizeof(double));

      nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  fg_psij0 = (PHI(*xj-((double)u)/n,0));
	  fg_psij1 = exp(2.0*(n*(*xj) - u)/ths->b[0]);
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_trafo_1d_compute(fj, g, psij_const, xj, n, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
    {
      psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
      K=ths->K;
      ip_s=K/(m+1);

      psij_const[2*m+1]=0;

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
	{
	  nfft_uo(ths,j,&u,&o,0);

	  ip_y = fabs(n*(*xj) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[l] = ths->psi[abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_trafo_1d_compute(fj, g, psij_const, xj, n, m);
	}
      nfft_free(psij_const);
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
  for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
    {
      nfft_uo(ths,j,&u,&o,0);

      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(*xj-((double)((u+l)))/n,0));

      nfft_trafo_1d_compute(fj, g, psij_const, xj, n, m);
    }
  nfft_free(psij_const);
}

static void nfft_adjoint_1d_B(nfft_plan *ths)
{
  int n,u,o,j,M,l,m, *psi_index_g,K,ip_s,ip_u;
  double _Complex *fj,*g;
  double *psij, *psij_const, *xj, ip_y, ip_w;
  double *fg_exp_l, fg_psij0, fg_psij1, fg_psij2;

  n=ths->n[0];
  M=ths->M_total;
  m=ths->m;

  g=ths->g;
  memset(g,0,ths->n_total*sizeof(double _Complex));

  if(ths->nfft_flags & PRE_FULL_PSI)
    {
      psi_index_g=ths->psi_index_g;
      for(j=0, fj=ths->f, psij=ths->psi; j<M; j++, fj++)
        for(l=0; l<=2*m+1; l++)
	  g[*psi_index_g++] += (*psij++) * (*fj);
      return;
    } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
    {
      for(j=0,fj=ths->f,xj=ths->x; j<M; j++,fj++,xj++)
	nfft_adjoint_1d_compute(fj, g, ths->psi+j*(2*m+2), xj, n, m);
      return;
    } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
    {
      psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc((2*m+2)*sizeof(double));

      nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
	{
	  fg_psij0 = ths->psi[2*j];
	  fg_psij1 = ths->psi[2*j+1];
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_adjoint_1d_compute(fj, g, psij_const, xj, n, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
    {
      psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc((2*m+2)*sizeof(double));

      nfft_1d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  fg_psij0 = (PHI(*xj-((double)u)/n,0));
	  fg_psij1 = exp(2.0*(n*(*xj) - u)/ths->b[0]);
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_adjoint_1d_compute(fj, g, psij_const, xj, n, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
    {
      psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
      K=ths->K;
      ip_s=K/(m+1);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
	{
	  nfft_uo(ths,j,&u,&o,0);

	  ip_y = fabs(n*(*xj) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[l] = ths->psi[abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_adjoint_1d_compute(fj, g, psij_const, xj, n, m);
	}
      nfft_free(psij_const);
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  psij_const=(double*)nfft_malloc((2*m+2)*sizeof(double));
  for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj++)
    {
      nfft_uo(ths,j,&u,&o,0);

      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(*xj-((double)((u+l)))/n,0));

      nfft_adjoint_1d_compute(fj, g, psij_const, xj, n, m);
    }
  nfft_free(psij_const);
}

void nfft_trafo_1d(nfft_plan *ths)
{
  int k,n,N;
  double _Complex *g_hat1,*g_hat2,*f_hat1,*f_hat2;
  double *c_phi_inv1, *c_phi_inv2;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N=ths->N[0];
  n=ths->n[0];

  f_hat1=ths->f_hat;
  f_hat2=&ths->f_hat[N/2];
  g_hat1=&ths->g_hat[n-N/2];
  g_hat2=ths->g_hat;

  TIC(0)
  memset(ths->g_hat,0,ths->n_total*sizeof(double _Complex));
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv1=ths->c_phi_inv[0];
      c_phi_inv2=&ths->c_phi_inv[0][N/2];
      for(k=0;k<N/2;k++)
	{
	  (*g_hat1++) = (*f_hat1++) * (*c_phi_inv1++);
	  (*g_hat2++) = (*f_hat2++) * (*c_phi_inv2++);
	}
    }
  else
    for(k=0;k<N/2;k++)
      {
	(*g_hat1++) = (*f_hat1++) / (PHI_HUT(k-N/2,0));
	(*g_hat2++) = (*f_hat2++) / (PHI_HUT(k,0));
      }

  TOC(0)

  TIC_FFTW(1)
  fftw_execute(ths->my_fftw_plan1);
  TOC_FFTW(1);

  TIC(2);
  nfft_trafo_1d_B(ths);
  TOC(2);
}

void nfft_adjoint_1d(nfft_plan *ths)
{
  int k,n,N;
  double _Complex *g_hat1,*g_hat2,*f_hat1,*f_hat2;
  double *c_phi_inv1, *c_phi_inv2;

  N=ths->N[0];
  n=ths->n[0];

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  f_hat1=ths->f_hat;
  f_hat2=&ths->f_hat[N/2];
  g_hat1=&ths->g_hat[n-N/2];
  g_hat2=ths->g_hat;

  TIC(2)
  nfft_adjoint_1d_B(ths);
  TOC(0)

  TIC_FFTW(1)
  fftw_execute(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv1=ths->c_phi_inv[0];
      c_phi_inv2=&ths->c_phi_inv[0][N/2];
      for(k=0;k<N/2;k++)
	{
	  (*f_hat1++) = (*g_hat1++) * (*c_phi_inv1++);
	  (*f_hat2++) = (*g_hat2++) * (*c_phi_inv2++);
	}
    }
  else
    for(k=0;k<N/2;k++)
      {
	(*f_hat1++) = (*g_hat1++) / (PHI_HUT(k-N/2,0));
	(*f_hat2++) = (*g_hat2++) / (PHI_HUT(k,0));
      }
  TOC(0)
}


/* ############################################################ SPECIFIC VERSIONS FOR d=2 */

static void nfft_2d_init_fg_exp_l(double *fg_exp_l, const int m, const double b)
{
  int l;
  double fg_exp_b0, fg_exp_b1, fg_exp_b2, fg_exp_b0_sq;

  fg_exp_b0 = exp(-1.0/b);
  fg_exp_b0_sq = fg_exp_b0*fg_exp_b0;
  fg_exp_b1 = 1.0;
  fg_exp_b2 = 1.0;
  fg_exp_l[0] = 1.0;
  for(l=1; l <= 2*m+1; l++)
    {
      fg_exp_b2 = fg_exp_b1*fg_exp_b0;
      fg_exp_b1 *= fg_exp_b0_sq;
      fg_exp_l[l] = fg_exp_l[l-1]*fg_exp_b2;
    }
}

static void nfft_trafo_2d_compute(double _Complex *fj, const double _Complex *g,
				  const double *psij_const0, const double *psij_const1,
				  const double *xj0, const double *xj1,
				  const int n0, const int n1, const int m)
{
  int u0,o0,l0,u1,o1,l1;
  const double _Complex *gj;
  const double *psij0,*psij1;

  psij0=psij_const0;
  psij1=psij_const1;

  nfft_uo2(&u0,&o0,*xj0, n0, m);
  nfft_uo2(&u1,&o1,*xj1, n1, m);

  nfft_vpr_double(psij_const0,2*m+2,"psij0");
  nfft_vpr_double(psij_const1,2*m+2,"psij1");

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

static void nfft_adjoint_2d_compute(const double _Complex *fj, double _Complex *g,
				    const double *psij_const0, const double *psij_const1,
				    const double *xj0, const double *xj1,
				    const int n0, const int n1, const int m)
{
  int u0,o0,l0,u1,o1,l1;
  double _Complex *gj;
  const double *psij0,*psij1;

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
  int n0,N0,n1,N1,u,o,j,M,l,m, *psi_index_g,K,ip_s,ip_u;
  double _Complex *fj,*g;
  double *psij, *psij_const, *xj, ip_y, ip_w;

  double *fg_exp_l, fg_psij0, fg_psij1, fg_psij2;

  N0=ths->N[0];
  n0=ths->n[0];
  N1=ths->N[1];
  n1=ths->n[1];
  M=ths->M_total;
  m=ths->m;

  g=ths->g;

  if(ths->nfft_flags & PRE_FULL_PSI)
    {
      psi_index_g=ths->psi_index_g;
      for(j=0, fj=ths->f, psij=ths->psi; j<M; j++, fj++)
        for(l=1, (*fj)=(*psij++) * g[(*psi_index_g++)]; l<(2*m+2)*(2*m+2); l++)
	  (*fj) += (*psij++) * g[(*psi_index_g++)];
      return;
    } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
    {
      for(j=0,fj=ths->f,xj=ths->x; j<M; j++,fj++,xj+=2)
	nfft_trafo_2d_compute(fj, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), xj, xj+1, n0, n1, m);
      return;
    } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
    {
      psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));

      nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
	{
	  fg_psij0 = ths->psi[2*j*2];
	  fg_psij1 = ths->psi[2*j*2+1];
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  fg_psij0 = ths->psi[2*(j*2+1)];
	  fg_psij1 = ths->psi[2*(j*2+1)+1];
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  nfft_trafo_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
    {
      psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));

      nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  fg_psij0 = (PHI(*xj-((double)u)/n0,0));
	  fg_psij1 = exp(2.0*(n0*(*xj) - u)/ths->b[0]);
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_uo(ths,j,&u,&o,1);
	  fg_psij0 = (PHI(*(xj+1)-((double)u)/n1,1));
	  fg_psij1 = exp(2.0*(n1*(*(xj+1)) - u)/ths->b[1]);
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  nfft_trafo_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
    {
      psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
      K=ths->K;
      ip_s=K/(m+1);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  ip_y = fabs(n0*(*(xj+0)) - u)*ip_s;
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[l] = ths->psi[abs(ip_u-l*ip_s)]*(1.0-ip_w) + ths->psi[abs(ip_u-l*ip_s+1)]*(ip_w);

	  printf("idx=%d\tpsi=%e\n",abs(ip_u-(2*m+1)*ip_s+1),ths->psi[abs(ip_u-(2*m+1)*ip_s+1)]);

	  nfft_uo(ths,j,&u,&o,1);
	  ip_y = fabs(n1*(*(xj+1)) - u)*ip_s;
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[2*m+2+l] = ths->psi[(K+1)+abs(ip_u-l*ip_s)]*(1.0-ip_w) + ths->psi[(K+1)+abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_trafo_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
	}
      nfft_free(psij_const);
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
  for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
    {
      nfft_uo(ths,j,&u,&o,0);
      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(*xj-((double)((u+l)))/n0,0));

      nfft_uo(ths,j,&u,&o,1);
      for(l=0;l<=2*m+1;l++)
	psij_const[2*m+2+l]=(PHI(*(xj+1)-((double)((u+l)))/n1,1));

      nfft_trafo_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
    }
  nfft_free(psij_const);
}

static void nfft_adjoint_2d_B(nfft_plan *ths)
{
  int n0,N0,n1,N1,u,o,j,M,l,m, *psi_index_g,K,ip_s,ip_u;
  double _Complex *fj,*g;
  double *psij, *psij_const, *xj ,ip_y, ip_w;

  double *fg_exp_l, fg_psij0, fg_psij1, fg_psij2;

  N0=ths->N[0];
  n0=ths->n[0];
  N1=ths->N[1];
  n1=ths->n[1];
  M=ths->M_total;
  m=ths->m;

  g=ths->g;
  memset(g,0,ths->n_total*sizeof(double _Complex));

  if(ths->nfft_flags & PRE_FULL_PSI)
    {
      psi_index_g=ths->psi_index_g;
      for(j=0, fj=ths->f, psij=ths->psi; j<M; j++, fj++)
	  for(l=0; l<(2*m+2)*(2*m+2); l++)
	      g[(*psi_index_g++)] += (*psij++) * (*fj);
      return;
    } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
    {
      for(j=0,fj=ths->f,xj=ths->x; j<M; j++,fj++,xj+=2)
	nfft_adjoint_2d_compute(fj, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), xj, xj+1, n0, n1, m);
      return;
    } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
    {
      psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));

      nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
	{
	  fg_psij0 = ths->psi[2*j*2];
	  fg_psij1 = ths->psi[2*j*2+1];
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  fg_psij0 = ths->psi[2*(j*2+1)];
	  fg_psij1 = ths->psi[2*(j*2+1)+1];
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  nfft_adjoint_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
    {
      psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));

      nfft_2d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_2d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  fg_psij0 = (PHI(*xj-((double)u)/n0,0));
	  fg_psij1 = exp(2.0*(n0*(*xj) - u)/ths->b[0]);
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_uo(ths,j,&u,&o,1);
	  fg_psij0 = (PHI(*(xj+1)-((double)u)/n1,1));
	  fg_psij1 = exp(2.0*(n1*(*(xj+1)) - u)/ths->b[1]);
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  nfft_adjoint_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
    {
      psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
      K=ths->K;
      ip_s=K/(m+1);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  ip_y = fabs(n0*(*(xj+0)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[l] = ths->psi[abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_uo(ths,j,&u,&o,1);
	  ip_y = fabs(n1*(*(xj+1)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[2*m+2+l] = ths->psi[(K+1)+abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[(K+1)+abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_adjoint_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
	}
      nfft_free(psij_const);
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  psij_const=(double*)nfft_malloc(2*(2*m+2)*sizeof(double));
  for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
    {
      nfft_uo(ths,j,&u,&o,0);
      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(*xj-((double)((u+l)))/n0,0));

      nfft_uo(ths,j,&u,&o,1);
      for(l=0;l<=2*m+1;l++)
	psij_const[2*m+2+l]=(PHI(*(xj+1)-((double)((u+l)))/n1,1));

      nfft_adjoint_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj, xj+1, n0, n1, m);
    }
  nfft_free(psij_const);
}


void nfft_trafo_2d(nfft_plan *ths)
{
  int k0,k1,n0,n1,N0,N1;
  double _Complex *g_hat,*f_hat;
  double *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  double ck01, ck02, ck11, ck12;
  double _Complex *g_hat11,*f_hat11,*g_hat21,*f_hat21,*g_hat12,*f_hat12,*g_hat22,*f_hat22;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];

  f_hat=ths->f_hat;
  g_hat=ths->g_hat;

  TIC(0)
  memset(ths->g_hat,0,ths->n_total*sizeof(double _Complex));
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

      for(k0=0;k0<N0/2;k0++)
	{
	  ck01=(*c_phi_inv01++);
	  ck02=(*c_phi_inv02++);

	  c_phi_inv11=ths->c_phi_inv[1];
	  c_phi_inv12=&ths->c_phi_inv[1][N1/2];

	  g_hat11=g_hat + (n0-N0/2+k0)*n1+n1-N1/2;
	  f_hat11=f_hat + k0*N1;
          g_hat21=g_hat + k0*n1+n1-N1/2;
          f_hat21=f_hat + (N0/2+k0)*N1;
          g_hat12=g_hat + (n0-N0/2+k0)*n1;
          f_hat12=f_hat + k0*N1+N1/2;
	  g_hat22=g_hat + k0*n1;
	  f_hat22=f_hat + (N0/2+k0)*N1+N1/2;
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
    }
  else
    for(k0=0;k0<N0/2;k0++)
      {
	ck01=1./(PHI_HUT(k0-N0/2,0));
	ck02=1./(PHI_HUT(k0,0));
	for(k1=0;k1<N1/2;k1++)
	  {
	    ck11=1./(PHI_HUT(k1-N1/2,1));
	    ck12=1./(PHI_HUT(k1,1));
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
  double _Complex *g_hat,*f_hat;
  double *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  double ck01, ck02, ck11, ck12;
  double _Complex *g_hat11,*f_hat11,*g_hat21,*f_hat21,*g_hat12,*f_hat12,*g_hat22,*f_hat22;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];

  f_hat=ths->f_hat;
  g_hat=ths->g_hat;

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
  
      for(k0=0;k0<N0/2;k0++)
	{
	  ck01=(*c_phi_inv01++);
	  ck02=(*c_phi_inv02++);

	  c_phi_inv11=ths->c_phi_inv[1];
	  c_phi_inv12=&ths->c_phi_inv[1][N1/2];
	  g_hat11=g_hat + (n0-N0/2+k0)*n1+n1-N1/2;
	  f_hat11=f_hat + k0*N1;
          g_hat21=g_hat + k0*n1+n1-N1/2;
          f_hat21=f_hat + (N0/2+k0)*N1;
          g_hat12=g_hat + (n0-N0/2+k0)*n1;
          f_hat12=f_hat + k0*N1+N1/2;
	  g_hat22=g_hat + k0*n1;
	  f_hat22=f_hat + (N0/2+k0)*N1+N1/2;
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
    }
  else
    for(k0=0;k0<N0/2;k0++)
      {
	ck01=1./(PHI_HUT(k0-N0/2,0));
	ck02=1./(PHI_HUT(k0,0));
	for(k1=0;k1<N1/2;k1++)
	  {
	    ck11=1./(PHI_HUT(k1-N1/2,1));
	    ck12=1./(PHI_HUT(k1,1));
	    f_hat[k0*N1+k1]             = g_hat[(n0-N0/2+k0)*n1+n1-N1/2+k1] * ck01 * ck11;
	    f_hat[(N0/2+k0)*N1+k1]      = g_hat[k0*n1+n1-N1/2+k1]           * ck02 * ck11;
	    f_hat[k0*N1+N1/2+k1]        = g_hat[(n0-N0/2+k0)*n1+k1]         * ck01 * ck12;
	    f_hat[(N0/2+k0)*N1+N1/2+k1] = g_hat[k0*n1+k1]                   * ck02 * ck12;
	  }
      }
  TOC(0)
}

/* ############################################################ SPECIFIC VERSIONS FOR d=3 */

static void nfft_3d_init_fg_exp_l(double *fg_exp_l, const int m, const double b)
{
  int l;
  double fg_exp_b0, fg_exp_b1, fg_exp_b2, fg_exp_b0_sq;

  fg_exp_b0 = exp(-1.0/b);
  fg_exp_b0_sq = fg_exp_b0*fg_exp_b0;
  fg_exp_b1 = 1.0;
  fg_exp_b2 = 1.0;
  fg_exp_l[0] = 1.0;
  for(l=1; l <= 2*m+1; l++)
    {
      fg_exp_b2 = fg_exp_b1*fg_exp_b0;
      fg_exp_b1 *= fg_exp_b0_sq;
      fg_exp_l[l] = fg_exp_l[l-1]*fg_exp_b2;
    }
}

static void nfft_trafo_3d_compute(double _Complex *fj, const double _Complex *g,
				  const double *psij_const0, const double *psij_const1, const double *psij_const2,
				  const double *xj0, const double *xj1, const double *xj2,
				  const int n0, const int n1, const int n2, const int m)
{
  int u0,o0,l0,u1,o1,l1,u2,o2,l2;
  const double _Complex *gj;
  const double *psij0,*psij1,*psij2;

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

static void nfft_adjoint_3d_compute(const double _Complex *fj, double _Complex *g,
				    const double *psij_const0, const double *psij_const1, const double *psij_const2,
				    const double *xj0, const double *xj1, const double *xj2,
				    const int n0, const int n1, const int n2, const int m)
{
  int u0,o0,l0,u1,o1,l1,u2,o2,l2;
  double _Complex *gj;
  const double *psij0,*psij1,*psij2;

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
  int n0,N0,n1,N1,n2,N2,u,o,j,M,l,m, *psi_index_g,K,ip_s,ip_u;
  double _Complex *fj,*g;
  double *psij, *psij_const, *xj, ip_y, ip_w;

  double *fg_exp_l, fg_psij0, fg_psij1, fg_psij2;

  N0=ths->N[0];
  n0=ths->n[0];
  N1=ths->N[1];
  n1=ths->n[1];
  N2=ths->N[2];
  n2=ths->n[2];
  M=ths->M_total;
  m=ths->m;

  g=ths->g;

  if(ths->nfft_flags & PRE_FULL_PSI)
    {
      psi_index_g=ths->psi_index_g;
      for(j=0, fj=ths->f, psij=ths->psi; j<M; j++, fj++)
        for(l=1, (*fj)=(*psij++) * g[(*psi_index_g++)]; l<(2*m+2)*(2*m+2)*(2*m+2); l++)
	  (*fj) += (*psij++) * g[(*psi_index_g++)];
      return;
    } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
    {
      for(j=0,fj=ths->f,xj=ths->x; j<M; j++,fj++,xj+=3)
	nfft_trafo_3d_compute(fj, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), xj, xj+1, xj+2, n0, n1, n2, m);
      return;
    } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
    {
      psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));

      nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
	{
	  fg_psij0 = ths->psi[2*j*3];
	  fg_psij1 = ths->psi[2*j*3+1];
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  fg_psij0 = ths->psi[2*(j*3+1)];
	  fg_psij1 = ths->psi[2*(j*3+1)+1];
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  fg_psij0 = ths->psi[2*(j*3+2)];
	  fg_psij1 = ths->psi[2*(j*3+2)+1];
	  fg_psij2 = 1.0;
	  psij_const[2*(2*m+2)] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
	    }

	  nfft_trafo_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
    {
      psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));

      nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  fg_psij0 = (PHI(*xj-((double)u)/n0,0));
	  fg_psij1 = exp(2.0*(n0*(*xj) - u)/ths->b[0]);
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_uo(ths,j,&u,&o,1);
	  fg_psij0 = (PHI(*(xj+1)-((double)u)/n1,1));
	  fg_psij1 = exp(2.0*(n1*(*(xj+1)) - u)/ths->b[1]);
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  nfft_uo(ths,j,&u,&o,2);
	  fg_psij0 = (PHI(*(xj+2)-((double)u)/n2,2));
	  fg_psij1 = exp(2.0*(n2*(*(xj+2)) - u)/ths->b[2]);
	  fg_psij2 = 1.0;
	  psij_const[2*(2*m+2)] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
	    }

	  nfft_trafo_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
    {
      psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
      K=ths->K;
      ip_s=K/(m+1);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  ip_y = fabs(n0*(*(xj+0)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[l] = ths->psi[abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_uo(ths,j,&u,&o,1);
	  ip_y = fabs(n1*(*(xj+1)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[2*m+2+l] = ths->psi[(K+1)+abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[(K+1)+abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_uo(ths,j,&u,&o,2);
	  ip_y = fabs(n2*(*(xj+2)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[2*(K+1)+abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_trafo_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
	}
      nfft_free(psij_const);
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
  for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
    {
      nfft_uo(ths,j,&u,&o,0);
      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(*xj-((double)((u+l)))/n0,0));

      nfft_uo(ths,j,&u,&o,1);
      for(l=0;l<=2*m+1;l++)
	psij_const[2*m+2+l]=(PHI(*(xj+1)-((double)((u+l)))/n1,1));

      nfft_uo(ths,j,&u,&o,2);
      for(l=0;l<=2*m+1;l++)
	psij_const[2*(2*m+2)+l]=(PHI(*(xj+2)-((double)((u+l)))/n2,2));

      nfft_trafo_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
    }
  nfft_free(psij_const);
}

static void nfft_adjoint_3d_B(nfft_plan *ths)
{
  int n0,N0,n1,N1,n2,N2,u,o,j,M,l,m, *psi_index_g,K,ip_s,ip_u;
  double _Complex *fj,*g;
  double *psij, *psij_const, *xj, ip_y, ip_w;

  double *fg_exp_l, fg_psij0, fg_psij1, fg_psij2;

  N0=ths->N[0];
  n0=ths->n[0];
  N1=ths->N[1];
  n1=ths->n[1];
  N2=ths->N[2];
  n2=ths->n[2];
  M=ths->M_total;
  m=ths->m;

  g=ths->g;
  memset(g,0,ths->n_total*sizeof(double _Complex));

  if(ths->nfft_flags & PRE_FULL_PSI)
    {
      psi_index_g=ths->psi_index_g;
      for(j=0, fj=ths->f, psij=ths->psi; j<M; j++, fj++)
        for(l=0; l<(2*m+2)*(2*m+2)*(2*m+2); l++)
	  g[(*psi_index_g++)] += (*psij++) * (*fj);
      return;
    } /* if(PRE_FULL_PSI) */

  if(ths->nfft_flags & PRE_PSI)
    {
      for(j=0,fj=ths->f,xj=ths->x; j<M; j++,fj++,xj+=3)
	nfft_adjoint_3d_compute(fj, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), xj, xj+1, xj+2, n0, n1, n2, m);
      return;
    } /* if(PRE_PSI) */

  if(ths->nfft_flags & PRE_FG_PSI)
    {
      psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));

      nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
	{
	  fg_psij0 = ths->psi[2*j*3];
	  fg_psij1 = ths->psi[2*j*3+1];
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  fg_psij0 = ths->psi[2*(j*3+1)];
	  fg_psij1 = ths->psi[2*(j*3+1)+1];
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  fg_psij0 = ths->psi[2*(j*3+2)];
	  fg_psij1 = ths->psi[2*(j*3+2)+1];
	  fg_psij2 = 1.0;
	  psij_const[2*(2*m+2)] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
	    }

	  nfft_adjoint_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(PRE_FG_PSI) */

  if(ths->nfft_flags & FG_PSI)
    {
      psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
      fg_exp_l=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));

      nfft_3d_init_fg_exp_l(fg_exp_l, m, ths->b[0]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*m+2, m, ths->b[1]);
      nfft_3d_init_fg_exp_l(fg_exp_l+2*(2*m+2), m, ths->b[2]);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  fg_psij0 = (PHI(*xj-((double)u)/n0,0));
	  fg_psij1 = exp(2.0*(n0*(*xj) - u)/ths->b[0]);
	  fg_psij2 = 1.0;
	  psij_const[0] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[l] = fg_psij0*fg_psij2*fg_exp_l[l];
	    }

	  nfft_uo(ths,j,&u,&o,1);
	  fg_psij0 = (PHI(*(xj+1)-((double)u)/n1,1));
	  fg_psij1 = exp(2.0*(n1*(*(xj+1)) - u)/ths->b[1]);
	  fg_psij2 = 1.0;
	  psij_const[2*m+2] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*m+2+l] = fg_psij0*fg_psij2*fg_exp_l[2*m+2+l];
	    }

	  nfft_uo(ths,j,&u,&o,2);
	  fg_psij0 = (PHI(*(xj+2)-((double)u)/n2,2));
	  fg_psij1 = exp(2.0*(n2*(*(xj+2)) - u)/ths->b[2]);
	  fg_psij2 = 1.0;
	  psij_const[2*(2*m+2)] = fg_psij0;
	  for(l=1; l<=2*m+1; l++)
	    {
	      fg_psij2 *= fg_psij1;
	      psij_const[2*(2*m+2)+l] = fg_psij0*fg_psij2*fg_exp_l[2*(2*m+2)+l];
	    }

	  nfft_adjoint_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
	}
      nfft_free(fg_exp_l);
      nfft_free(psij_const);
      return;
    } /* if(FG_PSI) */

  if(ths->nfft_flags & PRE_LIN_PSI)
    {
      psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
      K=ths->K;
      ip_s=K/(m+1);

      for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
	{
	  nfft_uo(ths,j,&u,&o,0);
	  ip_y = fabs(n0*(*(xj+0)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[l] = ths->psi[abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_uo(ths,j,&u,&o,1);
	  ip_y = fabs(n1*(*(xj+1)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[2*m+2+l] = ths->psi[(K+1)+abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[(K+1)+abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_uo(ths,j,&u,&o,2);
	  ip_y = fabs(n2*(*(xj+2)) - u)*((double)K)/(m+1);
	  ip_u = LRINT(floor(ip_y));
	  ip_w = ip_y-ip_u;
	  for(l=0; l < 2*m+2; l++)
	    psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+abs(ip_u-l*ip_s)]*(1.0-ip_w) +
	      ths->psi[2*(K+1)+abs(ip_u-l*ip_s+1)]*(ip_w);

	  nfft_adjoint_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
	}
      nfft_free(psij_const);
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  psij_const=(double*)nfft_malloc(3*(2*m+2)*sizeof(double));
  for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=3)
    {
      nfft_uo(ths,j,&u,&o,0);
      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(*xj-((double)((u+l)))/n0,0));

      nfft_uo(ths,j,&u,&o,1);
      for(l=0;l<=2*m+1;l++)
	psij_const[2*m+2+l]=(PHI(*(xj+1)-((double)((u+l)))/n1,1));

      nfft_uo(ths,j,&u,&o,2);
      for(l=0;l<=2*m+1;l++)
	psij_const[2*(2*m+2)+l]=(PHI(*(xj+2)-((double)((u+l)))/n2,2));

      nfft_adjoint_3d_compute(fj, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xj, xj+1, xj+2, n0, n1, n2, m);
    }
  nfft_free(psij_const);
}

void nfft_trafo_3d(nfft_plan *ths)
{
  int k0,k1,k2,n0,n1,n2,N0,N1,N2;
  double _Complex *g_hat,*f_hat;
  double *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12, *c_phi_inv21, *c_phi_inv22;
  double ck01, ck02, ck11, ck12, ck21, ck22;
  double _Complex *g_hat111,*f_hat111,*g_hat211,*f_hat211,*g_hat121,*f_hat121,*g_hat221,*f_hat221;
  double _Complex *g_hat112,*f_hat112,*g_hat212,*f_hat212,*g_hat122,*f_hat122,*g_hat222,*f_hat222;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  N2=ths->N[2];
  n0=ths->n[0];
  n1=ths->n[1];
  n2=ths->n[2];

  f_hat=ths->f_hat;
  g_hat=ths->g_hat;

  TIC(0)
  memset(ths->g_hat,0,ths->n_total*sizeof(double _Complex));
  if(ths->nfft_flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

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

	      g_hat111=g_hat + ((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2+n2-N2/2;
	      f_hat111=f_hat + (k0*N1+k1)*N2;
	      g_hat211=g_hat + (k0*n1+n1-N1/2+k1)*n2+n2-N2/2;
	      f_hat211=f_hat + ((N0/2+k0)*N1+k1)*N2;
	      g_hat121=g_hat + ((n0-N0/2+k0)*n1+k1)*n2+n2-N2/2;
	      f_hat121=f_hat + (k0*N1+N1/2+k1)*N2;
	      g_hat221=g_hat + (k0*n1+k1)*n2+n2-N2/2;
	      f_hat221=f_hat + ((N0/2+k0)*N1+N1/2+k1)*N2;

	      g_hat112=g_hat + ((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2;
	      f_hat112=f_hat + (k0*N1+k1)*N2+N2/2;
	      g_hat212=g_hat + (k0*n1+n1-N1/2+k1)*n2;
	      f_hat212=f_hat + ((N0/2+k0)*N1+k1)*N2+N2/2;
	      g_hat122=g_hat + ((n0-N0/2+k0)*n1+k1)*n2;
	      f_hat122=f_hat + (k0*N1+N1/2+k1)*N2+N2/2;
	      g_hat222=g_hat + (k0*n1+k1)*n2;
	      f_hat222=f_hat + ((N0/2+k0)*N1+N1/2+k1)*N2+N2/2;

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
    }
  else
    for(k0=0;k0<N0/2;k0++)
      {
	ck01=1./(PHI_HUT(k0-N0/2,0));
	ck02=1./(PHI_HUT(k0,0));
	for(k1=0;k1<N1/2;k1++)
	  {
	    ck11=1./(PHI_HUT(k1-N1/2,1));
	    ck12=1./(PHI_HUT(k1,1));

	    for(k2=0;k2<N2/2;k2++)
	      {
		ck21=1./(PHI_HUT(k2-N2/2,2));
		ck22=1./(PHI_HUT(k2,2));

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
  double _Complex *g_hat,*f_hat;
  double *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12, *c_phi_inv21, *c_phi_inv22;
  double ck01, ck02, ck11, ck12, ck21, ck22;
  double _Complex *g_hat111,*f_hat111,*g_hat211,*f_hat211,*g_hat121,*f_hat121,*g_hat221,*f_hat221;
  double _Complex *g_hat112,*f_hat112,*g_hat212,*f_hat212,*g_hat122,*f_hat122,*g_hat222,*f_hat222;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  N2=ths->N[2];
  n0=ths->n[0];
  n1=ths->n[1];
  n2=ths->n[2];

  f_hat=ths->f_hat;
  g_hat=ths->g_hat;

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

	      g_hat111=g_hat + ((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2+n2-N2/2;
	      f_hat111=f_hat + (k0*N1+k1)*N2;
	      g_hat211=g_hat + (k0*n1+n1-N1/2+k1)*n2+n2-N2/2;
	      f_hat211=f_hat + ((N0/2+k0)*N1+k1)*N2;
	      g_hat121=g_hat + ((n0-N0/2+k0)*n1+k1)*n2+n2-N2/2;
	      f_hat121=f_hat + (k0*N1+N1/2+k1)*N2;
	      g_hat221=g_hat + (k0*n1+k1)*n2+n2-N2/2;
	      f_hat221=f_hat + ((N0/2+k0)*N1+N1/2+k1)*N2;

	      g_hat112=g_hat + ((n0-N0/2+k0)*n1+n1-N1/2+k1)*n2;
	      f_hat112=f_hat + (k0*N1+k1)*N2+N2/2;
	      g_hat212=g_hat + (k0*n1+n1-N1/2+k1)*n2;
	      f_hat212=f_hat + ((N0/2+k0)*N1+k1)*N2+N2/2;
	      g_hat122=g_hat + ((n0-N0/2+k0)*n1+k1)*n2;
	      f_hat122=f_hat + (k0*N1+N1/2+k1)*N2+N2/2;
	      g_hat222=g_hat + (k0*n1+k1)*n2;
	      f_hat222=f_hat + ((N0/2+k0)*N1+N1/2+k1)*N2+N2/2;

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
    }
  else
    for(k0=0;k0<N0/2;k0++)
      {
	ck01=1./(PHI_HUT(k0-N0/2,0));
	ck02=1./(PHI_HUT(k0,0));
	for(k1=0;k1<N1/2;k1++)
	  {
	    ck11=1./(PHI_HUT(k1-N1/2,1));
	    ck12=1./(PHI_HUT(k1,1));

	    for(k2=0;k2<N2/2;k2++)
	      {
		ck21=1./(PHI_HUT(k2-N2/2,2));
		ck22=1./(PHI_HUT(k2,2));

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

  ths->c_phi_inv = (double**) nfft_malloc(ths->d*sizeof(double*));

  for(t=0; t<ths->d; t++)
    {
      ths->c_phi_inv[t]= (double*)nfft_malloc(ths->N[t]*sizeof(double));
      for(ks[t]=0; ks[t]<ths->N[t]; ks[t]++)
	ths->c_phi_inv[t][ks[t]]= 1.0/(PHI_HUT(ks[t]-ths->N[t]/2,t));
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
  double step;                          /**< step size in [0,(m+1)/n]        */

  for (t=0; t<ths->d; t++)
    {
      step=((double)(ths->m+1))/(ths->K*ths->n[t]);
      for(j=0;j<=ths->K;j++)
	{
	  ths->psi[(ths->K+1)*t + j] = PHI(j*step,t);
	} /* for(j) */
    } /* for(t) */
}

static void nfft_precompute_fg_psi(nfft_plan *ths)
{
  int t;                                /**< index over all dimensions       */
  int j;                                /**< index over all nodes            */
  int u, o;                             /**< depends on x_j                  */

  for (t=0; t<ths->d; t++)
    for(j=0;j<ths->M_total;j++)
      {
	nfft_uo(ths,j,&u,&o,t);

        ths->psi[2*(j*ths->d+t)]=
            (PHI((ths->x[j*ths->d+t]-((double)u)/ths->n[t]),t));

        ths->psi[2*(j*ths->d+t)+1]=
            exp(2.0*(ths->n[t]*ths->x[j*ths->d+t] - u) / ths->b[t]);
      } /* for(j) */
  /* for(t) */
} /* nfft_precompute_fg_psi */

void nfft_precompute_psi(nfft_plan *ths)
{
  int t;                                /**< index over all dimensions       */
  int j;                                /**< index over all nodes            */
  int l;                                /**< index u<=l<=o                   */
  int lj;                               /**< index 0<=lj<u+o+1               */
  int u, o;                             /**< depends on x_j                  */

  for (t=0; t<ths->d; t++)
    for(j=0;j<ths->M_total;j++)
      {
	nfft_uo(ths,j,&u,&o,t);

	for(l=u, lj=0; l <= o; l++, lj++)
	  ths->psi[(j*ths->d+t)*(2*ths->m+2)+lj]=
	    (PHI((ths->x[j*ths->d+t]-((double)l)/ths->n[t]),t));
      } /* for(j) */
  /* for(t) */
} /* nfft_precompute_psi */

void nfft_precompute_full_psi(nfft_plan *ths)
{
  int t,t2;                             /**< index over all dimensions       */
  int j;                                /**< index over all nodes            */
  int l_L;                              /**< plain index 0<=l_L<lprod        */
  int l[ths->d];                        /**< multi index u<=l<=o             */
  int lj[ths->d];                       /**< multi index 0<=lj<u+o+1         */
  int ll_plain[ths->d+1];               /**< postfix plain index             */
  int lprod;                            /**< 'bandwidth' of matrix B         */
  int u[ths->d], o[ths->d];             /**< depends on x_j                  */

  double phi_prod[ths->d+1];

  int ix,ix_old;

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

  ths->N_total=nfft_prod_int(ths->N, ths->d);
  ths->n_total=nfft_prod_int(ths->n, ths->d);

  ths->sigma = (double*) nfft_malloc(ths->d*sizeof(double));
  for(t = 0;t < ths->d; t++)
    ths->sigma[t] = ((double)ths->n[t])/ths->N[t];

  WINDOW_HELP_INIT;

  if(ths->nfft_flags & MALLOC_X)
    ths->x = (double*)nfft_malloc(ths->d*ths->M_total*sizeof(double));

  if(ths->nfft_flags & MALLOC_F_HAT)
    ths->f_hat = (double _Complex*)nfft_malloc(ths->N_total*sizeof(double _Complex));

  if(ths->nfft_flags & MALLOC_F)
    ths->f = (double _Complex*)nfft_malloc(ths->M_total*sizeof(double _Complex));

  if(ths->nfft_flags & PRE_PHI_HUT)
    nfft_precompute_phi_hut(ths);

  if(ths->nfft_flags & PRE_LIN_PSI)
  {
      ths->K=(1U<< 10)*(ths->m+1);
      ths->psi = (double*) nfft_malloc((ths->K+1)*ths->d*sizeof(double));
  }

  if(ths->nfft_flags & PRE_FG_PSI)
    ths->psi = (double*) nfft_malloc(ths->M_total*ths->d*2*sizeof(double));

  if(ths->nfft_flags & PRE_PSI)
    ths->psi = (double*) nfft_malloc(ths->M_total*ths->d*
				     (2*ths->m+2)*sizeof(double));

  if(ths->nfft_flags & PRE_FULL_PSI)
  {
      for(t=0,lprod = 1; t<ths->d; t++)
	  lprod *= 2*ths->m+2;

      ths->psi = (double*) nfft_malloc(ths->M_total*lprod*sizeof(double));

      ths->psi_index_f = (int*) nfft_malloc(ths->M_total*sizeof(int));
      ths->psi_index_g = (int*) nfft_malloc(ths->M_total*lprod*sizeof(int));
  }

  if(ths->nfft_flags & FFTW_INIT)
  {
    ths->g1=(double _Complex*)nfft_malloc(ths->n_total*sizeof(double _Complex));

    if(ths->nfft_flags & FFT_OUT_OF_PLACE)
      ths->g2 = (double _Complex*) nfft_malloc(ths->n_total*sizeof(double _Complex));
    else
      ths->g2 = ths->g1;

    ths->my_fftw_plan1 = fftw_plan_dft(ths->d, ths->n, ths->g1, ths->g2, FFTW_FORWARD, ths->fftw_flags);
    ths->my_fftw_plan2 = fftw_plan_dft(ths->d, ths->n, ths->g2, ths->g1,
      FFTW_BACKWARD, ths->fftw_flags);
  }

  ths->mv_trafo = (void (*) (void* ))nfft_trafo;
  ths->mv_adjoint = (void (*) (void* ))nfft_adjoint;
}

void nfft_init(nfft_plan *ths, int d, int *N, int M_total)
{
  int t;                                /**< index over all dimensions       */

  ths->d = d;

  ths->N=(int*) nfft_malloc(d*sizeof(int));
  for(t = 0;t < d; t++)
    ths->N[t] = N[t];

  ths->M_total = M_total;

  ths->n = (int*) nfft_malloc(d*sizeof(int));
  for(t = 0;t < d; t++)
    ths->n[t] = 2*nfft_next_power_of_2(ths->N[t]);

  WINDOW_HELP_ESTIMATE_m;

  ths->nfft_flags = PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                    FFTW_INIT| FFT_OUT_OF_PLACE;
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
  nfft_init(ths,1,N,M_total);
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

void nfft_check(nfft_plan *ths)
{
  int j;

  for(j=0;j<ths->M_total*ths->d;j++)
    if((ths->x[j]<-0.5) || (ths->x[j]>=0.5))
      fprintf(stderr,"nfft_check: ths->x[%d]=%e out of range [-0.5,0.5)\n",
	      j,ths->x[j]);

  for(j=0;j<ths->d;j++)
    {
      if(ths->sigma[j]<=1)
	fprintf(stderr,"nfft_check: oversampling factor too small\n");

      if(ths->N[j]<=ths->m)
	fprintf(stderr,
		"nfft_check: polynomial degree N is smaller than cut-off m\n");

      if(ths->N[j]%2==1)
	fprintf(stderr,"nfft_check: polynomial degree N has to be even\n");
    }
}

void nfft_finalize(nfft_plan *ths)
{
  int t; /* index over dimensions */

  if(ths->nfft_flags & FFTW_INIT)
  {
    fftw_destroy_plan(ths->my_fftw_plan2);
    fftw_destroy_plan(ths->my_fftw_plan1);

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
