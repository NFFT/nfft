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


#define MACRO_nndft_init_result_trafo memset(f,0,ths->M_total*sizeof(double _Complex));
#define MACRO_nndft_init_result_conjugated MACRO_nndft_init_result_trafo
#define MACRO_nndft_init_result_adjoint memset(f_hat,0,ths->N_total*sizeof(double _Complex));
#define MACRO_nndft_init_result_transposed MACRO_nndft_init_result_adjoint

#define MACRO_nndft_sign_trafo      (-2.0*KPI)
#define MACRO_nndft_sign_conjugated (+2.0*KPI)
#define MACRO_nndft_sign_adjoint    (+2.0*KPI)
#define MACRO_nndft_sign_transposed (-2.0*KPI)

#define MACRO_nndft_compute_trafo (*fj) += (*f_hat_k)*cexp(+ _Complex_I*omega);

#define MACRO_nndft_compute_conjugated MACRO_nndft_compute_trafo

#define MACRO_nndft_compute_adjoint (*f_hat_k) += (*fj)*cexp(+ _Complex_I*omega);

#define MACRO_nndft_compute_transposed MACRO_nndft_compute_adjoint

#define MACRO_nndft(which_one)                                                \
void nnfft_ ## which_one ## _direct (nnfft_plan *ths)                                    \
{                                                                             \
  int j;                               /**< index over all nodes (time)     */\
  int t;                               /**< index for dimensions            */\
  int l;                               /**< index over all nodes (fourier)  */\
  double _Complex *f_hat, *f;          /**< dito                            */\
  double _Complex *f_hat_k;            /**< actual Fourier coefficient      */\
  double _Complex *fj;                 /**< actual sample                   */\
  double omega;                        /**< sign times 2*pi*k*x             */\
                                                                              \
  f_hat=ths->f_hat; f=ths->f;                                                 \
                                                                              \
  MACRO_nndft_init_result_ ## which_one                                       \
                                                                              \
  for(j=0, fj=f; j<ths->M_total; j++, fj++)                                   \
  {                                                                           \
    for(l=0, f_hat_k=f_hat; l<ths->N_total; l++, f_hat_k++)                   \
    {                                                                         \
      omega=0.0;                                                              \
      for(t = 0; t<ths->d; t++)                                               \
        omega+=ths->v[l*ths->d+t] * ths->x[j*ths->d+t] * ths->N[t];           \
                                                                              \
      omega*= MACRO_nndft_sign_ ## which_one;                                 \
                                                                              \
      MACRO_nndft_compute_ ## which_one                                       \
                                                                              \
     } /* for(l) */                                                           \
   } /* for(j) */                                                             \
} /* nndft_trafo */                                                           \

MACRO_nndft(trafo)
MACRO_nndft(adjoint)

/** computes 2m+2 indices for the matrix B
 */
static void nnfft_uo(nnfft_plan *ths,int j,int *up,int *op,int act_dim)
{
  double c;
  int u,o;

  c = ths->v[j*ths->d+act_dim] * ths->n[act_dim];

  u = c; o = c;
  if(c < 0)
    u = u-1;
  else
    o = o+1;

  u = u - (ths->m); o = o + (ths->m);

  up[0]=u; op[0]=o;
}

/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$B, B^{\rm T}\f$
 */
#define MACRO_nnfft_B_init_result_A memset(f,0,ths->N_total*sizeof(double _Complex));
#define MACRO_nnfft_B_init_result_T memset(g,0,ths->aN1_total*sizeof(double _Complex));

#define MACRO_nnfft_B_PRE_FULL_PSI_compute_A {                                \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]];                            \
}

#define MACRO_nnfft_B_PRE_FULL_PSI_compute_T {                                \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj);                            \
}

#define MACRO_nnfft_B_compute_A {                                             \
  (*fj) += phi_prod[ths->d] * g[ll_plain[ths->d]];                            \
}

#define MACRO_nnfft_B_compute_T {                                             \
  g[ll_plain[ths->d]] += phi_prod[ths->d] * (*fj);                            \
}

#define MACRO_with_PRE_LIN_PSI (ths->psi[(ths->K+1)*t2+y_u[t2]]*              \
                                (y_u[t2]+1-y[t2]) +                           \
                                ths->psi[(ths->K+1)*t2+y_u[t2]+1]*            \
                                (y[t2]-y_u[t2]))
#define MACRO_with_PRE_PSI     ths->psi[(j*ths->d+t2)*(2*ths->m+2)+lj[t2]]
#define MACRO_without_PRE_PSI  PHI(ths->n[t2], -ths->v[j*ths->d+t2]+                      \
                               ((double)l[t2])/ths->N1[t2], t2)

#define MACRO_init_uo_l_lj_t {                                                \
  for(t = ths->d-1; t>=0; t--)                                                \
    {                                                                         \
      nnfft_uo(ths,j,&u[t],&o[t],t);                                          \
      l[t] = u[t];                                                            \
      lj[t] = 0;                                                              \
    } /* for(t) */                                                            \
  t++;                                                                        \
}

#define MACRO_update_with_PRE_PSI_LIN {                                       \
  for(t2=t; t2<ths->d; t2++)                                                  \
    {                                                                         \
      y[t2] = fabs(((-ths->N1[t2]*ths->v[j*ths->d+t2]+(double)l[t2])          \
          * ((double)ths->K))/(ths->m+1));                                    \
      y_u[t2] = (int)y[t2];                                                   \
    } /* for(t2) */                                                           \
}

#define MACRO_update_phi_prod_ll_plain(which_one) {                           \
  for(t2=t; t2<ths->d; t2++)                                                  \
    {                                                                         \
      phi_prod[t2+1]=phi_prod[t2]* MACRO_ ## which_one;                       \
      ll_plain[t2+1]=ll_plain[t2]*ths->aN1[t2] +                              \
                     (l[t2]+ths->aN1[t2]*3/2)%ths->aN1[t2];                   \
      /* 3/2 because of the (not needed) fftshift and to be in [0 aN1[t2]]?!*/\
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

#define MACRO_nnfft_B(which_one)                                              \
static inline void nnfft_B_ ## which_one (nnfft_plan *ths)                    \
{                                                                             \
  int lprod;                           /**< 'regular bandwidth' of matrix B */\
  int u[ths->d], o[ths->d];            /**< multi band with respect to x_j  */\
  int t, t2;                           /**< index dimensions                */\
  int j;                               /**< index nodes                     */\
  int l_L, ix;                         /**< index one row of B              */\
  int l[ths->d];                       /**< multi index u<=l<=o             */\
  int lj[ths->d];                      /**< multi index 0<=lj<u+o+1         */\
  int ll_plain[ths->d+1];              /**< postfix plain index in g        */\
  double phi_prod[ths->d+1];           /**< postfix product of PHI          */\
  double _Complex *f, *g;              /**< local copy                      */\
  double _Complex *fj;                 /**< local copy                      */\
  double y[ths->d];                                                           \
  int y_u[ths->d];                                                            \
                                                                              \
  f=ths->f_hat; g=ths->F;                                                     \
                                                                              \
  MACRO_nnfft_B_init_result_ ## which_one                                     \
                                                                              \
  if(ths->nnfft_flags & PRE_FULL_PSI)                                         \
    {                                                                         \
      for(ix=0, j=0, fj=f; j<ths->N_total; j++,fj++)                          \
        for(l_L=0; l_L<ths->psi_index_f[j]; l_L++, ix++)                      \
          MACRO_nnfft_B_PRE_FULL_PSI_compute_ ## which_one;                   \
      return;                                                                 \
    }                                                                         \
                                                                              \
  phi_prod[0]=1;                                                              \
  ll_plain[0]=0;                                                              \
                                                                              \
  for(t=0,lprod = 1; t<ths->d; t++)                                           \
    lprod *= (2*ths->m+2);                                                    \
                                                                              \
  if(ths->nnfft_flags & PRE_PSI)                                              \
    {                                                                         \
      for(j=0, fj=f; j<ths->N_total; j++, fj++)                               \
        {                                                                     \
          MACRO_init_uo_l_lj_t;                                               \
                                                                              \
          for(l_L=0; l_L<lprod; l_L++)                                        \
            {                                                                 \
              MACRO_update_phi_prod_ll_plain(with_PRE_PSI);                   \
                                                                              \
              MACRO_nnfft_B_compute_ ## which_one;                            \
                                                                              \
              MACRO_count_uo_l_lj_t;                                          \
            } /* for(l_L) */                                                  \
        } /* for(j) */                                                        \
      return;                                                                 \
    } /* if(PRE_PSI) */                                                       \
                                                                              \
  if(ths->nnfft_flags & PRE_LIN_PSI)                                          \
    {                                                                         \
      for(j=0, fj=f; j<ths->N_total; j++, fj++)                               \
        {                                                                     \
          MACRO_init_uo_l_lj_t;                                               \
                                                                              \
          for(l_L=0; l_L<lprod; l_L++)                                        \
            {                                                                 \
              MACRO_update_with_PRE_PSI_LIN;                                  \
                                                                              \
              MACRO_update_phi_prod_ll_plain(with_PRE_LIN_PSI);               \
                                                                              \
              MACRO_nnfft_B_compute_ ## which_one;                            \
                                                                              \
              MACRO_count_uo_l_lj_t;                                          \
            } /* for(l_L) */                                                  \
        } /* for(j) */                                                        \
      return;                                                                 \
    } /* if(PRE_LIN_PSI) */                                                   \
                                                                              \
  /* no precomputed psi at all */                                             \
  for(j=0, fj=f; j<ths->N_total; j++, fj++)                                   \
    {                                                                         \
                                                                              \
      MACRO_init_uo_l_lj_t;                                                   \
                                                                              \
      for(l_L=0; l_L<lprod; l_L++)                                            \
        {                                                                     \
          MACRO_update_phi_prod_ll_plain(without_PRE_PSI);                    \
                                                                              \
          MACRO_nnfft_B_compute_ ## which_one;                                \
                                                                              \
          MACRO_count_uo_l_lj_t;                                              \
        } /* for(l_L) */                                                      \
    } /* for(j) */                                                            \
} /* nnfft_B */

MACRO_nnfft_B(A)
MACRO_nnfft_B(T)

static inline void nnfft_D (nnfft_plan *ths){
  int j,t;
  double tmp;

  if(ths->nnfft_flags & PRE_PHI_HUT)
  {
      for(j=0; j<ths->M_total; j++)
	  ths->f[j] *= ths->c_phi_inv[j];
  }
  else
  {
      for(j=0; j<ths->M_total; j++)
      {
	  tmp = 1.0;
	  /* multiply with N1, because x was modified */
	  for(t=0; t<ths->d; t++)
	      tmp*= 1.0 /((PHI_HUT(ths->n[t], ths->x[ths->d*j + t]*((double)ths->N[t]),t)) );
	  ths->f[j] *= tmp;
      }
  }
}

/** user routines
 */
void nnfft_trafo(nnfft_plan *ths)
{
  int j,t;

  nnfft_B_T(ths);

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] / ((double)ths->sigma[t]);
    }
  }


  /* allows for external swaps of ths->f */
  ths->direct_plan->f = ths->f;

  nfft_trafo(ths->direct_plan);

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }

  nnfft_D(ths);

} /* nnfft_trafo */

void nnfft_adjoint(nnfft_plan *ths)
{
  int j,t;

  nnfft_D(ths);

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] / ((double)ths->sigma[t]);
    }
  }

  /* allows for external swaps of ths->f */
  ths->direct_plan->f=ths->f;

  nfft_adjoint(ths->direct_plan);

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }

  nnfft_B_A(ths);
} /* nnfft_adjoint */

/** initialisation of direct transform
 */
void nnfft_precompute_phi_hut(nnfft_plan *ths)
{
  int j;                                /**< index over all frequencies       */
  int t;                                /**< index over all dimensions        */
  double tmp;

  ths->c_phi_inv= (double*)nfft_malloc(ths->M_total*sizeof(double));

  for(j=0; j<ths->M_total; j++)
    {
      tmp = 1.0;
      for(t=0; t<ths->d; t++)
        tmp*= 1.0 /(PHI_HUT(ths->n[t],ths->x[ths->d*j + t]*((double)ths->N[t]),t));
      ths->c_phi_inv[j]=tmp;
    }
} /* nnfft_phi_hut */


/** create a lookup table
 */
void nnfft_precompute_lin_psi(nnfft_plan *ths)
{
  int t;                                /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  double step;                          /**< step size in [0,(m+1)/n]         */

  nfft_precompute_lin_psi(ths->direct_plan);

  for (t=0; t<ths->d; t++)
    {
      step=((double)(ths->m+1))/(ths->K*ths->N1[t]);
      for(j=0;j<=ths->K;j++)
        {
          ths->psi[(ths->K+1)*t + j] = PHI(ths->n[t],j*step,t);
        } /* for(j) */
    } /* for(t) */
}

void nnfft_precompute_psi(nnfft_plan *ths)
{
  int t;                                /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int l;                                /**< index u<=l<=o                    */
  int lj;                               /**< index 0<=lj<u+o+1                */
  int u, o;                             /**< depends on v_j                   */

  for (t=0; t<ths->d; t++)
    for(j=0;j<ths->N_total;j++)
      {
        nnfft_uo(ths,j,&u,&o,t);

        for(l=u, lj=0; l <= o; l++, lj++)
          ths->psi[(j*ths->d+t)*(2*ths->m+2)+lj]=
            (PHI(ths->n[t],(-ths->v[j*ths->d+t]+((double)l)/((double)ths->N1[t])),t));
      } /* for(j) */

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] / ((double)ths->sigma[t]);
    }
  }

  nfft_precompute_psi(ths->direct_plan);

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }
  /* for(t) */
} /* nfft_precompute_psi */



/**
 * computes all entries of B explicitly
 */
void nnfft_precompute_full_psi(nnfft_plan *ths)
{
  int t,t2;                             /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int l_L;                              /**< plain index 0<=l_L<lprod         */
  int l[ths->d];                       /**< multi index u<=l<=o              */
  int lj[ths->d];                      /**< multi index 0<=lj<u+o+1          */
  int ll_plain[ths->d+1];              /**< postfix plain index              */
  int lprod;                            /**< 'bandwidth' of matrix B          */
  int u[ths->d], o[ths->d];           /**< depends on x_j                   */

  double phi_prod[ths->d+1];

  int ix,ix_old;

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] / ((double)ths->sigma[t]);
    }
  }

  nnfft_precompute_psi(ths);

  nfft_precompute_full_psi(ths->direct_plan);

  for(j=0;j<ths->M_total;j++) {
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }

  phi_prod[0]=1;
  ll_plain[0]=0;

  for(t=0,lprod = 1; t<ths->d; t++)
    lprod *= 2*ths->m+2;

  for(j=0,ix=0,ix_old=0; j<ths->N_total; j++)
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

void nnfft_precompute_one_psi(nnfft_plan *ths)
{
  if(ths->nnfft_flags & PRE_PSI)
    nnfft_precompute_psi(ths);
  if(ths->nnfft_flags & PRE_FULL_PSI)
    nnfft_precompute_full_psi(ths);
  if(ths->nnfft_flags & PRE_LIN_PSI)
    nnfft_precompute_lin_psi(ths);
  /** precompute phi_hut, the entries of the matrix D */
  if(ths->nnfft_flags & PRE_PHI_HUT)
	  nnfft_precompute_phi_hut(ths);
}

static void nnfft_init_help(nnfft_plan *ths, int m2, unsigned nfft_flags, unsigned fftw_flags)
{
  int t;                                /**< index over all dimensions       */
  int lprod;                            /**< 'bandwidth' of matrix B         */
  int N2[ths->d];

  ths->aN1 = (int*) nfft_malloc(ths->d*sizeof(int));

  ths->a = (double*) nfft_malloc(ths->d*sizeof(double));

  ths->sigma = (double*) nfft_malloc(ths->d*sizeof(double));

  ths->n = ths->N1;

  ths->aN1_total=1;

  for(t = 0; t<ths->d; t++) {
    ths->a[t] = 1.0 + (2.0*((double)ths->m))/((double)ths->N1[t]);
    ths->aN1[t] = ths->a[t] * ((double)ths->N1[t]);
    /* aN1 should be even */
    if(ths->aN1[t]%2 != 0)
      ths->aN1[t] = ths->aN1[t] +1;

    ths->aN1_total*=ths->aN1[t];
    ths->sigma[t] = ((double) ths->N1[t] )/((double) ths->N[t]);;
    
    /* take the same oversampling factor in the inner NFFT */
    N2[t] = ceil(ths->sigma[t]*(ths->aN1[t]));
    
    /* N2 should be even */
    if(N2[t]%2 != 0)
      N2[t] = N2[t] +1;
  }

  WINDOW_HELP_INIT

  if(ths->nnfft_flags & MALLOC_X)
    ths->x = (double*)nfft_malloc(ths->d*ths->M_total*sizeof(double));
  if(ths->nnfft_flags & MALLOC_F){
    ths->f=(double _Complex*)nfft_malloc(ths->M_total*sizeof(double _Complex));
  }
  if(ths->nnfft_flags & MALLOC_V)
    ths->v = (double*)nfft_malloc(ths->d*ths->N_total*sizeof(double));
  if(ths->nnfft_flags & MALLOC_F_HAT)
    ths->f_hat = (double _Complex*)nfft_malloc(ths->N_total*sizeof(double _Complex));

  //BUGFIX SUSE 2
  /** precompute phi_hut, the entries of the matrix D */
//  if(ths->nnfft_flags & PRE_PHI_HUT)
//	  nnfft_precompute_phi_hut(ths);

  if(ths->nnfft_flags & PRE_LIN_PSI)
  {
    ths->K=(1U<< 10)*(ths->m+1);
    ths->psi = (double*) nfft_malloc((ths->K+1)*ths->d*sizeof(double));
  }

  if(ths->nnfft_flags & PRE_PSI){
    ths->psi = (double*)nfft_malloc(ths->N_total*ths->d*(2*ths->m+2)*sizeof(double));
  }

  if(ths->nnfft_flags & PRE_FULL_PSI)
  {
      for(t=0,lprod = 1; t<ths->d; t++)
          lprod *= 2*ths->m+2;

      ths->psi = (double*)nfft_malloc(ths->N_total*lprod*sizeof(double));

      ths->psi_index_f = (int*) nfft_malloc(ths->N_total*sizeof(int));
      ths->psi_index_g = (int*) nfft_malloc(ths->N_total*lprod*sizeof(int));
  }
  ths->direct_plan = (nfft_plan*)nfft_malloc(sizeof(nfft_plan));
  nfft_init_guru(ths->direct_plan, ths->d, ths->aN1, ths->M_total, N2, m2,
		 nfft_flags, fftw_flags);
  ths->direct_plan->x = ths->x;

  ths->direct_plan->f = ths->f;
  ths->F = ths->direct_plan->f_hat;

  ths->mv_trafo = (void (*) (void* ))nnfft_trafo;
  ths->mv_adjoint = (void (*) (void* ))nnfft_adjoint;
}

void nnfft_init_guru(nnfft_plan *ths, int d, int N_total, int M_total, int *N, int *N1,
		     int m, unsigned nnfft_flags)
{
  int t;                             /**< index over all dimensions        */

  unsigned nfft_flags;
  unsigned fftw_flags;

  ths->d= d;
  ths->M_total= M_total;
  ths->N_total= N_total;
  ths->m= m;
  ths->nnfft_flags= nnfft_flags;
  fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;
  nfft_flags= PRE_PHI_HUT| MALLOC_F_HAT| FFTW_INIT|
      ((d == 1) ? FFT_OUT_OF_PLACE : 0U) | NFFT_OMP_BLOCKWISE_ADJOINT;

  if(ths->nnfft_flags & PRE_PSI)
    nfft_flags = nfft_flags | PRE_PSI;

  if(ths->nnfft_flags & PRE_FULL_PSI)
    nfft_flags = nfft_flags | PRE_FULL_PSI;

  if(ths->nnfft_flags & PRE_LIN_PSI)
    nfft_flags = nfft_flags | PRE_LIN_PSI;

  ths->N = (int*) nfft_malloc(ths->d*sizeof(int));
  ths->N1 = (int*) nfft_malloc(ths->d*sizeof(int));

  for(t=0; t<d; t++) {
    ths->N[t] = N[t];
    ths->N1[t] = N1[t];
  }
  nnfft_init_help(ths,m,nfft_flags,fftw_flags);
}

void nnfft_init(nnfft_plan *ths, int d, int N_total, int M_total, int *N)
{
  int t;                            /**< index over all dimensions        */

  unsigned nfft_flags;
  unsigned fftw_flags;
  ths->d = d;
  ths->M_total = M_total;
  ths->N_total = N_total;
  /* m should be greater to get the same accuracy as the nfft */
/* Was soll dieser Ausdruck machen? Es handelt sich um eine Ganzzahl!

  WINDOW_HELP_ESTIMATE_m;
*/
  //BUGFIX SUSE 1
ths->m=WINDOW_HELP_ESTIMATE_m;


  ths->N = (int*) nfft_malloc(ths->d*sizeof(int));
  ths->N1 = (int*) nfft_malloc(ths->d*sizeof(int));

  for(t=0; t<d; t++) {
    ths->N[t] = N[t];
    /* the standard oversampling factor in the nnfft is 1.5 */
    ths->N1[t] = ceil(1.5*ths->N[t]);

    /* N1 should be even */
    if(ths->N1[t]%2 != 0)
      ths->N1[t] = ths->N1[t] +1;
  }

  ths->nnfft_flags=PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F;
  nfft_flags= PRE_PSI| PRE_PHI_HUT| MALLOC_F_HAT| FFTW_INIT|
      ((d == 1) ? FFT_OUT_OF_PLACE : 0U)| NFFT_OMP_BLOCKWISE_ADJOINT;

  fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;
  nnfft_init_help(ths,ths->m,nfft_flags,fftw_flags);
}

void nnfft_init_1d(nnfft_plan *ths,int N1, int M_total)
{
  nnfft_init(ths,1,N1,M_total,&N1);
}

void nnfft_finalize(nnfft_plan *ths)
{
  nfft_finalize(ths->direct_plan);
  nfft_free(ths->direct_plan);

  nfft_free(ths->aN1);
  nfft_free(ths->N);
  nfft_free(ths->N1);

  if(ths->nnfft_flags & PRE_FULL_PSI)
    {
      nfft_free(ths->psi_index_g);
      nfft_free(ths->psi_index_f);
      nfft_free(ths->psi);
    }

  if(ths->nnfft_flags & PRE_PSI)
    nfft_free(ths->psi);

  if(ths->nnfft_flags & PRE_LIN_PSI)
    nfft_free(ths->psi);

  if(ths->nnfft_flags & PRE_PHI_HUT)
    nfft_free(ths->c_phi_inv);

  if(ths->nnfft_flags & MALLOC_F)
    nfft_free(ths->f);

  if(ths->nnfft_flags & MALLOC_F_HAT)
    nfft_free(ths->f_hat);

  if(ths->nnfft_flags & MALLOC_X)
    nfft_free(ths->x);

  if(ths->nnfft_flags & MALLOC_V)
    nfft_free(ths->v);
}
