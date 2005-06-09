#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "options.h"
#include "window_defines.h"

#include "nfft3.h"


#define MACRO_nndft_init_result_trafo memset(f,0,ths->M_total*sizeof(complex));
#define MACRO_nndft_init_result_conjugated MACRO_nndft_init_result_trafo
#define MACRO_nndft_init_result_adjoint memset(f_hat,0,ths->N_total*               \
                                              sizeof(complex));
#define MACRO_nndft_init_result_transposed MACRO_nndft_init_result_adjoint

#define MACRO_nndft_sign_trafo      (+2.0*PI)
#define MACRO_nndft_sign_conjugated (-2.0*PI)
#define MACRO_nndft_sign_adjoint    (-2.0*PI)
#define MACRO_nndft_sign_transposed (+2.0*PI)

#define MACRO_nndft_compute_trafo (*fj) += (*f_hat_k)*cexp(-I*omega);

#define MACRO_nndft_compute_conjugated MACRO_nndft_compute_trafo

#define MACRO_nndft_compute_adjoint (*f_hat_k) += (*fj)*cexp(+I*omega);

#define MACRO_nndft_compute_transposed MACRO_nndft_compute_adjoint

#define MACRO_nndft(which_one)                                                 \
void nndft_ ## which_one (nnfft_plan *ths)                                    \
{                                                                              \
  int j;                                /**< index over all nodes (time)     */\
  int t;                                /**< index for dimensions            */\
  int l;                                /**< index over all nodes (fourier)  */\
  complex *f_hat, *f;                   /**< dito                            */\
  complex *f_hat_k;                     /**< actual Fourier coefficient      */\
  complex *fj;                          /**< actual sample                   */\
  double omega;                         /**< sign times 2*pi*k*x             */\
                                                                               \
  f_hat=ths->f_hat; f=ths->f;                                                \
                                                                               \
  MACRO_nndft_init_result_ ## which_one                                        \
                                                                               \
  for(j=0, fj=f; j<ths->M_total; j++, fj++)                                     \
  {                                                                            \
    for(l=0, f_hat_k=f_hat; l<ths->N_total; l++, f_hat_k++)                    \
    {                                                                          \
      omega=0.0;                                                               \
      for(t = 0; t<ths->d; t++)                                               \
        omega+=ths->v[l*ths->d+t] * ths->x[j*ths->d+t] * ths->N[t];           \
                                                                               \
      omega*= MACRO_nndft_sign_ ## which_one;                                  \
                                                                               \
      MACRO_nndft_compute_ ## which_one                                        \
                                                                               \
     } /* for(l) */                                                            \
   } /* for(j) */                                                              \
} /* nndft_trafo */                                                            \

MACRO_nndft(trafo)
MACRO_nndft(conjugated)
MACRO_nndft(adjoint)
MACRO_nndft(transposed)

/** computes 2m+2 indices for the matrix B
 */
void nnfft_uo(nnfft_plan *ths,int j,int *up,int *op,int act_dim)
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
  /* make the matrix B nonperiodic */

  up[0]=u; op[0]=o;
}

/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$B, B^{\rm T}\f$
 */
#define MACRO_nnfft_B_init_result_A memset(f,0,ths->M_total*sizeof(complex));
#define MACRO_nnfft_B_init_result_T memset(g,0,ths->aN1_total*sizeof(complex));

#define MACRO_nnfft_B_PRE_FULL_PSI_compute_A {                                  \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]];                           \
}

#define MACRO_nnfft_B_PRE_FULL_PSI_compute_T {                                  \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj);                           \
}

#define MACRO_nnfft_B_compute_A {                                               \
  (*fj) += phi_prod[ths->d] * g[ll_plain[ths->d]];                        \
}

#define MACRO_nnfft_B_compute_T {                                               \
  g[ll_plain[ths->d]] += phi_prod[ths->d] * (*fj);                        \
}

  /* Gewicht, d.h., Nachkommaanteil y-y_u im Speicher halten!!! */
#define MACRO_with_PRE_LIN_PSI (ths->psi[(ths->K+1)*t2+y_u[t2]]*             \
                                (y_u[t2]+1-y[t2]) +                            \
                                ths->psi[(ths->K+1)*t2+y_u[t2]+1]*           \
                                (y[t2]-y_u[t2])) 
#define MACRO_with_PRE_PSI     ths->psi[(j*ths->d+t2)*(2*ths->m+2)+lj[t2]]
#define MACRO_without_PRE_PSI  PHI(-ths->v[j*ths->d+t2]+                      \
                               ((double)l[t2])/ths->N1[t2], t2)

#define MACRO_init_uo_l_lj_t {                                                 \
  for(t = ths->d-1; t>=0; t--)                                                \
    {                                                                          \
      nnfft_uo(ths,j,&u[t],&o[t],t);                                           \
      l[t] = u[t];                                                             \
      lj[t] = 0;                                                               \
    } /* for(t) */                                                             \
  t++;                                                                         \
}

#define MACRO_update_with_PRE_PSI_LIN {                                        \
  for(t2=t; t2<ths->d; t2++)                                                  \
    {                                                                          \
      y[t2] = fabs(((-ths->N1[t2]*ths->v[j*ths->d+t2]+(double)l[t2])          \
          * ((double)ths->K))/(ths->m+1));                                   \
      y_u[t2] = (int)y[t2];                                                    \
    } /* for(t2) */                                                            \
}

#define MACRO_update_phi_prod_ll_plain(which_one) {                            \
  for(t2=t; t2<ths->d; t2++)                                                  \
    {                                                                          \
      phi_prod[t2+1]=phi_prod[t2]* MACRO_ ## which_one;                        \
      ll_plain[t2+1]=ll_plain[t2]*ths->aN1[t2] +(l[t2]+ths->aN1[t2]*3/2)%ths->aN1[t2];\
    } /* for(t2) */                                                            \
}

#define MACRO_count_uo_l_lj_t {                                                \
  for(t = ths->d-1; (t>0)&&(l[t]==o[t]); t--)                                 \
    {                                                                          \
      l[t] = u[t];                                                             \
      lj[t] = 0;                                                               \
    } /* for(t) */                                                             \
                                                                               \
  l[t]++;                                                                      \
  lj[t]++;                                                                     \
}

#define MACRO_nnfft_B(which_one)                                                \
inline void nnfft_B_ ## which_one (nnfft_plan *ths)                           \
{                                                                              \
  int lprod;                            /**< 'regular bandwidth' of matrix B */\
  int u[ths->d], o[ths->d];           /**< multi band with respect to x_j  */\
  int t, t2;                            /**< index dimensions                */\
  int j;                                /**< index nodes                     */\
  int l_L, ix;                          /**< index one row of B              */\
  int l[ths->d];                       /**< multi index u<=l<=o             */\
  int lj[ths->d];                      /**< multi index 0<=lj<u+o+1         */\
  int ll_plain[ths->d+1];              /**< postfix plain index in g        */\
  double phi_prod[ths->d+1];           /**< postfix product of PHI          */\
  complex *f, *g;                  /**< local copy                      */\
  complex *fj;                     /**< local copy                      */\
  double y[ths->d];                                                           \
  int y_u[ths->d];                                                            \
                                                                               \
  f=ths->f_hat; g=ths->F;                                                    \
                                                                               \
  MACRO_nnfft_B_init_result_ ## which_one                                                 \
                                                                               \
  if(ths->nnfft_flags & PRE_FULL_PSI)        \
    {                                                                          \
      for(ix=0, j=0, fj=f; j<ths->N_total; j++,fj++)\
        for(l_L=0; l_L<ths->psi_index_f[j]; l_L++, ix++)                      \
          MACRO_nnfft_B_PRE_FULL_PSI_compute_ ## which_one;                                \
      return;                                                                  \
    }                                                                          \
                                                                               \
  phi_prod[0]=1;                                                               \
  ll_plain[0]=0;                                                               \
                                                                               \
  for(t=0,lprod = 1; t<ths->d; t++)                                           \
    lprod *= (2*ths->m+2);                                                    \
                                                                               \
  if(ths->nnfft_flags & PRE_PSI)                                               \
    {                                                                          \
      for(j=0, fj=f; j<ths->N_total; j++, fj++)     \
        {                                                                      \
          MACRO_init_uo_l_lj_t;                                                \
                                                                               \
          for(l_L=0; l_L<lprod; l_L++)                                         \
            {                                                                  \
              MACRO_update_phi_prod_ll_plain(with_PRE_PSI);                    \
                                                                               \
              MACRO_nnfft_B_compute_ ## which_one;                              \
                                                                               \
              MACRO_count_uo_l_lj_t;                                           \
            } /* for(l_L) */                                                   \
        } /* for(j) */                                                         \
      return;                                                                  \
    } /* if(PRE_PSI) */                                                        \
                                                                               \
  if(ths->nnfft_flags & PRE_LIN_PSI)                                           \
    {                                                                          \
      for(j=0, fj=f; j<ths->N_total; j++, fj++)     \
        {                                                            \
          MACRO_init_uo_l_lj_t;                                                \
                                                                               \
          for(l_L=0; l_L<lprod; l_L++)                                         \
            {                                                                  \
              MACRO_update_with_PRE_PSI_LIN;                                   \
                                                                               \
              MACRO_update_phi_prod_ll_plain(with_PRE_LIN_PSI);                \
                                                                               \
              MACRO_nnfft_B_compute_ ## which_one;                              \
                                                                               \
              MACRO_count_uo_l_lj_t;                                           \
            } /* for(l_L) */                                                   \
        } /* for(j) */                                                         \
      return;                                                                  \
    } /* if(PRE_LIN_PSI) */                                                    \
                                                                               \
  /* no precomputed psi at all */                                              \
  for(j=0, fj=f; j<ths->N_total; j++, fj++)         \
    {                  \
                                                \
      MACRO_init_uo_l_lj_t;                                                    \
                                                                               \
      for(l_L=0; l_L<lprod; l_L++)                                             \
        {                                                                      \
          MACRO_update_phi_prod_ll_plain(without_PRE_PSI);                     \
                                                                              \
          MACRO_nnfft_B_compute_ ## which_one;                                  \
                                                                               \
          MACRO_count_uo_l_lj_t;                                               \
        } /* for(l_L) */                                                     \
    } /* for(j) */                                                            \
} /* nnfft_B */               

MACRO_nnfft_B(A)
MACRO_nnfft_B(T)

inline void nnfft_D (nnfft_plan *ths){
  int j,t;
  double tmp;
  
  if(ths->nnfft_flags & PRE_PHI_HUT) {
    for(j=0; j<ths->M_total; j++)
      ths->f[j] *= ths->c_phi_inv[j];
  } else {
    for(j=0; j<ths->M_total; j++)
    {
      tmp = 1.0;
      /* multiply with N1, because x was modified */
      for(t=0; t<ths->d; t++)
        tmp*= 1.0 /((PHI_HUT(ths->x[ths->d*j + t]*((double)ths->N[t]),t)) );
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
  
  ths->direct_plan->f = ths->f;
  nfft_trafo(ths->direct_plan);
  
  for(j=0;j<ths->M_total;j++) {  
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }
  
  nnfft_D(ths);
} /* nnfft_trafo */

void nnfft_conjugated(nnfft_plan *ths)
{
  int j,t;
  
  nnfft_B_T(ths);
  
  for(j=0;j<ths->M_total;j++) {  
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] / ((double)ths->sigma[t]);
    }
  }
  
  ths->direct_plan->f = ths->f;
  nfft_conjugated(ths->direct_plan);

  for(j=0;j<ths->M_total;j++) {  
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }  
  
  nnfft_D(ths);
} /* nnfft_conjugated */

void nnfft_adjoint(nnfft_plan *ths)
{
  int j,t;
  
  nnfft_D(ths);
  
  for(j=0;j<ths->M_total;j++) {  
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] / ((double)ths->sigma[t]);
    }
  }
  
  ths->direct_plan->f = ths->f;
  nfft_adjoint(ths->direct_plan);
  
  for(j=0;j<ths->M_total;j++) {  
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }  
  
  nnfft_B_A(ths);
} /* nnfft_adjoint */

void nnfft_transposed(nnfft_plan *ths)
{
  int j,t;
  
  nnfft_D(ths);
  
  for(j=0;j<ths->M_total;j++) {  
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] / ((double)ths->sigma[t]);
    }
  }
  
  ths->direct_plan->f = ths->f;
  nfft_transposed(ths->direct_plan);
  
  for(j=0;j<ths->M_total;j++) {  
    for(t=0;t<ths->d;t++) {
      ths->x[j*ths->d+t]= ths->x[j*ths->d+t] * ((double)ths->sigma[t]);
    }
  }  
  
  nnfft_B_A(ths);
} /* nnfft_transposed */


/** initialisation of direct transform 
 */
void nnfft_precompute_phi_hut(nnfft_plan *ths)
{
  int j;                                /**< index over all frequencies       */
  int t;                                /**< index over all dimensions        */
  double tmp;

  ths->c_phi_inv= (double*)fftw_malloc(ths->M_total*sizeof(double));
  
  for(j=0; j<ths->M_total; j++)
    {
      tmp = 1.0;
      for(t=0; t<ths->d; t++)
        tmp*= 1.0 /(PHI_HUT(ths->x[ths->d*j + t]*((double)ths->N[t]),t));
      ths->c_phi_inv[j]=tmp;
    }
} /* nnfft_phi_hut */


/** create a lookup table, but NOT for each node
 *  good idea K=2^xx
 *  TODO: estimate K, call from init
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
          ths->psi[(ths->K+1)*t + j] = PHI(j*step,t);
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
            (PHI((-ths->v[j*ths->d+t]+((double)l)/((double)ths->N1[t])),t));
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

/** computes all entries of B explicitly */
void nnfft_precompute_full_psi(nnfft_plan *ths, double eps)
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

  int *index_g, *index_f;
  double *new_psi;
  int ix,ix_old,size_psi;
  
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

  size_psi=ths->N_total;
  index_f = (int*) malloc(ths->N_total*sizeof(int));
  index_g = (int*) malloc(size_psi*sizeof(int));
  new_psi = (double*) malloc(size_psi*sizeof(double));
  
  for(t=0,lprod = 1; t<ths->d; t++)
    {
      lprod *= 2*ths->m+2;
      eps=eps*(PHI(0,t));
    }
  
  for(ix=0,ix_old=0,j=0; j<ths->N_total; j++)
    {
      MACRO_init_uo_l_lj_t;
      
      for(l_L=0; l_L<lprod; l_L++)
        {
          MACRO_update_phi_prod_ll_plain(with_PRE_PSI);
          
          if(phi_prod[ths->d]>eps)
            {
              index_g[ix]=ll_plain[ths->d];
              new_psi[ix]=phi_prod[ths->d];
              ix++;
              if(ix==size_psi)
                {
                  size_psi+=ths->N_total;
                  index_g=(int*)realloc(index_g,size_psi*sizeof(int));
                  new_psi=(double*)realloc(new_psi,size_psi*sizeof(double));
                }
            }
          MACRO_count_uo_l_lj_t;
        } /* for(l_L) */
      
      index_f[j]=ix-ix_old;
      ix_old=ix;
    } /* for(j) */
  
  free(ths->psi);
  size_psi=ix;
  ths->size_psi=size_psi;
  index_g=(int*)realloc(index_g,size_psi*sizeof(int));
  new_psi=(double*)realloc(new_psi,size_psi*sizeof(double));
  
  ths->psi         =new_psi;
  ths->psi_index_g =index_g; 
  ths->psi_index_f=index_f;
}

void nnfft_init_help(nnfft_plan *ths, int m2, int *N2, unsigned nfft_flags, unsigned fftw_flags)
{
  int t;                                /**< index over all dimensions        */
  
  ths->aN1 = (int*) fftw_malloc(ths->d*sizeof(int));
  
  ths->a = (double*) fftw_malloc(ths->d*sizeof(double));
  
  ths->sigma = (double*) fftw_malloc(ths->d*sizeof(double));
  
  ths->n = ths->N1;
  
  ths->aN1_total=1;
  
  for(t = 0; t<ths->d; t++) {
    ths->a[t] = 1.0 + (2.0*((double)ths->m))/((double)ths->N1[t]);
    ths->aN1[t] = ths->a[t] * ((double)ths->N1[t]);
    if(ths->aN1[t]%2 != 0)
      ths->aN1[t] = ths->aN1[t] +1;
      
    ths->aN1_total*=ths->aN1[t];
    ths->sigma[t] = ((double) ths->N1[t] )/((double) ths->N[t]);;
  }
  
  WINDOW_HELP_INIT
  
  if(ths->nnfft_flags & MALLOC_X)
    ths->x = (double*)fftw_malloc(ths->d*ths->M_total*
                                        sizeof(double));
                                        
  if(ths->nnfft_flags & MALLOC_V)
    ths->v = (double*)fftw_malloc(ths->d*ths->N_total*
                                        sizeof(double));

  if(ths->nnfft_flags & MALLOC_F_HAT)
    ths->f_hat = (complex*)fftw_malloc(ths->N_total*
                                                  sizeof(complex));
  if(ths->nnfft_flags & MALLOC_F)
    ths->f=(complex*)fftw_malloc(ths->M_total*sizeof(complex));
    
  if(ths->nnfft_flags & PRE_LIN_PSI)
  {
    ths->K=100000; /* estimate is badly needed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    ths->psi = (double*) fftw_malloc((ths->K+1)*ths->d*sizeof(double));
  }
    
  /* NO FFTW_MALLOC HERE */
  if(ths->nnfft_flags & PRE_PSI || ths->nnfft_flags & PRE_FULL_PSI)
    ths->psi = (double*) malloc(ths->N_total*ths->d*
                                           (2*ths->m+2)*sizeof(double));
                                           
  ths->direct_plan = (nfft_plan*) malloc(sizeof(nfft_plan));
    
  nfft_init_guru(ths->direct_plan, ths->d, ths->aN1, ths->M_total, N2, m2, 
	               nfft_flags, fftw_flags);
                        
  ths->direct_plan->x = ths->x;
  ths->direct_plan->f = ths->f;
  ths->F = ths->direct_plan->f_hat;
}

void nnfft_init_guru(nnfft_plan *ths, int d, int N_total, int M_total, int *N, int *N1,
                        int m, unsigned nnfft_flags)
{
  int t;                             /**< index over all dimensions        */
  
  int N2[d];
  unsigned nfft_flags;
  unsigned fftw_flags;
  
  ths->d= d;
  ths->M_total= M_total;  
  ths->N_total= N_total;
  ths->m= m;
  ths->nnfft_flags= nnfft_flags;
  fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;
  nfft_flags= PRE_PHI_HUT| MALLOC_F_HAT| FFTW_INIT| FFT_OUT_OF_PLACE;
  
  if(ths->nnfft_flags & PRE_PSI)
    nfft_flags = nfft_flags | PRE_PSI;
    
  if(ths->nnfft_flags & PRE_FULL_PSI)
    nfft_flags = nfft_flags | PRE_FULL_PSI;
    
  if(ths->nnfft_flags & PRE_LIN_PSI)
    nfft_flags = nfft_flags | PRE_LIN_PSI;
  
  ths->N = (int*) fftw_malloc(ths->d*sizeof(int));
  ths->N1 = (int*) fftw_malloc(ths->d*sizeof(int));
  
  for(t=0; t<d; t++) {
    ths->N[t] = N[t];
    ths->N1[t] = N1[t];    
    N2[t] = 2*next_power_of_2((int)((double)ths->N1[t])* 
           (1.0+2.0*((double)ths->m)/((double)ths->N1[t])));
  }
  nnfft_init_help(ths,m,N2,nfft_flags,fftw_flags);  
}

void nnfft_init(nnfft_plan *ths, int d, int N_total, int M_total, int *N)
{
  int t;                            /**< index over all dimensions        */

  int N2[d];
  unsigned nfft_flags;
  unsigned fftw_flags;

  ths->d = d;
  ths->M_total = M_total;
  ths->N_total = N_total;

  WINDOW_HELP_ESTIMATE_m;
  
  ths->N = (int*) fftw_malloc(ths->d*sizeof(int));  
  ths->N1 = (int*) fftw_malloc(ths->d*sizeof(int));
  
  for(t=0; t<d; t++) {
    ths->N[t] = N[t];
    ths->N1[t] = 2*next_power_of_2(ths->N[t]);
    N2[t] = 2*next_power_of_2((int)((double)ths->N1[t])* 
            (1.0+2.0*((double)ths->m)/((double)ths->N1[t])));
  }
  ths->nnfft_flags=PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F;
  nfft_flags= PRE_PSI| PRE_PHI_HUT| MALLOC_F_HAT| FFTW_INIT| FFT_OUT_OF_PLACE;
  
  fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;
  
  nnfft_init_help(ths,ths->m,N2,nfft_flags,fftw_flags);    
}

void nnfft_finalize(nnfft_plan *ths)
{ 
  nfft_finalize(ths->direct_plan);
  
  free(ths->direct_plan);

  free(ths->aN1);
  free(ths->N);
  free(ths->N1);
  
  /* NO FFTW_FREE HERE */
  if((ths->nnfft_flags & PRE_PSI) || (ths->nnfft_flags & PRE_LIN_PSI))
    {
      if(ths->nnfft_flags & PRE_FULL_PSI)
        {
          free(ths->psi_index_g);
          free(ths->psi_index_f);
        }

      free(ths->psi);
    }
      
  if(ths->nnfft_flags & PRE_PHI_HUT)
    fftw_free(ths->c_phi_inv);
  
  if(ths->nnfft_flags & MALLOC_F)
    fftw_free(ths->f);

  if(ths->nnfft_flags & MALLOC_F_HAT)
    fftw_free(ths->f_hat);

  if(ths->nnfft_flags & MALLOC_X)
    fftw_free(ths->x);
    
  if(ths->nnfft_flags & MALLOC_V)
    fftw_free(ths->v);
}
