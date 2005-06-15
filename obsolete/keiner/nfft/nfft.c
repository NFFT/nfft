/**
 * Library.
 * Includes simple and fast computation of the NDFT (direct problem) as well
 * as (iterative) solution to the inverse problem.
 * authors: D. Potts, S. Kunis (c) 2002,2003
 */

#include "nfft.h"
#include "window_defines.h"

/** direct computation of non equispaced fourier transforms
 *  ndft_trafo, ndft_conjugated, ndft_adjoint, ndft_transposed
 *  require O(M N^d) arithemtical operations
 *
 * direct computation of the ndft_trafo and ndft_conjugated, formula (1.1)
 * ndft_trafo:
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 * ndft_conjugated:
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 *
 * direct computation of the ndft_adjoint and ndft_transposed, formula (1.2)
 * ndft_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 * ndft_transposed:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(-2(pi) k x[j])
 */

/** macros and small sub routines for the direct transforms
 */
#define MACRO_ndft_init_result_trafo memset(f,0,this->M*sizeof(fftw_complex));
#define MACRO_ndft_init_result_conjugated MACRO_ndft_init_result_trafo
#define MACRO_ndft_init_result_adjoint memset(f_hat,0,this->N_L*               \
                                              sizeof(fftw_complex));
#define MACRO_ndft_init_result_transposed MACRO_ndft_init_result_adjoint

#define MACRO_ndft_sign_trafo      +2*PI*this->x[j*this->d+t]
#define MACRO_ndft_sign_conjugated -2*PI*this->x[j*this->d+t]
#define MACRO_ndft_sign_adjoint    +2*PI*this->x[j*this->d+t]
#define MACRO_ndft_sign_transposed -2*PI*this->x[j*this->d+t]

#define MACRO_init_k_N_Omega_x(which_one) {                                    \
for(t=0; t<this->d; t++)                                                       \
  {                                                                            \
    k[t]=-this->N[t]/2;                                                        \
    x[t]= MACRO_ndft_sign_ ## which_one;                                       \
    Omega[t+1]=k[t]*x[t]+Omega[t];                                             \
  }                                                                            \
omega=Omega[this->d];                                                          \
}                                                                              \

#define MACRO_count_k_N_Omega {                                                \
for(t = this->d-1; (t >= 1) && (k[t] == this->N[t]/2-1); t--)                  \
  k[t]-= this->N[t]-1;                                                         \
                                                                               \
k[t]++;                                                                        \
                                                                               \
for(t2 = t; t2<this->d; t2++)                                                  \
  Omega[t2+1]=k[t2]*x[t2]+Omega[t2];                                           \
                                                                               \
omega=Omega[this->d];                                                          \
}

#define MACRO_ndft_compute_trafo {                                             \
  (*fj0)+=(*f_hat_k0)*cos_omega+(*f_hat_k1)*sin_omega;                         \
  (*fj1)+=(*f_hat_k1)*cos_omega-(*f_hat_k0)*sin_omega;                         \
}
#define MACRO_ndft_compute_conjugated MACRO_ndft_compute_trafo

#define MACRO_ndft_compute_adjoint {                                           \
  (*f_hat_k0)+= (*fj0)*cos_omega-(*fj1)*sin_omega;                             \
  (*f_hat_k1)+= (*fj0)*sin_omega+(*fj1)*cos_omega;                             \
}
#define MACRO_ndft_compute_transposed MACRO_ndft_compute_adjoint

#define MACRO_ndft(which_one)                                                  \
void ndft_ ## which_one (nfft_plan *this)                                      \
{                                                                              \
  int j;                                /**< index over all nodes            */\
  int t,t2;                             /**< index for dimensions            */\
  int k_L;                              /**< plain index for summation       */\
  fftw_complex *f_hat, *f;              /**< dito                            */\
  double *f_hat_k0, *f_hat_k1;          /**< actual Fourier coefficient      */\
  double *fj0,*fj1;                     /**< actual sample                   */\
  double x[this->d];                    /**< actual node x[d*j+t]            */\
  int k[this->d];                       /**< multi index for summation       */\
  double omega, Omega[this->d+1];       /**< sign times 2*pi*k*x             */\
  double sin_omega,cos_omega;           /**< just for better reading         */\
                                                                               \
  f_hat=this->f_hat; f=this->f;                                                \
                                                                               \
  MACRO_ndft_init_result_ ## which_one                                         \
                                                                               \
  if(this->d==1)                                                               \
    {                                                                          \
  /* univariate case (due to performance) */                                   \
      t=0;                                                                     \
      for(j=0, fj0= &f[0][0], fj1= &f[0][1];                                   \
          j<this->M; j++, fj0+=2, fj1+=2)                                      \
        {                                                                      \
	  for(k_L=0, f_hat_k0= &f_hat[0][0], f_hat_k1= &f_hat[0][1];           \
              k_L<this->N_L; k_L++, f_hat_k0+=2, f_hat_k1+=2)                  \
	    {                                                                  \
	      omega=(k_L-this->N_L/2)* MACRO_ndft_sign_ ## which_one;          \
	      cos_omega= cos(omega);                                           \
	      sin_omega= sin(omega);                                           \
              MACRO_ndft_compute_ ## which_one;                                \
	    }                                                                  \
        }                                                                      \
    }                                                                          \
  else                                                                         \
    {                                                                          \
  /* multivariate case */                                                      \
      Omega[0]=0;                                                              \
      for(j=0, fj0=&f[0][0], fj1=&f[0][1];                                     \
          j<this->M; j++, fj0+=2, fj1+=2)                                      \
        {                                                                      \
          MACRO_init_k_N_Omega_x(which_one);                                   \
          for(k_L=0, f_hat_k0= &f_hat[0][0], f_hat_k1= &f_hat[0][1];           \
              k_L<this->N_L; k_L++, f_hat_k0+=2, f_hat_k1+=2)                  \
	    {                                                                  \
              cos_omega= cos(omega);                                           \
              sin_omega= sin(omega);                                           \
                                                                               \
              MACRO_ndft_compute_ ## which_one;                                \
                                                                               \
	      MACRO_count_k_N_Omega;                                           \
	    } /* for(k_L) */                                                   \
        } /* for(j) */                                                         \
    } /* else */                                                               \
} /* ndft_trafo */


/** user routines
 */
MACRO_ndft(trafo)
MACRO_ndft(conjugated)
MACRO_ndft(adjoint)
MACRO_ndft(transposed)


/** fast computation of non equispaced fourier transforms
 *  require O(N^d log(N) + M) arithemtical operations
 *
 * fast computation of the nfft_trafo and nfft_conjugated, formula (1.1)
 * nfft_trafo:
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 * nfft_conjugated:
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint and nfft_transposed, formula (1.2)
 * nfft_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 * nfft_transposed:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(-2(pi) k x[j])
 */

/** macros and small sub routines for the fast transforms
 */

/** computes 2m+2 indices for the matrix B
 */
void nfft_uo(nfft_plan *this,int j,int *up,int *op,int act_dim)
{
  double c;
  int u,o;
	
  c = this->x[j*this->d+act_dim] * this->n[act_dim];
  u = c; o = c;
  if(c < 0)                  
    u = u-1;                  
  else
    o = o+1;
  
  u = u - (this->m); o = o + (this->m);

  up[0]=u; op[0]=o;
}

#define MACRO_nfft_D_compute_A {                                               \
 g_hat[k_plain[this->d]][0]= f_hat[ks_plain[this->d]][0]*c_phi_inv_k[this->d]; \
 g_hat[k_plain[this->d]][1]= f_hat[ks_plain[this->d]][1]*c_phi_inv_k[this->d]; \
}

#define MACRO_nfft_D_compute_T {                                               \
 f_hat[ks_plain[this->d]][0]= g_hat[k_plain[this->d]][0]*c_phi_inv_k[this->d]; \
 f_hat[ks_plain[this->d]][1]= g_hat[k_plain[this->d]][1]*c_phi_inv_k[this->d]; \
}

#define MACRO_nfft_D_init_result_A  memset(g_hat,0,this->n_L*                  \
                                          sizeof(fftw_complex));
#define MACRO_nfft_D_init_result_T memset(f_hat,0,this->N_L*                   \
                                          sizeof(fftw_complex));

#define MACRO_with_PRE_PHI_HUT * this->c_phi_inv[t2][ks[t2]];
#define MACRO_without_PRE_PHI_HUT / (PHI_HUT(ks[t2]-this->N[t2]/2,t2));

#define MACRO_init_k_ks {                                                      \
  for(t = this->d-1; t>=0; t--)                                                \
    {                                                                          \
      kp[t]= 0;                                                                \
      k[t] = 0;                                                                \
      ks[t] = this->N[t]/2;                                                    \
    }                                                                          \
  t++;                                                                         \
}

#define MACRO_update_c_phi_inv_k(which_one) {                                  \
  for(t2=0; t2<this->d; t2++)                                                  \
    {                                                                          \
      c_phi_inv_k[t2+1]= c_phi_inv_k[t2] MACRO_ ##which_one;                   \
      ks_plain[t2+1]= ks_plain[t2]*this->N[t2]+ks[t2];                         \
      k_plain[t2+1]= k_plain[t2]*this->n[t2]+k[t2];                            \
    }                                                                          \
}

#define MACRO_count_k_ks {                                                     \
  for(t=this->d-1; (t>0)&& (kp[t]==this->N[t]-1); t--)                         \
    {                                                                          \
      kp[t]= 0;                                                                \
      k[t]= 0;                                                                 \
      ks[t]= this->N[t]/2;                                                     \
    }                                                                          \
                                                                               \
  kp[t]++; k[t]++; ks[t]++;                                                    \
  if(kp[t]==this->N[t]/2)                                                      \
    {                                                                          \
      k[t]= this->n[t]-this->N[t]/2;                                           \
      ks[t]= 0;                                                                \
    }                                                                          \
}                                                                              \


/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$D, D^T\f$
 */
#define MACRO_nfft_D(which_one)                                                \
inline void nfft_D_ ## which_one (nfft_plan *this)                             \
{                                                                              \
  int t, t2, t_temp;                            /**< index dimensions                */\
  int k_L;                              /**< plain index                     */\
  int kp[this->d];                      /**< multi index (simple)            */\
  int k[this->d];                       /**< multi index in g_hat            */\
  int ks[this->d];                      /**< multi index in f_hat, c_phi_inv */\
  double c_phi_inv_k[this->d+1];        /**< postfix product of PHI_HUT      */\
  int k_plain[this->d+1];               /**< postfix plain index             */\
  int ks_plain[this->d+1];              /**< postfix plain index             */\
  fftw_complex *f_hat, *g_hat;          /**< local copy                      */\
                                                                               \
  f_hat=this->f_hat; g_hat=this->g_hat;                                        \
  MACRO_nfft_D_init_result_ ## which_one;                                      \
                                                                               \
  c_phi_inv_k[0]=1;                                                            \
  k_plain[0]=0;                                                                \
  ks_plain[0]=0;                                                               \
                                                                               \
  if(this->nfft_flags & PRE_PHI_HUT)                                           \
    {                                                                          \
      MACRO_init_k_ks;                                                         \
                                     \
      t_temp = t;																															\
      for(k_L=0; k_L<this->N_L; k_L++)                                         \
	{                                                                      \
          MACRO_update_c_phi_inv_k(with_PRE_PHI_HUT);                          \
                                                                               \
	  MACRO_nfft_D_compute_ ## which_one;                                  \
	                                                                       \
	  MACRO_count_k_ks;                                                    \
	} /* for(k_L) */                                                       \
    } /* if(PRE_PHI_HUT) */                                                    \
  else                                                                         \
    {                                                                          \
      MACRO_init_k_ks;                                                         \
                                                                               \
      t_temp = t;																															\
      for(k_L=0; k_L<this->N_L; k_L++)                                         \
	{                                                                      \
          MACRO_update_c_phi_inv_k(without_PRE_PHI_HUT);                       \
                                                                               \
	  MACRO_nfft_D_compute_ ## which_one;                                  \
	                                                                       \
	  MACRO_count_k_ks;                                                    \
	} /* for(k_L) */                                                       \
    } /* else(PRE_PHI_HUT) */                                                  \
} /* nfft_D */


MACRO_nfft_D(A)
MACRO_nfft_D(T)

/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$B, B^{\rm T}\f$
 */ 
#define MACRO_nfft_B_init_result_A  memset(f,0,this->M*sizeof(fftw_complex));
#define MACRO_nfft_B_init_result_T memset(g,0,this->n_L*sizeof(fftw_complex));

#define MACRO_nfft_B_PRE_FULL_PSI_compute_A {                                  \
  (*fj0)+=this->psi[ix]*g[this->psi_index_g[ix]][0];                           \
  (*fj1)+=this->psi[ix]*g[this->psi_index_g[ix]][1];                           \
}

#define MACRO_nfft_B_PRE_FULL_PSI_compute_T {                                  \
  g[this->psi_index_g[ix]][0]+=this->psi[ix]*(*fj0);                           \
  g[this->psi_index_g[ix]][1]+=this->psi[ix]*(*fj1);                           \
}

#define MACRO_nfft_B_compute_A {                                               \
  (*fj0) += phi_prod[this->d]* g[ll_plain[this->d]][0];                        \
  (*fj1) += phi_prod[this->d]* g[ll_plain[this->d]][1];                        \
}

#define MACRO_nfft_B_compute_T {                                               \
  g[ll_plain[this->d]][0] += phi_prod[this->d]* (*fj0);                        \
  g[ll_plain[this->d]][1] += phi_prod[this->d]* (*fj1);                        \
}

#define MACRO_with_PRE_PSI     this->psi[(j*this->d+t2)*(2*this->m+2)+lj[t2]]
#define MACRO_without_PRE_PSI  PHI(this->x[j*this->d+t2]-                      \
                               ((double)l[t2])/this->n[t2], t2)

#define MACRO_init_uo_l_lj_t {                                                 \
  for(t = this->d-1; t>=0; t--)                                                \
    {                                                                          \
      nfft_uo(this,j,&u[t],&o[t],t);                                           \
      l[t] = u[t];                                                             \
      lj[t] = 0;                                                               \
    } /* for(t) */                                                             \
  t++;                                                                         \
}

#define MACRO_update_phi_prod_ll_plain(which_one) {                            \
  for(t2=t; t2<this->d; t2++)                                                  \
    {                                                                          \
      phi_prod[t2+1]=phi_prod[t2]* MACRO_ ## which_one;                        \
      ll_plain[t2+1]=ll_plain[t2]*this->n[t2] +(l[t2]+this->n[t2])%this->n[t2];\
    } /* for(t2) */                                                            \
}

#define MACRO_count_uo_l_lj_t {                                                \
  for(t = this->d-1; (t>0)&&(l[t]==o[t]); t--)                                 \
    {                                                                          \
      l[t] = u[t];                                                             \
      lj[t] = 0;                                                               \
    } /* for(t) */                                                             \
                                                                               \
  l[t]++;                                                                      \
  lj[t]++;                                                                     \
}

#define MACRO_nfft_B(which_one)                                                \
inline void nfft_B_ ## which_one (nfft_plan *this)                             \
{                                                                              \
  int lprod;                            /**< 'regular bandwidth' of matrix B */\
  int u[this->d], o[this->d];           /**< multi band with respect to x_j  */\
  int t, t2;                            /**< index dimensions                */\
  int j;                                /**< index nodes                     */\
  int l_L, ix;                          /**< index one row of B              */\
  int l[this->d];                       /**< multi index u<=l<=o             */\
  int lj[this->d];                      /**< multi index 0<=lj<u+o+1         */\
  int ll_plain[this->d+1];              /**< postfix plain index in g        */\
  double phi_prod[this->d+1];           /**< postfix product of PHI          */\
  fftw_complex *f, *g;                  /**< local copy                      */\
  double *fj0,*fj1;                     /**< local copy                      */\
                                                                               \
  f=this->f; g=this->g;                                                        \
                                                                               \
  MACRO_nfft_B_init_result_ ## which_one;                                      \
                                                                               \
  if((this->nfft_flags & PRE_PSI)&&(this->nfft_flags & PRE_FULL_PSI))          \
    {                                                                          \
      for(ix=0, j=0, fj0=&f[0][0], fj1=&f[0][1]; j<this->M; j++,fj0+=2,fj1+=2) \
        for(l_L=0; l_L<this->psi_index_f[j]; l_L++, ix++)                      \
	  MACRO_nfft_B_PRE_FULL_PSI_compute_ ## which_one;                     \
      return;                                                                  \
    }                                                                          \
                                                                               \
  phi_prod[0]=1;                                                               \
  ll_plain[0]=0;                                                               \
                                                                               \
  for(t=0,lprod = 1; t<this->d; t++)                                           \
    lprod *= (2*this->m+2);                                                    \
                                                                               \
  if(this->nfft_flags & PRE_PSI)                                               \
    {                                                                          \
      for(j=0, fj0=&f[0][0], fj1=&f[0][1]; j<this->M; j++, fj0+=2, fj1+=2)     \
	{                                                                      \
          MACRO_init_uo_l_lj_t;                                                \
                                                                               \
	  for(l_L=0; l_L<lprod; l_L++)                                         \
	    {                                                                  \
              MACRO_update_phi_prod_ll_plain(with_PRE_PSI);                    \
                                                                               \
	      MACRO_nfft_B_compute_ ## which_one;                              \
		                                                               \
	      MACRO_count_uo_l_lj_t;                                           \
            } /* for(l_L) */                                                   \
	} /* for(j) */                                                         \
      return;                                                                  \
    } /* if(PRE_PSI) */                                                        \
                                                                               \
  /* no precomputed psi at all */                                              \
  for(j=0, fj0=&f[0][0], fj1=&f[0][1]; j<this->M; j++, fj0+=2, fj1+=2)         \
    {                                                                          \
      MACRO_init_uo_l_lj_t;                                                    \
	                                                                       \
      for(l_L=0; l_L<lprod; l_L++)                                             \
     	{                                                                      \
          MACRO_update_phi_prod_ll_plain(without_PRE_PSI);                     \
                                                                               \
          MACRO_nfft_B_compute_ ## which_one;                                  \
		                                                               \
          MACRO_count_uo_l_lj_t;                                               \
	} /* for(l_L) */                                                       \
    } /* for(j) */                                                             \
} /* nfft_B */                                                                 \

MACRO_nfft_B(A)
MACRO_nfft_B(T)


/** user routines
 */
void nfft_trafo(nfft_plan *this)
{
  /* use this->my_fftw_plan1 */
  this->g_hat=this->g1;
  this->g=this->g2;
 
  /** form \f$ \hat g_k = \frac{\hat{f}_k}{c_k\left(\phi\right)} \text{ for }
   *  k \in I_N \f$
   */
  //printf("before T1\n");
	//fflush(stdout);
  T1;
  //printf("before D_A\n");
	//fflush(stdout);	
  nfft_D_A(this);
  //printf("before T_2\n");
	//fflush(stdout);
  T2(1);

  /** compute by d-variate discrete Fourier transform
   *  \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   *  \text{ for } l \in I_n \f$
   */
  T1;
  fftw_execute(this->my_fftw_plan1);
  T2(2);

  /** set \f$ f_j = \sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
   *  \text{ for } j=0,\hdots,M-1 \f$
   */
  T1;
  nfft_B_A(this);
  T2(3);
} /* nfft_trafo */

void nfft_conjugated(nfft_plan *this)
{
  /* use this->my_fftw_plan2 */
  this->g_hat=this->g2;
  this->g=this->g1;
 
  /** form \f$ \hat g_k = \frac{\hat{f}_k}{c_k\left(\phi\right)} \text{ for }
   *  k \in I_N \f$
   */  
  T1;
  nfft_D_A(this);
  T2(1);

  /** compute by d-variate discrete Fourier transform
   *  \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{+2\pi {\rm i} \frac{kl}{n}}
   *  \text{ for } l \in I_n \f$
   */
  T1;
  fftw_execute(this->my_fftw_plan2);
  T2(2);

  /** set \f$ f_j = \sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
   *  \text{ for } j=0,\hdots,M-1 \f$
   */
  T1;
  nfft_B_A(this);
  T2(3);
} /* nfft_conjugated */

void nfft_adjoint(nfft_plan *this)
{
  /* use this->my_fftw_plan2 */
  this->g_hat=this->g1;
  this->g=this->g2;
  
  /** set \f$ g_l = \sum_{j=0}^{M-1} f_j \psi\left(x_j-\frac{l}{n}\right)
   *  \text{ for } l \in I_n,m(x_j) \f$
   */
  T1;
  nfft_B_T(this);
  T2(1);
 
  /** compute by d-variate discrete Fourier transform
   *  \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{+2\pi {\rm i} \frac{kl}{n}}
   *  \text{ for }  k \in I_N\f$
   */
  T1;
  fftw_execute(this->my_fftw_plan2);
  T2(2);
  
  /** form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
   *  k \in I_N \f$
   */
  T1;
  nfft_D_T(this);
  T2(3);
} /* nfft_adjoint */

void nfft_transposed(nfft_plan *this)
{
  /* use this->my_fftw_plan1 */
  this->g_hat=this->g2;
  this->g=this->g1;

  /** set \f$ g_l = \sum_{j=0}^{M-1} f_j \psi\left(x_j-\frac{l}{n}\right)
   *  \text{ for } l \in I_n,m(x_j) \f$
   */
  T1;
  nfft_B_T(this);
  T2(1);
 
  /** compute by d-variate discrete Fourier transform
   *  \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   *  \text{ for }  k \in I_N\f$
   */ 
  T1;
  fftw_execute(this->my_fftw_plan1);
  T2(2);
  
  /** form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
   *  k \in I_N \f$
   */
  T1;
  nfft_D_T(this);
  T2(3);
} /* nfft_transposed */


/** initialisation of direct transform 
 */
void nfft_precompute_phi_hut(nfft_plan *this)
{
  int ks[this->d];                      /**< index over all frequencies       */
  int t;                                /**< index over all dimensions        */

  this->c_phi_inv = (double**) fftw_malloc(this->d*sizeof(double*));

  for(t=0; t<this->d; t++)
    {
      this->c_phi_inv[t]= (double*)fftw_malloc(this->N[t]*sizeof(double));
      for(ks[t]=0; ks[t]<this->N[t]; ks[t]++)  
	this->c_phi_inv[t][ks[t]]= 1.0/(PHI_HUT(ks[t]-this->N[t]/2,t));
    }
} /* nfft_phi_hut */

void nfft_precompute_psi(nfft_plan *this)
{
  int t;                                /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int l;                                /**< index u<=l<=o                    */
  int lj;                               /**< index 0<=lj<u+o+1                */
  int u, o;                             /**< depends on x_j                   */
  
  for (t=0; t<this->d; t++)
    for(j=0;j<this->M;j++)
      {
	nfft_uo(this,j,&u,&o,t);
	
	for(l=u, lj=0; l <= o; l++, lj++)
	  this->psi[(j*this->d+t)*(2*this->m+2)+lj]=
	    (PHI((this->x[j*this->d+t]-((double)l)/this->n[t]),t));
      } /* for(j) */
  /* for(t) */
} /* nfft_precompute_psi */

/** more memory usage, a bit faster */
void nfft_full_psi(nfft_plan *this, double eps)
{
  int t,t2;                             /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int l_L;                              /**< plain index 0<=l_L<lprod         */
  int l[this->d];                       /**< multi index u<=l<=o              */
  int lj[this->d];                      /**< multi index 0<=lj<u+o+1          */
  int ll_plain[this->d+1];              /**< postfix plain index              */
  int lprod;                            /**< 'bandwidth' of matrix B          */
  int u[this->d], o[this->d];           /**< depends on x_j                   */
  
  double phi_prod[this->d+1];

  int *index_g, *index_f;
  double *new_psi;
  int ix,ix_old,size_psi;

  phi_prod[0]=1;
  ll_plain[0]=0;

  if(this->nfft_flags & PRE_PSI)
    {
      size_psi=this->M;
      index_f = (int*) malloc(this->M*sizeof(int));
      index_g = (int*) malloc(size_psi*sizeof(int));
      new_psi = (double*) malloc(size_psi*sizeof(double));
      
      for(t=0,lprod = 1; t<this->d; t++)
	{
	  lprod *= 2*this->m+2;
	  eps=eps*(PHI(0,t));
	}
      
      for(ix=0,ix_old=0,j=0; j<this->M; j++)
	{
	  MACRO_init_uo_l_lj_t;
	  
	  for(l_L=0; l_L<lprod; l_L++)
            {
	      MACRO_update_phi_prod_ll_plain(with_PRE_PSI);

              if(phi_prod[this->d]>eps)
		{
		  index_g[ix]=ll_plain[this->d];
		  new_psi[ix]=phi_prod[this->d];
		  ix++;
		  if(ix==size_psi)
		    {
		      size_psi+=this->M;
		      index_g=(int*)realloc(index_g,size_psi*sizeof(int));
		      new_psi=(double*)realloc(new_psi,size_psi*sizeof(double));
		    }
		}
	      MACRO_count_uo_l_lj_t;
            } /* for(l_L) */

	  
	  index_f[j]=ix-ix_old;
	  ix_old=ix;
	} /* for(j) */

      free(this->psi);
      size_psi=ix;
      this->size_psi=size_psi;
      index_g=(int*)realloc(index_g,size_psi*sizeof(int));
      new_psi=(double*)realloc(new_psi,size_psi*sizeof(double));

      this->psi         =new_psi;
      this->psi_index_g =index_g; 
      this->psi_index_f=index_f;
      
      //printf("size_psi / (lprod*M)=%e\n",size_psi / ((double)lprod*this->M));
    } /* if(PRE_PSI) */
}


void nfft_init_help(nfft_plan *this)
{
  int t;                                /**< index over all dimensions        */

  this->N_L=nfft_prod_int(this->N, this->d);
  this->n_L=nfft_prod_int(this->n, this->d);

  this->sigma = (double*) fftw_malloc(this->d*sizeof(double));
  for(t = 0;t < this->d; t++)
    this->sigma[t] = ((double)this->n[t])/this->N[t];
  
  WINDOW_HELP_INIT;

  if(this->nfft_flags & MALLOC_X)
    this->x = (double*)fftw_malloc(this->d*this->M*
					sizeof(double));

  if(this->nfft_flags & MALLOC_F_HAT)
    this->f_hat = (fftw_complex*)fftw_malloc(this->N_L*
						  sizeof(fftw_complex));
  if(this->nfft_flags & MALLOC_F)
    this->f=(fftw_complex*)fftw_malloc(this->M*sizeof(fftw_complex));

  if(this->nfft_flags & PRE_PHI_HUT)
    nfft_precompute_phi_hut(this);

  /* NO FFTW_MALLOC HERE */
  if(this->nfft_flags & PRE_PSI)
    this->psi = (double*) malloc(this->M*this->d*
					   (2*this->m+2)*sizeof(double));
  
  this->g1=(fftw_complex*)fftw_malloc(this->n_L*sizeof(fftw_complex));

  if(this->nfft_flags & FFT_OUT_OF_PLACE)
    this->g2 = (fftw_complex*) fftw_malloc(this->n_L*
						sizeof(fftw_complex));
  else
    this->g2 = this->g1;
  
  this->my_fftw_plan1 = 
    fftw_plan_dft(this->d, this->n, this->g1, this->g2,
		  FFTW_FORWARD, this->fftw_flags);
  this->my_fftw_plan2 = 
    fftw_plan_dft(this->d, this->n, this->g2, this->g1,
		  FFTW_BACKWARD, this->fftw_flags);
}

void nfft_init(nfft_plan *this, int d, int *N, int M)
{
  int t;                                /**< index over all dimensions        */

  this->d = d;

  this->N=(int*) fftw_malloc(d*sizeof(int));
  for(t = 0;t < d; t++)
    this->N[t] = N[t];

  this->M = M;

  this->n = (int*) fftw_malloc(d*sizeof(int));
  for(t = 0;t < d; t++)
    this->n[t] = 2*next_power_of_2(this->N[t]);

  WINDOW_HELP_ESTIMATE_m;

  this->nfft_flags= PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F| FFT_OUT_OF_PLACE;
  this->fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;

  nfft_init_help(this);    
}

void nfft_init_specific(nfft_plan *this, int d, int *N, int M, int *n,
			int m, unsigned nfft_flags, unsigned fftw_flags)
{
  int t;                                /**< index over all dimensions        */

  this->d =d;
  this->N= (int*) fftw_malloc(this->d*sizeof(int));
  for(t=0; t<d; t++)
    this->N[t]= N[t];
  this->M= M;
  this->n= (int*) fftw_malloc(this->d*sizeof(int));
  for(t=0; t<d; t++)
    this->n[t]= n[t];
  this->m= m;
  this->nfft_flags= nfft_flags;
  this->fftw_flags= fftw_flags;

  nfft_init_help(this);  
}

void nfft_init_1d(nfft_plan *this, int N1, int M)
{
  int N[1];

  N[0]=N1;
  nfft_init(this,1,N,M);
}

void nfft_init_2d(nfft_plan *this, int N1, int N2, int M)
{
  int N[2];

  N[0]=N1;
  N[1]=N2;
  nfft_init(this,2,N,M);
}

void nfft_init_3d(nfft_plan *this, int N1, int N2, int N3, int M)
{
  int N[3];

  N[0]=N1;
  N[1]=N2;
  N[2]=N3;
  nfft_init(this,3,N,M);
}

void nfft_finalize(nfft_plan *this)
{
  int t; /* index over dimensions */

  fftw_destroy_plan(this->my_fftw_plan2);
  fftw_destroy_plan(this->my_fftw_plan1);

  if(this->nfft_flags & FFT_OUT_OF_PLACE)
    fftw_free(this->g2);

  fftw_free(this->g1);

  /* NO FFTW_FREE HERE */
  if(this->nfft_flags & PRE_PSI)
    {
      if(this->nfft_flags & PRE_FULL_PSI)
	{
	  free(this->psi_index_g);
	  free(this->psi_index_f);
	}

      free(this->psi);
    }
      
  if(this->nfft_flags & PRE_PHI_HUT)
    {
      for(t=0; t<this->d; t++)
        fftw_free(this->c_phi_inv[t]);
      fftw_free(this->c_phi_inv);
    }

  if(this->nfft_flags & MALLOC_F)
    fftw_free(this->f);

  if(this->nfft_flags & MALLOC_F_HAT)
    fftw_free(this->f_hat);

  if(this->nfft_flags & MALLOC_X)
  fftw_free(this->x);
 
  WINDOW_HELP_FINALIZE;
 
  fftw_free(this->sigma);
  fftw_free(this->n);
  fftw_free(this->N);
}


/** inverse nfft 
 * -----------------------------------------------------------------------------
 * -----------------------------------------------------------------------------
 */
void infft_init_specific(infft_plan *this_iplan, nfft_plan *direct_plan,
			 int infft_flags)
{
  this_iplan->direct_plan = direct_plan;
  this_iplan->infft_flags = infft_flags;
  
  this_iplan->given_f = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));

  this_iplan->r_iter = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));

  this_iplan->f_hat_iter = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));

  this_iplan->p_hat_iter = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));

  if(this_iplan->infft_flags & LANDWEBER)
    {
      this_iplan->z_hat_iter = this_iplan->p_hat_iter;
    }

  if(this_iplan->infft_flags & STEEPEST_DESCENT)
    {
      this_iplan->z_hat_iter = this_iplan->p_hat_iter;

      this_iplan->v_iter = (fftw_complex*)
	fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));
    }

  if(this_iplan->infft_flags & CGNR_E)
    {
      this_iplan->z_hat_iter = (fftw_complex*)
	fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));

      this_iplan->v_iter = (fftw_complex*)
	fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));
    }

  if(this_iplan->infft_flags & CGNE_R)
    {
      this_iplan->z_hat_iter = this_iplan->p_hat_iter;
    }

  if(this_iplan->infft_flags & ITERATE_2nd)
    this_iplan->f_hat_iter_2nd = (fftw_complex*) 
      fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));

  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->w = 
      (double*) fftw_malloc(this_iplan->direct_plan->M*sizeof(double));
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->w_hat = 
      (double*) fftw_malloc(this_iplan->direct_plan->N_L*sizeof(double));
}

void infft_init(infft_plan *this_iplan, nfft_plan *direct_plan)
{
  infft_init_specific(this_iplan, direct_plan, CGNR_E);
}

void infft_before_loop_help(infft_plan *this_iplan)
{
  /** step 2
   *  overwrites this_iplan->direct_plan->f_hat 
   *  overwrites this_iplan->r_iter
   */
  copyc(this_iplan->direct_plan->f_hat, this_iplan->f_hat_iter,
	this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->r_iter, this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->r_iter, this_iplan->direct_plan->f);

  updatec_axpy(this_iplan->r_iter, -1.0, this_iplan->given_f,
	       this_iplan->direct_plan->M);

  if((!(this_iplan->infft_flags & LANDWEBER)) ||
     (this_iplan->infft_flags & NORMS_FOR_LANDWEBER))
    {
      if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
	this_iplan->dot_r_iter =
	  dotproductc_w(this_iplan->r_iter, this_iplan->w,
			this_iplan->direct_plan->M);
      else
	this_iplan->dot_r_iter =
	  dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);
    }
  
  /** step 3
   *  overwrites this_iplan->direct_plan->f
   *  overwrites this_iplan->z_hat_iter resp. this_iplan->z_hat_iter
   */ 

  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w, this_iplan->r_iter,
	    this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
	  this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat);
  
  if((!(this_iplan->infft_flags & LANDWEBER)) ||
     (this_iplan->infft_flags & NORMS_FOR_LANDWEBER))
    {
      if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
	this_iplan->dot_z_hat_iter = 
	  dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat,
			this_iplan->direct_plan->N_L);
      else
	this_iplan->dot_z_hat_iter =
	  dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
    }

  if(this_iplan->infft_flags & CGNE_R)
    this_iplan->dot_p_hat_iter = this_iplan->dot_z_hat_iter;

} /* void infft_before_loop_help */

void infft_before_loop(infft_plan *this_iplan)
{
  infft_before_loop_help(this_iplan);

  if(this_iplan->infft_flags & CGNR_E)
    {
      /** step 4-6
       *  overwrites this_iplan->f_hat_iter_2nd
       */
      if(this_iplan->infft_flags & ITERATE_2nd)
	copyc(this_iplan->f_hat_iter_2nd, this_iplan->f_hat_iter,
	      this_iplan->direct_plan->N_L);
     
      /** step 7
       *  overwrites this_iplan->p_hat_iter
       */
      copyc(this_iplan->p_hat_iter, this_iplan->z_hat_iter,
	    this_iplan->direct_plan->N_L);
    }

  if(this_iplan->infft_flags & CGNE_R)
    {
      /** step 4-7
       *  overwrites this_iplan->f_hat_iter_2nd
       */
      if(this_iplan->infft_flags & ITERATE_2nd)
	{
	  this_iplan->gamma_iter=1.0;
	  copyc(this_iplan->f_hat_iter_2nd, this_iplan->f_hat_iter,
		this_iplan->direct_plan->N_L);
	}
    }
} /* void infft_before_loop */

void infft_loop_one_step_landweber(infft_plan *this_iplan)
{
  /** step 5
   *  updates this_iplan->f_hat_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
		  this_iplan->w_hat, this_iplan->z_hat_iter,
		  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter,
		this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  
  /** step 6
   *  original residual, not the updated residual,
   *  overwrites this_iplan->r_iter
   *  overwrites this_iplan->direct_plan->f_hat 
   */
  copyc(this_iplan->direct_plan->f_hat, this_iplan->f_hat_iter,
	this_iplan->direct_plan->N_L);

  SWAPC(this_iplan->r_iter,this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->r_iter,this_iplan->direct_plan->f);
  
  updatec_axpy(this_iplan->r_iter, -1.0, this_iplan->given_f, 
	       this_iplan->direct_plan->M);

  if(this_iplan->infft_flags & NORMS_FOR_LANDWEBER)
    {
      if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
	this_iplan->dot_r_iter = dotproductc_w(this_iplan->r_iter,this_iplan->w,
					       this_iplan->direct_plan->M);
      else
	this_iplan->dot_r_iter =
	  dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);
    }

  /** step 7
   *  overwrites this_iplan->direct_plan->f 
   *  overwrites this_iplan->z_hat_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w,
	    this_iplan->r_iter, this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
	  this_iplan->direct_plan->M);
    
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);

  if(this_iplan->infft_flags & NORMS_FOR_LANDWEBER)
    {
      if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
	this_iplan->dot_z_hat_iter = 
	  dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat,
			this_iplan->direct_plan->N_L);
      else
	this_iplan->dot_z_hat_iter =
	  dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
    }
} /* void infft_loop_one_step_landweber */

void infft_loop_one_step_steepest_descent(infft_plan *this_iplan)
{
  /** step 5
   *  overwrites this_iplan->direct_plan->f_hat 
   *  overwrites this_iplan->v_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    copyc_w(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
	    this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  else
    copyc(this_iplan->direct_plan->f_hat, this_iplan->z_hat_iter,
	  this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_v_iter = dotproductc_w(this_iplan->v_iter, this_iplan->w,
					   this_iplan->direct_plan->M);
  else
    this_iplan->dot_v_iter =
      dotproductc(this_iplan->v_iter, this_iplan->direct_plan->M);
  
  /** step 6
   */
  this_iplan->alpha_iter =
    this_iplan->dot_z_hat_iter / this_iplan->dot_v_iter;

  /** step 7
   *  updates this_iplan->f_hat_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
		  this_iplan->w_hat,this_iplan->z_hat_iter,
		  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter,
		this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);

  /** step 8
   *  updates this_iplan->r_iter
   */
  updatec_xpay(this_iplan->r_iter, -this_iplan->alpha_iter, this_iplan->v_iter,
	      this_iplan->direct_plan->M);

  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter = dotproductc_w(this_iplan->r_iter, this_iplan->w,
					   this_iplan->direct_plan->M);
  else
    this_iplan->dot_r_iter =
      dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);

  /** step 9
   *  overwrites this_iplan->direct_plan->f
   *  overwrites this_iplan->z_hat_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w, this_iplan->r_iter,
	    this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
	  this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->dot_z_hat_iter =
      dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat, 
		    this_iplan->direct_plan->N_L);
  else
    this_iplan->dot_z_hat_iter =
      dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
} /* void infft_loop_one_step_steepest_descent */

void infft_loop_one_step_cgnr_e(infft_plan *this_iplan)
{
  /** step 9
   *  overwrites this_iplan->direct_plan->f_hat 
   *  overwrites this_iplan->v_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    copyc_w(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
	    this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);
  else
    copyc(this_iplan->direct_plan->f_hat, this_iplan->p_hat_iter,
	  this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_v_iter = dotproductc_w(this_iplan->v_iter, this_iplan->w,
					   this_iplan->direct_plan->M);
  else
    this_iplan->dot_v_iter =
      dotproductc(this_iplan->v_iter, this_iplan->direct_plan->M);
  
  /** step 10
   */
  this_iplan->alpha_iter =
    this_iplan->dot_z_hat_iter / this_iplan->dot_v_iter;

  /** step 11
   *  updates this_iplan->f_hat_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
		  this_iplan->w_hat, this_iplan->p_hat_iter,
		  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter, 
		this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);

  /** step 12-15
   *  updates this_iplan->f_hat_iter_2nd
   */
  if(this_iplan->infft_flags & ITERATE_2nd)
    {
      this_iplan->alpha_iter_2nd =
	this_iplan->dot_r_iter / this_iplan->dot_z_hat_iter;
      
      if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
	updatec_xpawy(this_iplan->f_hat_iter_2nd, this_iplan->alpha_iter_2nd,
		      this_iplan->w_hat, this_iplan->z_hat_iter,
		      this_iplan->direct_plan->N_L);
      else
	updatec_xpay(this_iplan->f_hat_iter_2nd, this_iplan->alpha_iter_2nd,
		    this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
    }
  
  /** step 16
   *  updates this_iplan->r_iter
   */
  updatec_xpay(this_iplan->r_iter, -this_iplan->alpha_iter, this_iplan->v_iter,
	      this_iplan->direct_plan->M);
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter = dotproductc_w(this_iplan->r_iter, this_iplan->w,
					   this_iplan->direct_plan->M);
  else
    this_iplan->dot_r_iter =
      dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);

  /** step 17
   *  overwrites this_iplan->direct_plan->f
   *  overwrites this_iplan->direct_plan->r_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w, 
	    this_iplan->r_iter, this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
	  this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);

  this_iplan->dot_z_hat_iter_old = this_iplan->dot_z_hat_iter;
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->dot_z_hat_iter =
      dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat, 
		    this_iplan->direct_plan->N_L);
  else
    this_iplan->dot_z_hat_iter =
      dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  
  /** step 18
   */
  this_iplan->beta_iter =
    this_iplan->dot_z_hat_iter / this_iplan->dot_z_hat_iter_old;
  
  /** step 19
   *  updates this_iplan->p_hat_iter
   */
  updatec_axpy(this_iplan->p_hat_iter, this_iplan->beta_iter, 
	       this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
} /* void infft_loop_one_step_cgnr_e */

void infft_loop_one_step_cgne_r(infft_plan *this_iplan)
{
  /** step 9
   */
  this_iplan->alpha_iter =
    this_iplan->dot_r_iter / this_iplan->dot_p_hat_iter;
  
  /** step 10
   *  updates this_iplan->f_hat_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
		  this_iplan->w_hat, this_iplan->p_hat_iter,
		  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter,
		this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);
  
  /** step 11
   *  overwrites this_iplan->direct_plan->f_hat 
   *  overwrites this_iplan->direct_plan->f
   *  updates this_iplan->r_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    copyc_w(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
	    this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);
  else
    copyc(this_iplan->direct_plan->f_hat, this_iplan->p_hat_iter,
	  this_iplan->direct_plan->N_L);
  
  nfft_trafo(this_iplan->direct_plan);

  updatec_xpay(this_iplan->r_iter, -this_iplan->alpha_iter,
	      this_iplan->direct_plan->f, this_iplan->direct_plan->M);
  
  this_iplan->dot_r_iter_old = this_iplan->dot_r_iter;
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter =
      dotproductc_w(this_iplan->r_iter, this_iplan->w,
		    this_iplan->direct_plan->M);
  else
    this_iplan->dot_r_iter =
      dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);
  
  /** step 12
   */
  this_iplan->beta_iter =
    this_iplan->dot_r_iter / this_iplan->dot_r_iter_old;
  
  /** step 13-16
   *  updates this_iplan->f_hat_iter_2nd
   */
  if(this_iplan->infft_flags & ITERATE_2nd)
    {
      this_iplan->gamma_iter_old = this_iplan->gamma_iter;
      this_iplan->gamma_iter =
	this_iplan->beta_iter * this_iplan->gamma_iter_old + 1;
      
      updatec_axpby(this_iplan->f_hat_iter_2nd, 
		    this_iplan->beta_iter * this_iplan->gamma_iter_old / 
		    this_iplan->gamma_iter, this_iplan->f_hat_iter, 
		    1.0 / this_iplan->gamma_iter, this_iplan->direct_plan->N_L);
    }
  
  /** step 16
   *  overwrites this_iplan->direct_plan->f
   *  overwrites this_iplan->direct_plan->f_hat
   *  updates this_iplan->p_hat_iter
   */
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w,
	    this_iplan->r_iter, this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
	  this_iplan->direct_plan->M); 
  
  nfft_adjoint(this_iplan->direct_plan);
  
  updatec_axpy(this_iplan->p_hat_iter, this_iplan->beta_iter,
	       this_iplan->direct_plan->f_hat, this_iplan->direct_plan->N_L);
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->dot_p_hat_iter = 
      dotproductc_w(this_iplan->p_hat_iter, this_iplan->w_hat,
		    this_iplan->direct_plan->N_L);
  else
    this_iplan->dot_p_hat_iter =
      dotproductc(this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);
}

void infft_loop_one_step(infft_plan *this_iplan)
{
  if(this_iplan->infft_flags & LANDWEBER)
    infft_loop_one_step_landweber(this_iplan);

  if(this_iplan->infft_flags & STEEPEST_DESCENT)
    infft_loop_one_step_steepest_descent(this_iplan);
  
  if(this_iplan->infft_flags & CGNR_E)
    infft_loop_one_step_cgnr_e(this_iplan);
  
  if(this_iplan->infft_flags & CGNE_R)
    infft_loop_one_step_cgne_r(this_iplan);
}

void infft_finalize(infft_plan *this_iplan)
{
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    fftw_free(this_iplan->w);
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    fftw_free(this_iplan->w_hat);
  
  if(this_iplan->infft_flags & ITERATE_2nd)
    fftw_free(this_iplan->f_hat_iter_2nd);
  
  if(this_iplan->infft_flags & CGNR_E)
    {
      fftw_free(this_iplan->v_iter);
      fftw_free(this_iplan->z_hat_iter);
    }
  
  if(this_iplan->infft_flags & STEEPEST_DESCENT)
    fftw_free(this_iplan->v_iter);
  
  fftw_free(this_iplan->p_hat_iter);
  fftw_free(this_iplan->f_hat_iter);

  fftw_free(this_iplan->r_iter);
  fftw_free(this_iplan->given_f);
}
