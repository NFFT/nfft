#include "nnfft.h"
#include "window_defines.h"

#define MACRO_nndft_init_result_trafo memset(f,0,this->M_total*sizeof(fftw_complex));
#define MACRO_nndft_init_result_conjugated MACRO_nndft_init_result_trafo
#define MACRO_nndft_init_result_adjoint memset(f_hat,0,this->N_total*               \
                                              sizeof(fftw_complex));
#define MACRO_nndft_init_result_transposed MACRO_nndft_init_result_adjoint

#define MACRO_nndft_sign_trafo      (+2.0*PI)
#define MACRO_nndft_sign_conjugated (-2.0*PI)
#define MACRO_nndft_sign_adjoint    (-2.0*PI)
#define MACRO_nndft_sign_transposed (+2.0*PI)

#define MACRO_nndft_compute_trafo {                                             \
  (*fj0)+=(*f_hat_k0)*cos_omega+(*f_hat_k1)*sin_omega;                          \
  (*fj1)+=(*f_hat_k1)*cos_omega-(*f_hat_k0)*sin_omega;                          \
}
#define MACRO_nndft_compute_conjugated MACRO_nndft_compute_trafo

#define MACRO_nndft_compute_adjoint {                                           \
  (*f_hat_k0)+=(*fj0)*cos_omega+(*fj1)*sin_omega;                               \
  (*f_hat_k1)+=(*fj1)*cos_omega-(*fj0)*sin_omega;                               \
}
#define MACRO_nndft_compute_transposed MACRO_nndft_compute_adjoint

#define MACRO_nndft(which_one)                                                 \
void nndft_ ## which_one (nnfft_plan *this)                                    \
{                                                                              \
  int j;                                /**< index over all nodes (time)     */\
  int t;                                /**< index for dimensions            */\
  int l;                                /**< index over all nodes (fouier)   */\
  fftw_complex *f_hat, *f;              /**< dito                            */\
  double *f_hat_k0, *f_hat_k1;          /**< actual Fourier coefficient      */\
  double *fj0,*fj1;                     /**< actual sample                   */\
  double omega;                         /**< sign times 2*pi*k*x             */\
  double sin_omega,cos_omega;           /**< just for better reading         */\
                                                                               \
  f_hat=this->f_hat; f=this->f;                                                \
                                                                               \
  MACRO_nndft_init_result_ ## which_one                                        \
                                                                               \
  for(j=0, fj0=&f[0][0], fj1=&f[0][1];                                         \
          j<this->M_total; j++, fj0+=2, fj1+=2)                                     \
  {                                                                            \
    for(l=0, f_hat_k0= &f_hat[0][0], f_hat_k1= &f_hat[0][1];                   \
      l<this->N_total; l++, f_hat_k0+=2, f_hat_k1+=2)                               \
    {                                                                          \
      omega=0.0;                                                               \
      for(t = 0; t<this->d; t++)                                               \
        omega+=this->v[l*this->d+t] * this->x[j*this->d+t] * this->N[t];       \
                                                                               \
      omega*= MACRO_nndft_sign_ ## which_one;                                  \
                                                                               \
      cos_omega= cos(omega);                                                   \
      sin_omega= sin(omega);                                                   \
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
void nnfft_uo(nnfft_plan *this,int j,int *up,int *op,int act_dim)
{
  double c;
  int u,o;

  c = this->v[j*this->d+act_dim] * this->n[act_dim];
  u = c; o = c;
  if(c < 0)                  
    u = u-1;                  
  else
    o = o+1;
  
  u = u - (this->m); o = o + (this->m);
  /* make the matrix B nonperiodic */

  up[0]=u; op[0]=o;
}

/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$B, B^{\rm T}\f$
 */
#define MACRO_nnfft_B_init_result_A memset(f,0,this->M_total*sizeof(fftw_complex));
#define MACRO_nnfft_B_init_result_T memset(g,0,this->aN1_L*sizeof(fftw_complex));

#define MACRO_nnfft_B_PRE_FULL_PSI_compute_A {                                  \
  (*fj0)+=this->psi[ix]*g[this->psi_index_g[ix]][0];                           \
  (*fj1)+=this->psi[ix]*g[this->psi_index_g[ix]][1];                           \
}

#define MACRO_nnfft_B_PRE_FULL_PSI_compute_T {                                  \
  g[this->psi_index_g[ix]][0]+=this->psi[ix]*(*fj0);                           \
  g[this->psi_index_g[ix]][1]+=this->psi[ix]*(*fj1);                           \
}

#define MACRO_nnfft_B_compute_A {                                               \
  (*fj0) += phi_prod[this->d]* g[ll_plain[this->d]][0];                        \
  (*fj1) += phi_prod[this->d]* g[ll_plain[this->d]][1];                        \
}

#define MACRO_nnfft_B_compute_T {                                               \
  g[ll_plain[this->d]][0] += phi_prod[this->d]* (*fj0);                        \
  g[ll_plain[this->d]][1] += phi_prod[this->d]* (*fj1);                        \
}

  /* Gewicht, d.h., Nachkommaanteil y-y_u im Speicher halten!!! */
#define MACRO_with_PRE_LIN_PSI (this->psi[(this->K+1)*t2+y_u[t2]]*             \
                                (y_u[t2]+1-y[t2]) +                            \
                                this->psi[(this->K+1)*t2+y_u[t2]+1]*           \
                                (y[t2]-y_u[t2])) 
#define MACRO_with_PRE_PSI     this->psi[(j*this->d+t2)*(2*this->m+2)+lj[t2]]
#define MACRO_without_PRE_PSI  PHI(-this->v[j*this->d+t2]+                      \
                               ((double)l[t2])/this->N1[t2], t2)

#define MACRO_init_uo_l_lj_t {                                                 \
  for(t = this->d-1; t>=0; t--)                                                \
    {                                                                          \
      nnfft_uo(this,j,&u[t],&o[t],t);                                           \
      l[t] = u[t];                                                             \
      lj[t] = 0;                                                               \
    } /* for(t) */                                                             \
  t++;                                                                         \
}

#define MACRO_update_with_PRE_PSI_LIN {                                        \
  for(t2=t; t2<this->d; t2++)                                                  \
    {                                                                          \
      y[t2] = fabs(((-this->N1[t2]*this->v[j*this->d+t2]+(double)l[t2])          \
          * ((double)this->K))/(this->m+1));                                   \
      y_u[t2] = (int)y[t2];                                                    \
    } /* for(t2) */                                                            \
}

#define MACRO_update_phi_prod_ll_plain(which_one) {                            \
  for(t2=t; t2<this->d; t2++)                                                  \
    {                                                                          \
      phi_prod[t2+1]=phi_prod[t2]* MACRO_ ## which_one;                        \
      ll_plain[t2+1]=ll_plain[t2]*this->aN1[t2] +(l[t2]+this->aN1[t2]*3/2)%this->aN1[t2];\
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

#define MACRO_nnfft_B(which_one)                                                \
inline void nnfft_B_ ## which_one (nnfft_plan *this)                           \
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
  double y[this->d];                                                           \
  int y_u[this->d];                                                            \
                                                                               \
  f=this->f_hat; g=this->F;                                                    \
                                                                               \
  MACRO_nnfft_B_init_result_ ## which_one                                                 \
                                                                               \
  if((this->nnfft_flags & PRE_PSI)&&(this->nnfft_flags & PRE_FULL_PSI))        \
    {                                                                          \
      for(ix=0, j=0, fj0=&f[0][0], fj1=&f[0][1]; j<this->N_total; j++,fj0+=2,fj1+=2)\
        for(l_L=0; l_L<this->psi_index_f[j]; l_L++, ix++)                      \
          MACRO_nnfft_B_PRE_FULL_PSI_compute_ ## which_one;                                \
      return;                                                                  \
    }                                                                          \
                                                                               \
  phi_prod[0]=1;                                                               \
  ll_plain[0]=0;                                                               \
                                                                               \
  for(t=0,lprod = 1; t<this->d; t++)                                           \
    lprod *= (2*this->m+2);                                                    \
                                                                               \
  if(this->nnfft_flags & PRE_PSI)                                               \
    {                                                                          \
      for(j=0, fj0=&f[0][0], fj1=&f[0][1]; j<this->N_total; j++, fj0+=2, fj1+=2)     \
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
  if(this->nnfft_flags & PRE_LIN_PSI)                                           \
    {                                                                          \
      for(j=0, fj0=&f[0][0], fj1=&f[0][1]; j<this->N_total; j++, fj0+=2, fj1+=2)     \
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
  for(j=0, fj0=&f[0][0], fj1=&f[0][1]; j<this->N_total; j++, fj0+=2, fj1+=2)         \
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

inline void nnfft_D (nnfft_plan *this){
  int j,t;
  double tmp;
  
  if(this->nnfft_flags & PRE_PHI_HUT) {
    for(j=0; j<this->M_total; j++) {
      this->f[j][0] *= this->c_phi_inv[j];
      this->f[j][1] *= this->c_phi_inv[j];
    }
  } else {
    for(j=0; j<this->M_total; j++)
    {
      tmp = 1.0;
      /* multiply with N1, because x was modified */
      for(t=0; t<this->d; t++)
        tmp*= 1.0 /((PHI_HUT(this->x[this->d*j + t]*((double)this->N[t]),t)) );
      this->f[j][0] *= tmp ;
      this->f[j][1] *= tmp ;
    }
  }
}

/** user routines
 */
void nnfft_trafo(nnfft_plan *this)
{
  int j,t;

  nnfft_B_T(this);

  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] / ((double)this->sigma[t]);
    }
  }
  
  this->direct_plan->f = this->f;
  nfft_trafo(this->direct_plan);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] * ((double)this->sigma[t]);
    }
  }
  
  nnfft_D(this);
} /* nnfft_trafo */

void nnfft_conjugated(nnfft_plan *this)
{
  int j,t;
  
  nnfft_B_T(this);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] / ((double)this->sigma[t]);
    }
  }
  
  this->direct_plan->f = this->f;
  nfft_conjugated(this->direct_plan);

  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] * ((double)this->sigma[t]);
    }
  }  
  
  nnfft_D(this);
} /* nnfft_conjugated */

void nnfft_adjoint(nnfft_plan *this)
{
  int j,t;
  
  nnfft_D(this);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] / ((double)this->sigma[t]);
    }
  }
  
  this->direct_plan->f = this->f;
  nfft_adjoint(this->direct_plan);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] * ((double)this->sigma[t]);
    }
  }  
  
  nnfft_B_A(this);
} /* nnfft_adjoint */

void nnfft_transposed(nnfft_plan *this)
{
  int j,t;
  
  nnfft_D(this);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] / ((double)this->sigma[t]);
    }
  }
  
  this->direct_plan->f = this->f;
  nfft_transposed(this->direct_plan);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] * ((double)this->sigma[t]);
    }
  }  
  
  nnfft_B_A(this);
} /* nnfft_transposed */


/** initialisation of direct transform 
 */
void nnfft_precompute_phi_hut(nnfft_plan *this)
{
  int j;                                /**< index over all frequencies       */
  int t;                                /**< index over all dimensions        */
  double tmp;

  this->c_phi_inv= (double*)fftw_malloc(this->M_total*sizeof(double));
  
  for(j=0; j<this->M_total; j++)
    {
      tmp = 1.0;
      for(t=0; t<this->d; t++)
        tmp*= 1.0 /(PHI_HUT(this->x[this->d*j + t]*((double)this->N[t]),t));
      this->c_phi_inv[j]=tmp;
    }
} /* nnfft_phi_hut */


/** create a lookup table, but NOT for each node
 *  good idea K=2^xx
 *  TODO: estimate K, call from init
 */
void nnfft_precompute_lin_psi(nnfft_plan *this, int K)
{
  int t;                                /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */  
  double step;                          /**< step size in [0,(m+1)/n]         */

  this->K=K;
  this->psi = (double*) malloc((this->K+1)*this->d*sizeof(double));
  
  nfft_precompute_lin_psi(this->direct_plan,K);
  
  for (t=0; t<this->d; t++)
    {
      step=((double)(this->m+1))/(K*this->N1[t]);
      for(j=0;j<=K;j++)
        {
          this->psi[(K+1)*t + j] = PHI(j*step,t);
        } /* for(j) */
    } /* for(t) */
}

void nnfft_precompute_psi(nnfft_plan *this)
{
  int t;                                /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int l;                                /**< index u<=l<=o                    */
  int lj;                               /**< index 0<=lj<u+o+1                */
  int u, o;                             /**< depends on v_j                   */
  
  for (t=0; t<this->d; t++)
    for(j=0;j<this->N_total;j++)
      {
        nnfft_uo(this,j,&u,&o,t);
        
        for(l=u, lj=0; l <= o; l++, lj++)
          this->psi[(j*this->d+t)*(2*this->m+2)+lj]=
            (PHI((-this->v[j*this->d+t]+((double)l)/((double)this->N1[t])),t));
      } /* for(j) */
      
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] / ((double)this->sigma[t]);
    }
  }
  
  nfft_precompute_psi(this->direct_plan);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] * ((double)this->sigma[t]);
    }
  }
  /* for(t) */
} /* nfft_precompute_psi */

/** computes all entries of B explicitly */
void nnfft_full_psi(nnfft_plan *this, double eps)
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
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] / ((double)this->sigma[t]);
    }
  }
  
  nfft_full_psi(this->direct_plan,eps);
  
  for(j=0;j<this->M_total;j++) {  
    for(t=0;t<this->d;t++) {
      this->x[j*this->d+t]= this->x[j*this->d+t] * ((double)this->sigma[t]);
    }
  }
  
  phi_prod[0]=1;
  ll_plain[0]=0;

  size_psi=this->N_total;
  index_f = (int*) malloc(this->N_total*sizeof(int));
  index_g = (int*) malloc(size_psi*sizeof(int));
  new_psi = (double*) malloc(size_psi*sizeof(double));
  
  for(t=0,lprod = 1; t<this->d; t++)
    {
      lprod *= 2*this->m+2;
      eps=eps*(PHI(0,t));
    }
  
  for(ix=0,ix_old=0,j=0; j<this->N_total; j++)
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
                  size_psi+=this->N_total;
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
}

void nnfft_init_help(nnfft_plan *this, int m2, int *N2, unsigned nfft_flags, unsigned fftw_flags)
{
  int t;                                /**< index over all dimensions        */
  
  this->aN1 = (int*) fftw_malloc(this->d*sizeof(int));
  
  this->a = (double*) fftw_malloc(this->d*sizeof(double));
  
  this->sigma = (double*) fftw_malloc(this->d*sizeof(double));
  
  this->n = this->N1;
  
  this->aN1_L=1;
  
  for(t = 0; t<this->d; t++) {
    this->a[t] = 1.0 + (2.0*((double)this->m))/((double)this->N1[t]);
    this->aN1[t] = this->a[t] * ((double)this->N1[t]);
    if(this->aN1[t]%2 != 0)
      this->aN1[t] = this->aN1[t] +1;
      
    this->aN1_L*=this->aN1[t];
    this->sigma[t] = ((double) this->N1[t] )/((double) this->N[t]);;
  }
  
  WINDOW_HELP_INIT
  
  if(this->nnfft_flags & MALLOC_X)
    this->x = (double*)fftw_malloc(this->d*this->M_total*
                                        sizeof(double));
                                        
  if(this->nnfft_flags & MALLOC_V)
    this->v = (double*)fftw_malloc(this->d*this->N_total*
                                        sizeof(double));

  if(this->nnfft_flags & MALLOC_F_HAT)
    this->f_hat = (fftw_complex*)fftw_malloc(this->N_total*
                                                  sizeof(fftw_complex));
  if(this->nnfft_flags & MALLOC_F)
    this->f=(fftw_complex*)fftw_malloc(this->M_total*sizeof(fftw_complex));
    
  /* NO FFTW_MALLOC HERE */
  if(this->nnfft_flags & PRE_PSI)
    this->psi = (double*) malloc(this->N_total*this->d*
                                           (2*this->m+2)*sizeof(double));
                                           
  this->direct_plan = (nfft_plan*) malloc(sizeof(nfft_plan));
    
  nfft_init_specific(this->direct_plan, this->d, this->aN1, this->M_total, N2,
                        m2, nfft_flags, fftw_flags);
                        
  this->direct_plan->x = this->x;
  this->direct_plan->f = this->f;
  this->F = this->direct_plan->f_hat;
}

void nnfft_init_specific(nnfft_plan *this, int d, int N_total, int M_total, int *N, int *N1,
                        int m, unsigned nnfft_flags)
{
  int t;                             /**< index over all dimensions        */
  
  int N2[d];
  unsigned nfft_flags;
  unsigned fftw_flags;
  
  this->d= d;
  this->M_total= M_total;  
  this->N_total= N_total;
  this->m= m;
  this->nnfft_flags= nnfft_flags;
  fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;
  nfft_flags= PRE_PHI_HUT| MALLOC_F_HAT| FFTW_INIT| FFT_OUT_OF_PLACE;
  
  if(this->nnfft_flags & PRE_PSI)
    nfft_flags = nfft_flags | PRE_PSI;
    
  if(this->nnfft_flags & PRE_FULL_PSI)
    nfft_flags = nfft_flags | PRE_FULL_PSI;
    
  if(this->nnfft_flags & PRE_LIN_PSI)
    nfft_flags = nfft_flags | PRE_LIN_PSI;
  
  this->N = (int*) fftw_malloc(this->d*sizeof(int));
  this->N1 = (int*) fftw_malloc(this->d*sizeof(int));
  
  for(t=0; t<d; t++) {
    this->N[t] = N[t];
    this->N1[t] = N1[t];    
    N2[t] = 2*next_power_of_2((int)((double)this->N1[t])* 
           (1.0+2.0*((double)this->m)/((double)this->N1[t])));
  }
  nnfft_init_help(this,m,N2,nfft_flags,fftw_flags);  
}

void nnfft_init(nnfft_plan *this, int d, int N_total, int M_total, int *N)
{
  int t;                            /**< index over all dimensions        */

  int N2[d];
  unsigned nfft_flags;
  unsigned fftw_flags;

  this->d = d;
  this->M_total = M_total;
  this->N_total = N_total;

  WINDOW_HELP_ESTIMATE_m;
  
  this->N = (int*) fftw_malloc(this->d*sizeof(int));  
  this->N1 = (int*) fftw_malloc(this->d*sizeof(int));
  
  for(t=0; t<d; t++) {
    this->N[t] = N[t];
    this->N1[t] = 2*next_power_of_2(this->N[t]);
    N2[t] = 2*next_power_of_2((int)((double)this->N1[t])* 
            (1.0+2.0*((double)this->m)/((double)this->N1[t])));
  }
  this->nnfft_flags=PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F;
  nfft_flags= PRE_PSI| PRE_PHI_HUT| MALLOC_F_HAT| FFTW_INIT| FFT_OUT_OF_PLACE;
  
  fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;
  
  nnfft_init_help(this,this->m,N2,nfft_flags,fftw_flags);    
}

void nnfft_finalize(nnfft_plan *this)
{ 
  nfft_finalize(this->direct_plan);
  
  free(this->direct_plan);

  free(this->aN1);
  free(this->N);
  free(this->N1);
  
  /* NO FFTW_FREE HERE */
  if((this->nnfft_flags & PRE_PSI) || (this->nnfft_flags & PRE_LIN_PSI))
    {
      if(this->nnfft_flags & PRE_FULL_PSI)
        {
          free(this->psi_index_g);
          free(this->psi_index_f);
        }

      free(this->psi);
    }
      
  if(this->nnfft_flags & PRE_PHI_HUT)
    fftw_free(this->c_phi_inv);
  
  if(this->nnfft_flags & MALLOC_F)
    fftw_free(this->f);

  if(this->nnfft_flags & MALLOC_F_HAT)
    fftw_free(this->f_hat);

  if(this->nnfft_flags & MALLOC_X)
    fftw_free(this->x);
    
  if(this->nnfft_flags & MALLOC_V)
    fftw_free(this->v);
}
