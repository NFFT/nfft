/**
 * Library.
 * Includes simple and fast computation of the NFCT (direct problem)
 * authors: S. Klatt (c) 2004-2006
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"
#include "nfft3.h"

/** 
 * handy shortcuts
 **/
#define NFCT_DEFAULT_FLAGS   PRE_PHI_HUT|\
                             PRE_PSI|\
                             MALLOC_X|\
                             MALLOC_F_HAT|\
                             MALLOC_F|\
                             FFTW_INIT|\
                             FFT_OUT_OF_PLACE

#define FFTW_DEFAULT_FLAGS   FFTW_ESTIMATE|\
                             FFTW_DESTROY_INPUT

#define NFCT_SUMMANDS ( 2 * ths->m + 2)
#define NODE(p,r) ( ths->x[(p) * ths->d + (r)])

#define MACRO_ndct_init_result_trafo       \
  memset( f, 0, ths->M_total * sizeof( double));
#define MACRO_ndct_init_result_adjoint     \
  memset( f_hat, 0, ths->N_total * sizeof( double));


#define MACRO_nfct_D_init_result_A         \
  memset(g_hat,   0, nfct_prod_int(ths->n, ths->d) * sizeof( double));
#define MACRO_nfct_D_init_result_T         \
  memset(f_hat, 0, ths->N_total * sizeof( double));

#define MACRO_nfct_B_init_result_A         \
  memset(f, 0, ths->M_total * sizeof( double));
#define MACRO_nfct_B_init_result_T         \
  memset(g, 0, nfct_prod_int( ths->n, ths->d) * sizeof( double));


#define NFCT_PRE_WINFUN( d)  ths->N[d] = 2 * ths->N[d];             \
                             ths->n[d] = nfct_fftw_2N( ths->n[d]);

#define NFCT_POST_WINFUN( d) ths->N[d] = (int)(0.5 * ths->N[d]);    \
                             ths->n[d] = nfct_fftw_2N_rev( ths->n[d]);


#define NFCT_WINDOW_HELP_INIT  WINDOW_HELP_INIT


double nfct_phi_hut( nfct_plan *ths, int k, int d)
{
  NFCT_PRE_WINFUN( d);
  double phi_hut_tmp = PHI_HUT( k, d);
  NFCT_POST_WINFUN( d);

  return phi_hut_tmp;
}

double nfct_phi( nfct_plan *ths, double x, int d)
{
  NFCT_PRE_WINFUN( d); 
  double phi_tmp = PHI( x, d); 
  NFCT_POST_WINFUN( d); 

  return phi_tmp; 
} 

int nfct_fftw_2N( int n)
{
  return 2 * ( n - 1);
}

int nfct_fftw_2N_rev( int n)
{
  return (int)(0.5 * n) + 1;
}


#define MACRO_with_cos_vec       cos_vec[t][ka[t]]
#define MACRO_without_cos_vec    cos( 2.0 * PI * ka[t] * NODE(j,t))


#define MACRO_with_PRE_PHI_HUT      ths->c_phi_inv[t][kg[t]]; 
#define MACRO_compute_PHI_HUT_INV   (1.0 / (nfct_phi_hut( ths, kg[t], t))) 

#define MACRO_with_PRE_PSI   ths->psi[(j * ths->d + t) * NFCT_SUMMANDS + lc[t]] 
#define MACRO_compute_PSI    nfct_phi( ths, NODE(j,t) - (( double)(lc[t] + lb[t])) / (double)nfct_fftw_2N( ths->n[t]), t) 




/** direct computation of non equispaced cosine transforms
*  ndct_trafo,  ndct_adjoint
*  require O(M N^d) arithemtical operations
*
* direct computation of the ndct_trafo, formula (1.1)
* ndct_trafo:
* for j=0,...,M-1                                                             
*  f[j] = sum_{k in I_N^d} f_hat[k] * cos(2 (pi) k x[j])
*
* direct computation of the ndft_adjoint, formula (1.2)
* ndct_adjoint:
* for k in I_N^d
*  f_hat[k] = sum_{j=0}^{M-1} f[j] * cos(2 (pi) k x[j])
*/


#define MACRO_ndct_malloc__cos_vec                                      \
                                                                        \
  double **cos_vec;                                                     \
  cos_vec = (double**) malloc( ths->d * sizeof( double*));              \
  for( t = 0; t < ths->d; t++)                                        	\
    cos_vec[t] = (double*) malloc( ths->N[t] * sizeof( double));        \




#define MACRO_ndct_free__cos_vec                                        \
{                                                               	\
  /* free allocated memory */                                         	\
  for( t = 0; t < ths->d; t++)                                        	\
    free( cos_vec[t]);                    	                        \
  free( cos_vec);                                                     	\
}



#define MACRO_ndct_init__cos_vec                                        \
{                                                                       \
  for( t = 0; t < ths->d; t++)                                          \
  {                                                                     \
    cos_vec[t][0] = 1.0;                                                \
    cos_vec[t][1] = cos( 2.0 * PI * NODE(j,t));	                        \
    for( k = 2; k < ths->N[t]; k++)				        \
      cos_vec[t][k] = 2.0 * cos_vec[t][1] * cos_vec[t][k-1]           	\
                      - cos_vec[t][k-2];                              	\
  }                                                                     \
}



#define MACRO_ndct_init__k__cos_k( which_one)                           \
{                                                                       \
  cos_k[0] = 1.0;                                                       \
  for( t = 0; t < ths->d; t++)                                          \
    ka[t] = 0;                                                          \
                                                                        \
    for( t = 0; t < ths->d; t++)                                        \
    {                                                                   \
      cos_k[t+1] = cos_k[t] * MACRO_ ##which_one;                       \
    }                                                                   \
}



#define MACRO_ndct_count__k__cos_k( which_one)                          \
{                                                                       \
  ka[ths->d-1]++;                                                       \
  i = ths->d - 1;                                                       \
  while( ( ka[i] == ths->N[i]) && ( i>0))                               \
  {                                                                     \
    ka[i - 1]++;                                                        \
    ka[i] = 0;                                                          \
                                                                        \
    i--;                                                                \
  }                                                                     \
  for( t = i; t < ths->d; t++)                                          \
    cos_k[t+1] = cos_k[t] * MACRO_ ##which_one;                         \
}



#define MACRO_ndct_compute__trafo                                       \
{                                                                       \
  f[j] += f_hat[k] * cos_k[ths->d];	                                    \
}

#define MACRO_ndct_compute__adjoint                                     \
{                                                                       \
  f_hat[k] += f[j] * cos_k[ths->d];                                     \
}

/* slow (trafo) transform */
#define MACRO_ndct(which_one)                                           \
  void ndct_ ## which_one ( nfct_plan *ths)                             \
  {                                                                     \
    int j, k, t, i;                                                     \
    int ka[ths->d];                                                     \
    double cos_k[ths->d+1];                                             \
                                                                        \
    double *f     = ths->f;                                             \
    double *f_hat = ths->f_hat;                                         \
                                                                        \
    MACRO_ndct_init_result_ ## which_one;			                    \
    if( ths->d == 1)                                                    \
      for( j = 0; j < ths->M_total; j++)                                \
      {                                                                 \
        for( k = 0; k < ths->N[0]; k++)                                 \
       	{                                                               \
	  cos_k[ths->d] = cos( 2.0 * PI * k * NODE(j,0));               \
          MACRO_ndct_compute__ ## which_one;                            \
        }							        \
      }                                                                 \
    else   /*FIXME: remove slow slow ... */                             \
      if( 1 == 0)                                                       \
        /* slow ndct */                                                 \
        for( j = 0; j < ths->M_total; j++)                              \
        {                                                               \
          MACRO_ndct_init__k__cos_k(without_cos_vec);                   \
                                                                        \
          for( k = 0; k < ths->N_total; k++)                            \
          {                                                             \
            MACRO_ndct_compute__ ## which_one;                          \
                                                                        \
            MACRO_ndct_count__k__cos_k(without_cos_vec);                \
          }     					                \
        }                                                               \
      else                                                              \
        {                                                               \
          /* fast ndct_trafo */                                         \
          MACRO_ndct_malloc__cos_vec;                                   \
                                                                        \
          for( j = 0; j < ths->M_total; j++)                            \
          {                                                             \
                                                                        \
            MACRO_ndct_init__cos_vec;                                   \
                                                                        \
            MACRO_ndct_init__k__cos_k(with_cos_vec);                    \
                                                                        \
            for( k = 0; k < ths->N_total; k++)                          \
            {                                                           \
                                                                        \
              MACRO_ndct_compute__ ## which_one;                        \
                							\
              MACRO_ndct_count__k__cos_k(with_cos_vec);                 \
            }                                                           \
          }                                                             \
          MACRO_ndct_free__cos_vec;                                     \
        }                                                               \
} /* ndct_{trafo, adjoint} */


MACRO_ndct(trafo)
MACRO_ndct(adjoint)




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



#define MACRO_nfct__lower_boundary( j,act_dim)                                  \
{                                                                               \
  lb[(act_dim)] =                                                               \
    (int)( NODE((j),(act_dim)) * nfct_fftw_2N( ths->n[(act_dim)])) - ths->m;    \
}

#define MACRO_nfct_D_compute_A          					\
{                                                                               \
  g_hat[kg_plain[ths->d]] = f_hat[k_L] * c_phi_inv_k[ths->d];	                \
}

#define MACRO_nfct_D_compute_T                                                  \
{                                                                               \
  f_hat[k_L] = g_hat[kg_plain[ths->d]] * c_phi_inv_k[ths->d];                   \
}



#define MACRO_init__kg          						\
{                                                                               \
  for( t = 0; t < ths->d; t++)					                \
    kg[t] = 0;							                \
								                \
  i = 0;								        \
}


#define MACRO_count__kg                                                         \
{                                                                               \
                                                                                \
  kg[ths->d - 1]++;                                                             \
  i = ths->d - 1;                                                               \
  while( ( kg[i] == ths->N[i]) && ( i > 0))                                     \
  {                                                                             \
    kg[i - 1]++;                                                                \
    kg[i] = 0;                                                                  \
                                                                                \
    i--;                                                                        \
  }                                                                             \
}


#define MACRO_update__phi_inv_k__kg_plain( what_kind, which_phi)        	\
{                                                                               \
  for( t = i; t < ths->d; t++)                                                  \
  {                                                                             \
    MACRO__c_phi_inv_k__ ## what_kind( which_phi);                              \
    kg_plain[t+1] = kg_plain[t] * ths->n[t] + kg[t];                            \
  }                                                                             \
}


#define MACRO__c_phi_inv_k__A( which_phi)                                       \
{                                                                               \
  if( kg[t] == 0)                                                               \
  {                                                                             \
    c_phi_inv_k[t+1] = c_phi_inv_k[t] * MACRO_ ## which_phi;                 	\
  }                                                                             \
  else                                                                          \
  {                                                                             \
    c_phi_inv_k[t+1] = 0.5 * c_phi_inv_k[t] * MACRO_ ## which_phi;         	\
  }                                                                             \
}


#define MACRO__c_phi_inv_k__T( which_phi)                                       \
{                                                                               \
  c_phi_inv_k[t+1] = c_phi_inv_k[t] * MACRO_ ## which_phi;                  	\
}



#define MACRO_nfct_D(which_one)                                                 \
static inline void nfct_D_ ## which_one (nfct_plan *ths)                               \
{			                                                        \
                                                                                \
  int k_L;                        /**< plain index                    */        \
                                                                                \
  int i, t;                                                                     \
  int kg[ths->d];                 /**< multi index in g_hat,c_phi     */        \
  double c_phi_inv_k[ths->d+1];   /**< postfix product of PHI_HUT_INV */        \
  int kg_plain[ths->d+1];         /**< postfix plain index            */        \
                                                                                \
  double *g_hat, *f_hat;          /**< local copy                     */        \
                                                                                \
  g_hat = ths->g_hat;                                                           \
  f_hat = ths->f_hat;                                                           \
                                                                                \
  MACRO_nfct_D_init_result_ ## which_one                                        \
                                                                                \
  c_phi_inv_k[0] = 1.0;                                                         \
  kg_plain[0]    =   0;                                                         \
                                                                                \
  MACRO_init__kg;                                                               \
                                                                                \
  if( ths->nfct_flags & PRE_PHI_HUT)                                            \
    for( k_L = 0; k_L < ths->N_total; k_L++)                                    \
    {                                                                           \
      MACRO_update__phi_inv_k__kg_plain( which_one, with_PRE_PHI_HUT);          \
                                                                                \
      MACRO_nfct_D_compute_ ## which_one;                                       \
                                                                                \
      MACRO_count__kg;                                                          \
                                                                                \
    } /* for(k_L) */                                                            \
  else                                                                          \
    for( k_L = 0; k_L < ths->N_total; k_L++)                                    \
    {                                                                           \
      MACRO_update__phi_inv_k__kg_plain( which_one,compute_PHI_HUT_INV);        \
                                                                                \
      MACRO_nfct_D_compute_ ## which_one;                                       \
                                                                                \
      MACRO_count__kg;                                                          \
                                                                                \
    } /* for(k_L) */                                                            \
} /* nfct_D */

MACRO_nfct_D(A)
MACRO_nfct_D(T)







/** sub routines for the fast transforms
*  matrix vector multiplication with \f$B, B^{\rm T}\f$
*/ 
#define MACRO_nfct_B_PRE_FULL_PSI_compute_A                                     \
{                                                                               \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]];                              \
}

#define MACRO_nfct_B_PRE_FULL_PSI_compute_T                                     \
{                                                                               \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj);                              \
}



#define MACRO_nfct_B_compute_A                                                  \
{                                                                               \
  (*fj) += phi_tilde[ths->d] * g[lg_plain[ths->d]];                             \
}

#define MACRO_nfct_B_compute_T                                                  \
{                                                                               \
  g[lg_plain[ths->d]] += phi_tilde[ths->d] * (*fj);                             \
}


#define MACRO_compute_lg_offset__count_lg(i0)   				\
{                                                                               \
  /* determine index in g-array corresponding to lb[(i0)] */                    \
  if( lb[(i0)] < 0)                                                             \
    lg_offset[(i0)] =                                                           \
      (lb[(i0)] % nfct_fftw_2N( ths->n[(i0)])) + nfct_fftw_2N( ths->n[(i0)]);   \
  else                                                                          \
    lg_offset[(i0)] = lb[(i0)] % (nfct_fftw_2N( ths->n[(i0)]));                 \
    if( lg_offset[(i0)] >= ths->n[(i0)])                                        \
      lg_offset[(i0)] = -( nfct_fftw_2N( ths->n[(i0)]) - lg_offset[(i0)]);      \
}

#define MACRO_set__lg__to__lg_offset                                            \
{                                                                               \
  if( lg_offset[i] <= 0)                                                        \
  {                                                                             \
    lg[i] = -lg_offset[i];                                                      \
    count_lg[i] = -1;                                                           \
  }                                                                             \
  else                                                                          \
  {                                                                             \
    lg[i] = +lg_offset[i];                                                      \
    count_lg[i] = +1;                                                           \
  }                                                                             \
}

#define MACRO_count__lg(dim)                                                    \
{                                                                               \
                                                                                \
  /* turn around if we hit one of the boundaries */                             \
  if( (lg[(dim)] == 0) || (lg[(dim)] == ths->n[(dim)]-1))                       \
    count_lg[(dim)] *= -1;                                                      \
  /* move array index */                                                        \
  lg[(dim)] += count_lg[(dim)];                                                 \
}




#define MACRO_init_lb_lg_lc                                                     \
{                                                                               \
  for( i = 0; i < ths->d; i++)                                                  \
  {                                                                             \
    MACRO_nfct__lower_boundary( j, i);                                          \
                                                                                \
    MACRO_compute_lg_offset__count_lg( i);                                      \
    MACRO_set__lg__to__lg_offset;                                               \
                                                                                \
    /* counter for lg */                                                        \
    lc[i] = 0;                                                                  \
   }                                                                            \
                                                                                \
   i = 0;                                                                       \
}



#define MACRO_count__lg_lc                                                      \
{                                                                               \
  MACRO_count__lg( ths->d-1);                                                   \
                                                                                \
  lc[ths->d - 1]++;                                                             \
  i = ths->d - 1;                                                               \
  while( ( lc[i] == NFCT_SUMMANDS) && ( i > 0))                                 \
  {                                                                             \
    lc[i - 1]++;                                                                \
    lc[i] = 0;                                                                  \
                                                                                \
    /* ansonsten lg[i-1] verschieben */                                         \
    MACRO_count__lg( i - 1);                                                    \
    /* lg[i] = anfangswert */                                                   \
    MACRO_set__lg__to__lg_offset;                                               \
                                                                                \
    i--;                                                                        \
  }                                                                             \
}

#define  MACRO_update_phi_tilde_lg_plain( which_one, which_psi)                 \
{                                                                               \
  for( t = i; t < ths->d; t++)                                                  \
  {                                                                             \
    MACRO__phi_tilde__ ## which_one( which_psi);                                \
    lg_plain[t+1]  = lg_plain[t]  * ths->n[t] + lg[t];                          \
  }                                                                             \
}

#define MACRO__phi_tilde__A( which_psi)                                         \
{                                                                               \
  phi_tilde[t+1] = phi_tilde[t] * MACRO_ ## which_psi;                        	\
}

#define MACRO__phi_tilde__T( which_psi)                                         \
{                                                                               \
  if( lg[t] == 0 || lg[t] == ths->n[t] - 1)                                     \
  {                                                                             \
    phi_tilde[t+1] = phi_tilde[t] * MACRO_ ## which_psi;                        \
  }                                                                             \
  else                                                                          \
  {                                                                             \
    phi_tilde[t+1] = 0.5 * phi_tilde[t] * MACRO_ ## which_psi;           	\
  }                                                                             \
}



#define MACRO_nfct_B( which_one)                                                \
static inline void nfct_B_ ## which_one ( nfct_plan *ths)                              \
{ /* MACRO_nfct_B */                                                            \
  int lb[ths->d];               /**< multi band with respect to x_j */          \
  int j, t, i;                  /**< index nodes, help vars         */          \
  int lprod, l_L, ix;           /**< index one row of B             */          \
  int lc[ths->d];               /**< multi index 0<=lc<2m+2         */          \
  int lg[ths->d];               /**< real index of g in array       */          \
  int lg_offset[ths->d];        /**< offset in g according to u     */          \
  int count_lg[ths->d];         /**< count summands (2m+2)          */          \
  int lg_plain[ths->d+1];       /**< index of g in multi_array      */          \
  double *f, *g;                /**< local copy                     */          \
  double phi_tilde[ths->d+1];   /**< holds values for psi           */          \
  double *fj;                   /**< pointer to final result        */          \
                                                                                \
  f = ths->f; g = ths->g;                                                       \
                                                                                \
  MACRO_nfct_B_init_result_ ## which_one                                        \
                                                                                \
  /* both flags are set */                                                      \
  if(( ths->nfct_flags & PRE_PSI)&&( ths->nfct_flags & PRE_FULL_PSI))           \
  {                                                                             \
    for( ix = 0, j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1)             \
      for( l_L = 0; l_L < ths->psi_index_f[j]; l_L++, ix++)                     \
      {                                                                         \
        MACRO_nfct_B_PRE_FULL_PSI_compute_ ## which_one;                        \
      }                                                                         \
  }                                                                             \
  else                                                                          \
  {                                                                             \
    phi_tilde[0] = 1;                                                           \
    lg_plain[0]  = 0;                                                           \
                                                                                \
    for( t = 0, lprod = 1; t < ths->d; t++)                                     \
      lprod *= NFCT_SUMMANDS;                                                   \
                                                                                \
    /**/                                                                        \
    /* PRE_PSI flag is set */                                                   \
    /**/                                                                        \
    if( ths->nfct_flags & PRE_PSI)                                              \
      for( j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1)                   \
        {                                                                       \
          MACRO_init_lb_lg_lc;                                                  \
                                                                                \
          for( l_L = 0; l_L < lprod; l_L++)                                     \
          {                                                                     \
            MACRO_update_phi_tilde_lg_plain( which_one, with_PRE_PSI);          \
                                                                                \
            MACRO_nfct_B_compute_ ## which_one;                                 \
                                                                                \
            MACRO_count__lg_lc;                                                 \
          } /* for( l_L) */                                                     \
        } /* for( j) */                                                         \
                                                                                \
    /**/                                                                        \
    /*  no PSI flag is set */                                                   \
    /**/                                                                        \
    else                                                                        \
      for( j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1)                   \
      {                                                                         \
        MACRO_init_lb_lg_lc;                                                    \
                                                                                \
        for( l_L = 0; l_L < lprod; l_L++)                                       \
        {                                                                       \
          MACRO_update_phi_tilde_lg_plain( which_one,compute_PSI);              \
                                                                                \
          MACRO_nfct_B_compute_ ## which_one;                                   \
                                                                                \
          MACRO_count__lg_lc;                                                   \
                                                                                \
        } /* for(l_L) */                                                        \
      } /* for(j) */                                                            \
  } /* else( PRE_PSI && FULL_PRE_PSI) */                                        \
} /* nfct_B */

MACRO_nfct_B(A)
MACRO_nfct_B(T)




/** more memory usage, a bit faster */
#define MACRO_nfct_full_psi( which_one)                                         \
void nfct_full_psi__ ## which_one( nfct_plan *ths)                              \
{                                                                               \
  int t, i;                  /**< index over all dimensions        */           \
  int j;                     /**< index over all nodes             */           \
  int l_L;                   /**< plain index 0<=l_L<lprod         */           \
  int lc[ths->d];            /**< multi index 0<=lj<u+o+1          */           \
  int lg_plain[ths->d+1];    /**< postfix plain index              */           \
  int count_lg[ths->d];                                                         \
  int lg_offset[ths->d];                                                        \
  int lg[ths->d];                                                               \
  int lprod;                 /**< 'bandwidth' of matrix B          */           \
  int lb[ths->d];            /**< depends on x_j                   */           \
                                                                                \
  double phi_tilde[ths->d+1];                                                   \
  double eps = ths->nfct_full_psi_eps;                                          \
                                                                                \
  int *index_g, *index_f;                                                       \
  double *new_psi;                                                              \
  int ix, ix_old, size_psi;                                                     \
                                                                                \
  phi_tilde[0] = 1.0;                                                           \
  lg_plain[0]  =   0;                                                           \
                                                                                \
  if( ths->nfct_flags & PRE_PSI)                                                \
  {                                                                             \
    size_psi = ths->M_total;                                                    \
    index_f  =    (int*) malloc( ths->M_total  * sizeof( int));                 \
    index_g  =    (int*) malloc( size_psi * sizeof( int));                      \
    new_psi  = (double*) malloc( size_psi * sizeof( double));                   \
                                                                                \
    for( t = 0,lprod = 1; t < ths->d; t++)                                      \
    {                                                                           \
      lprod *= NFCT_SUMMANDS;                                                   \
      eps *= nfct_phi( ths, 0, t);                                              \
    }                                                                           \
                                                                                \
    for( ix = 0, ix_old = 0, j = 0; j < ths->M_total; j++)                      \
    {                                                                           \
      MACRO_init_lb_lg_lc;                                                      \
                                                                                \
      for( l_L = 0; l_L < lprod; l_L++)                                         \
      {                                                                         \
        MACRO_update_phi_tilde_lg_plain( which_one, with_PRE_PSI);              \
                                                                                \
        if( phi_tilde[ths->d] > eps)                                            \
        {                                       		                \
          index_g[ix] =  lg_plain[ths->d];                                      \
          new_psi[ix] = phi_tilde[ths->d];                                      \
                                                                                \
          ix++;                                                                 \
          if( ix == size_psi)                                                   \
          {                                                                     \
            size_psi += ths->M_total;                                           \
            index_g   =    (int*)realloc( index_g,                              \
                                          size_psi * sizeof( int));             \
            new_psi   = (double*)realloc( new_psi,                              \
                                          size_psi * sizeof( double));          \
          }                                                                     \
        }                                                                       \
                                                                                \
        MACRO_count__lg_lc;                                                     \
                                                                                \
      } /* for(l_L) */                                                          \
                                                                                \
      index_f[j] = ix - ix_old;                                                 \
      ix_old     = ix;                                                          \
                                                                                \
    } /* for(j) */                                                              \
                                                                                \
    free( ths->psi);                                                            \
    size_psi = ix;                                                              \
    ths->size_psi = size_psi;                                                   \
                                                                                \
    index_g =    (int*)realloc( index_g, size_psi * sizeof( int));              \
    new_psi = (double*)realloc( new_psi, size_psi * sizeof( double));           \
                                                                                \
    ths->psi         = new_psi;                                                 \
    ths->psi_index_g = index_g;                                                 \
    ths->psi_index_f = index_f;                                                 \
                                                                                \
  } /* if(PRE_PSI) */                                                           \
}

MACRO_nfct_full_psi( A)
MACRO_nfct_full_psi( T)





/** 
*  user routines
**/
void nfct_trafo( nfct_plan *ths)
{
  /**
   * use ths->my_fftw_r2r_plan
   *
   **/
  ths->g_hat = ths->g1;
  ths->g     = ths->g2;


  /** 
   * form \f$ \hat g_k = \frac{\hat f_k}{c_k\left(\phi\right)} \text{ for }
   * k \in I_N \f$
   *
   **/ 
  TIC(0)
  nfct_D_A( ths);
  TOC(0)


  /** 
   * compute by d-variate discrete Fourier transform
   * \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   * \text{ for } l \in I_n \f$
   *
   **/
  TIC(1)
  fftw_execute( ths->my_fftw_r2r_plan);
  TOC(1)


  if( ths->nfct_flags & PRE_FULL_PSI)
    nfct_full_psi__A( ths);


  /**
   *  set \f$ f_j = \sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
   *  \text{ for } j=0,\hdots,M-1 \f$
   *
   **/
  TIC(2)
  nfct_B_A( ths);
  TOC(2)

  if(ths->nfct_flags & PRE_FULL_PSI) {
    free( ths->psi_index_g);
    free( ths->psi_index_f);
  }

} /* nfct_trafo */




void nfct_adjoint( nfct_plan *ths)
{
  /**
   * use ths->my_fftw_plan 
   *
   **/
  ths->g_hat = ths->g2;
  ths->g     = ths->g1;

  if( ths->nfct_flags & PRE_FULL_PSI)
    nfct_full_psi__T( ths);

  /** 
   * set \f$ g_l = \sum_{j=0}^{M-1} f_j \psi\left(x_j-\frac{l}{n}\right)
   * \text{ for } l \in I_n,m(x_j) \f$
   *
   **/
  TIC(2)
  nfct_B_T( ths);
  TOC(2)

  if(ths->nfct_flags & PRE_FULL_PSI) {
    free( ths->psi_index_g);
    free( ths->psi_index_f);
  }

  /**
   * compute by d-variate discrete cosine transform
   * \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   * \text{ for }  k \in I_N\f$
   *
   **/ 
  TIC(1)
  fftw_execute( ths->my_fftw_r2r_plan);
  TOC(1)

  /**
   * form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
   * k \in I_N \f$
   *
   **/
  TIC(0)
  nfct_D_T( ths);
  TOC(0)

} /* nfct_adjoint */



/**
* initialization of direct transform
*
**/
void nfct_precompute_phi_hut( nfct_plan *ths)
{
  int kg[ths->d];                       /**< index over all frequencies       */
  int t;                                /**< index over all dimensions        */

  ths->c_phi_inv = (double**)fftw_malloc( ths->d * sizeof( double*));

  for( t = 0; t < ths->d; t++)
  {
    ths->c_phi_inv[t] = (double*)fftw_malloc( ths->N[t] * sizeof( double));

    for( kg[t] = 0; kg[t] < ths->N[t]; kg[t]++)
    {
      ths->c_phi_inv[t][kg[t]] = MACRO_compute_PHI_HUT_INV;
    }
  }
} /* nfct_phi_hut */



void nfct_precompute_psi( nfct_plan *ths)
{
  int t;                                /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int lc[ths->d];                       /**< index 0<=lj<u+o+1                */
  int lb[ths->d];                       /**< depends on x_j                   */

  for (t = 0; t < ths->d; t++)
  {
    for(j = 0; j < ths->M_total; j++)
    {

      MACRO_nfct__lower_boundary( j, t);

      for( lc[t] = 0; lc[t] < NFCT_SUMMANDS; lc[t]++)
	ths->psi[(j * ths->d + t) * NFCT_SUMMANDS + lc[t]] = MACRO_compute_PSI;

    } /* for(j) */
  } /* for(t) */
} /* nfct_precompute_psi */







void nfct_init_help( nfct_plan *ths)
{
  int t;                                /**< index over all dimensions        */

  ths->N_total = nfct_prod_int( ths->N, ths->d);

  ths->sigma   = (double*)fftw_malloc( ths->d * sizeof(double));
  
  for( t = 0; t < ths->d; t++)
    ths->sigma[t] = ((double)( ths->n[t] - 1)) / ths->N[t];

  /**
   * assign r2r transform kinds for each dimension
   **/
  ths->r2r_kind = (fftw_r2r_kind*)fftw_malloc ( ths->d * sizeof (fftw_r2r_kind));
  for (t = 0; t < ths->d; t++)
    ths->r2r_kind[t] = FFTW_REDFT00;


  NFCT_WINDOW_HELP_INIT;

  if(ths->nfct_flags & MALLOC_X)
    ths->x = (double*)fftw_malloc( ths->d * ths->M_total * sizeof( double));

  if(ths->nfct_flags & MALLOC_F_HAT)
    ths->f_hat = (double*)fftw_malloc( ths->N_total * sizeof( double));

  if(ths->nfct_flags & MALLOC_F)
    ths->f = (double*)fftw_malloc( ths->M_total * sizeof( double));

  if(ths->nfct_flags & PRE_PHI_HUT)
    nfct_precompute_phi_hut(ths);

  /* NO FFTW_MALLOC HERE */
  if(ths->nfct_flags & PRE_PSI)
  {
    ths->psi =
      (double*) malloc(ths->M_total * ths->d * NFCT_SUMMANDS * sizeof( double));

    /**
     * set default for full_psi_eps 
     **/
    ths->nfct_full_psi_eps = pow(10, -10);
  }
  
  if(ths->nfct_flags & FFTW_INIT)
  {
    ths->g1 =
      (double*)fftw_malloc( nfct_prod_int(ths->n, ths->d) * sizeof( double));

    if(ths->nfct_flags & FFT_OUT_OF_PLACE)
      ths->g2 =
	(double*) fftw_malloc( nfct_prod_int(ths->n, ths->d) * sizeof(double));
    else
      ths->g2 = ths->g1;
    
    ths->my_fftw_r2r_plan = 
      fftw_plan_r2r( ths->d, ths->n, ths->g1, ths->g2, 
		     ths->r2r_kind, ths->fftw_flags);
  }
}





void nfct_init( nfct_plan *ths, int d, int *N, int M_total)
{
  int t;

  ths->d       = d;
  ths->M_total = M_total;

  ths->N = (int*) fftw_malloc( ths->d * sizeof( int));

  for(t = 0;t < d; t++)
    ths->N[t] = N[t];

  ths->n = (int*) fftw_malloc( ths->d * sizeof( int));

  for( t = 0; t < d; t++)
    ths->n[t] = nfct_fftw_2N( nfft_next_power_of_2( ths->N[t]));

  WINDOW_HELP_ESTIMATE_m;

  ths->nfct_flags = NFCT_DEFAULT_FLAGS;
  ths->fftw_flags = FFTW_DEFAULT_FLAGS;

  nfct_init_help( ths);    
}


void nfct_init_m( nfct_plan *ths, int d, int *N, int M_total, int m)
{
  int t, n[d]; 
 
  for( t = 0; t < d; t++)  
    n[t] = nfct_fftw_2N( nfft_next_power_of_2( N[t]));  
 
  nfct_init_guru( ths, d, N, M_total, n, m, NFCT_DEFAULT_FLAGS, FFTW_DEFAULT_FLAGS); 
} 


void nfct_init_guru( nfct_plan *ths, int d, int *N,
		     int M_total, int *n, int m,
		     unsigned nfct_flags, unsigned fftw_flags)
{
  int t;             /**< index over all dimensions */

  ths->d = d;
  ths->M_total = M_total;

  ths->N = (int*)fftw_malloc( ths->d * sizeof( int));

  for( t = 0; t < d; t++)
    ths->N[t] = N[t];

  ths->n = (int*)fftw_malloc( ths->d * sizeof( int));

  for( t = 0; t < d; t++)
    ths->n[t] = n[t];

  ths->m = m;

  ths->nfct_flags = nfct_flags;
  ths->fftw_flags = fftw_flags;

  nfct_init_help( ths);  
}


void nfct_init_1d( nfct_plan *ths, int N0, int M_total)
{
  int N[1];

  N[0] = N0;
  nfct_init( ths, 1, N, M_total);
}

void nfct_init_2d( nfct_plan *ths, int N0, int N1, int M_total)
{
  int N[2];

  N[0] = N0;
  N[1] = N1;
  nfct_init( ths, 2, N, M_total);
}

void nfct_init_3d( nfct_plan *ths, int N0, int N1, int N2, int M_total)
{
  int N[3];

  N[0] = N0;
  N[1] = N1;
  N[2] = N2;
  nfct_init( ths, 3, N, M_total);
}


/**
 *  deallocate the memory of the plan
 **/
void nfct_finalize( nfct_plan *ths)
{
  /**
   * index over dimensions
   **/
  int t;

  if( ths->nfct_flags & FFTW_INIT)
  {
    fftw_destroy_plan( ths->my_fftw_r2r_plan);
     
    if(ths->nfct_flags & FFT_OUT_OF_PLACE)
      fftw_free( ths->g2);
      
    fftw_free( ths->g1);
  }

  /* NO FFTW_FREE HERE */
  if( ths->nfct_flags & PRE_PSI)
  {
    free( ths->psi);
  }

  if(ths->nfct_flags & PRE_PHI_HUT)
  {
    for( t = 0; t < ths->d; t++)
      fftw_free( ths->c_phi_inv[t]);
    fftw_free( ths->c_phi_inv);
  }

  if( ths->nfct_flags & MALLOC_F)
    fftw_free( ths->f);

  if( ths->nfct_flags & MALLOC_F_HAT)
    fftw_free( ths->f_hat);

  if( ths->nfct_flags & MALLOC_X)
  fftw_free( ths->x);
 
  WINDOW_HELP_FINALIZE;
 
  fftw_free( ths->N);
  fftw_free( ths->n);
  fftw_free( ths->sigma);

} /* nfct_finalize */

