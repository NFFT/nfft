/**
 * \file fpt.c
 * \brief Implementation file for the FPT module
 * \author Jens Keiner
 */

/* Include standard C headers. */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/* Include FPT module header. */
#include "nfft3.h"

/* Include NFFT3 utilities header. */
#include "util.h"

/* Some macros for index calculation. */

/** Computes the minimum degree at the top of a cascade. */
#define K_START_TILDE(x,y) (NFFT_MAX(NFFT_MIN(x,y-2),0))
/** Computes the maximum degree at the top of a cascade. */
#define K_END_TILDE(x,y) NFFT_MIN(x,y-1)
/**
 * Computes the index of the first block of four functions in a cascade
 * level.
 */
#define FIRST_L(x,y) ((int)floor((x)/(double)y))
/**
 * Computes the index of the last block of four functions in a cascade
 * level.
 */
#define LAST_L(x,y) ((int)ceil(((x)+1)/(double)y)-1)

#define N_TILDE(y) (y-1)

#define IS_SYMMETRIC(x,y,z) (x >= ((y-1.0)/z))
//#define IS_SYMMETRIC(x,y,z) (x >= ((y-1.0)/z))

#ifdef TEST_STAB
  #define NFFT_MAX(a,b) ((a>b)?(a):(b))
#endif

#define FPT_BREAK_EVEN 4

/**
 * Holds data for a single multiplication step in the cascade summation.
 */
typedef struct fpt_step_
{
  bool stable;                            /**< Indicates if the values
                                               contained represent a fast or
                                               a slow stabilized step.        */
  int N_stab;                             /**< TODO Add comment here.         */
  int t_stab;                             /**< TODO Add comment here.         */
  double **a11,**a12,**a21,**a22;         /**< The matrix components          */
  double *gamma;                          /**<                                */
} fpt_step;

/**
 * Holds data for a single cascade summation.
 */
typedef struct fpt_data_
{
  fpt_step **steps;                       /**< The cascade summation steps   */
  int k_start;
  double *alphaN;
  double *betaN;
  double *gammaN;
  double alpha_0;
  double beta_0;
  double gamma_m1;
  /* Data for direct transform. */
  double *alpha;
  double *beta;
  double *gamma;
} fpt_data;

/**
 * Holds data for a set of cascade summations.
 */
typedef struct fpt_set_s_
{
  int flags;                              /**< The flags                     */
  int M;                                  /**< The number of DPT transforms  */
  int N;                                  /**< The transform length. Must be
                                               a power of two.               */
  int t;                                  /**< The exponent of N             */
  fpt_data *dpt;                          /**< The DPT transform data        */
  double **xcvecs;                        /**< Array of pointers to arrays
                                               containing the Chebyshev
                                               nodes                         */
  double *xc;                             /**< Array for Chebychev-nodes.    */
  double complex *temp;                          /**< */
  double complex *work;                          /**< */
  double complex *result;                        /**< */
  double complex *vec3;
  double complex *vec4;
  double complex *z;
  fftw_plan *plans_dct3;                  /**< Transform plans for the fftw
                                               library                       */
  fftw_plan *plans_dct2;                  /**< Transform plans for the fftw
                                               library                       */
  fftw_r2r_kind *kinds;                   /**< Transform kinds for fftw
                                               library                       */
  fftw_r2r_kind *kindsr;                  /**< Transform kinds for fftw
                                               library                       */

  int *lengths; /**< Transform lengths for fftw library */

  /* Data for slow transforms. */
  double *xc_slow;
} fpt_set_s;

inline void abuvxpwy(double a, double b, double complex* u, double complex* x, double* v,
  double complex* y, double* w, int n)
{
  int l;
  double complex *u_ptr, *x_ptr, *y_ptr;
  double *v_ptr, *w_ptr;

  u_ptr = u;
  x_ptr = x;
  v_ptr = v;
  y_ptr = y;
  w_ptr = w;

  for (l = 0; l < n; l++)
  {
    *u++ = a * (b * (*v++) * (*x++) + (*w++) * (*y++));
  }
}

#define ABUVXPWY_SYMMETRIC(NAME,S1,S2) \
inline void NAME(double a, double b, double complex* u, double complex* x, double* v, \
  double complex* y, double* w, int n) \
{ \
  int l; \
  double complex *u_ptr, *x_ptr, *y_ptr; \
  double *v_ptr, *w_ptr; \
  \
  u_ptr = u; \
  x_ptr = x; \
  v_ptr = v; \
  y_ptr = y; \
  w_ptr = w; \
  \
  for (l = 0; l < n/2; l++) \
  { \
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++)); \
  } \
  v_ptr--; \
  w_ptr--; \
  for (l = 0; l < n/2; l++) \
  { \
    *u_ptr++ = a * (b * S1 * (*v_ptr--) * (*x_ptr++) + S2 * (*w_ptr--) * (*y_ptr++)); \
  } \
}

ABUVXPWY_SYMMETRIC(abuvxpwy_symmetric1,1.0,-1.0)
ABUVXPWY_SYMMETRIC(abuvxpwy_symmetric2,-1.0,1.0)

#define ABUVXPWY_SYMMETRIC_1(NAME,S1) \
inline void NAME(double a, double b, double complex* u, double complex* x, double* v, \
  double complex* y, double* w, int n) \
{ \
  int l; \
  double complex *u_ptr, *x_ptr, *y_ptr; \
  double *v_ptr, *w_ptr; \
  \
  u_ptr = u; \
  x_ptr = x; \
  v_ptr = v; \
  y_ptr = y; \
  w_ptr = w; \
  \
  for (l = 0; l < n/2; l++) \
  { \
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++)); \
  } \
  v_ptr--; \
  /*w_ptr--;*/ \
  for (l = 0; l < n/2; l++) \
  { \
    *u_ptr++ = a * (b * S1 * (*v_ptr--) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++)); \
  } \
}

ABUVXPWY_SYMMETRIC_1(abuvxpwy_symmetric1_1,1.0)
ABUVXPWY_SYMMETRIC_1(abuvxpwy_symmetric1_2,-1.0)

#define ABUVXPWY_SYMMETRIC_2(NAME,S1) \
inline void NAME(double a, double b, double complex* u, double complex* x, double* v, \
  double complex* y, double* w, int n) \
{ \
  int l; \
  double complex *u_ptr, *x_ptr, *y_ptr; \
  double *v_ptr, *w_ptr; \
  \
  u_ptr = u; \
  x_ptr = x; \
  v_ptr = v; \
  y_ptr = y; \
  w_ptr = w; \
  \
  for (l = 0; l < n/2; l++) \
  { \
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++)); \
  } \
  /*v_ptr--;*/ \
  w_ptr--; \
  for (l = 0; l < n/2; l++) \
  { \
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + S1 * (*w_ptr--) * (*y_ptr++)); \
  } \
}

ABUVXPWY_SYMMETRIC_2(abuvxpwy_symmetric2_1,1.0)
ABUVXPWY_SYMMETRIC_2(abuvxpwy_symmetric2_2,-1.0)

inline void auvxpwy(double a, double complex* u, double complex* x, double* v, double complex* y,
  double* w, int n)
{
  int l;
  double complex *u_ptr, *x_ptr, *y_ptr;
  double *v_ptr, *w_ptr;

  u_ptr = u;
  x_ptr = x;
  v_ptr = v;
  y_ptr = y;
  w_ptr = w;

  for (l = 0; l < n; l++)
  {
    /*fprintf(stderr,"u = %le, v = %le, x = %le, w = %le, y = %le\n",*u_ptr,*v_ptr,*x_ptr,*w_ptr,*y_ptr);*/
    *u_ptr++ = a * ((*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++));
  }
  /*fprintf(stderr,"\n");*/
}

#define AUVXPWY_SYMMETRIC(NAME,S1,S2) \
inline void NAME(double a, double complex* u, double complex* x, double* v, double complex* y, \
  double* w, int n) \
{ \
  int l; \
  double complex *u_ptr, *x_ptr, *y_ptr; \
  double *v_ptr, *w_ptr; \
\
  u_ptr = u; \
  x_ptr = x; \
  v_ptr = v; \
  y_ptr = y; \
  w_ptr = w; \
\
  /*for (l = 0; l < n; l++)*/ \
  /*{*/ \
    /*fprintf(stderr,"u = %le, v = %le, x = %le, w = %le, y = %le\n",*/ \
    /*  u_ptr[l],v_ptr[l],x_ptr[l],w_ptr[l],y_ptr[l]);*/ \
  /*}*/ \
  \
  \
  for (l = 0; l < n/2; l++) \
  { \
    /*fprintf(stderr,"u = %le, v = %le, x = %le, w = %le, y = %le\n",*u,*v,*x,*w,*y);*/ \
    *u_ptr++ = a * ((*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++)); \
  } \
  v_ptr--; \
  w_ptr--; \
  for (l = 0; l < n/2; l++) \
  { \
    /* fprintf(stderr,"u = %le, v = %le, x = %le, w = %le, y = %le\n",*u,*v,*x,*w,*y);*/ \
    *u_ptr++ = a * (S1 * (*v_ptr--) * (*x_ptr++) + S2 * (*w_ptr--) * (*y_ptr++)); \
  } \
  /*fprintf(stderr,"\n");*/ \
}

AUVXPWY_SYMMETRIC(auvxpwy_symmetric,1.0,-1.0)

#define FPT_DO_STEP(NAME,M1_FUNCTION,M2_FUNCTION) \
inline void NAME(double complex  *a, double complex *b, double *a11, double *a12, \
  double *a21, double *a22, double gamma, int tau, fpt_set set) \
{ \
  /** The length of the coefficient arrays. */ \
  int length = 1<<(tau+1); \
  /** Twice the length of the coefficient arrays. */ \
  double norm = 1.0/(length<<1); \
  \
  /* Compensate for factors introduced by a raw DCT-III. */ \
  a[0] *= 2.0; \
  b[0] *= 2.0; \
  \
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */ \
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)a,(double*)a); \
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)b,(double*)b); \
  \
  /*for (k = 0; k < length; k++)*/ \
  /*{*/ \
    /*fprintf(stderr,"fpt_do_step: a11 = %le, a12 = %le, a21 = %le, a22 = %le\n",*/ \
    /*  a11[k],a12[k],a21[k],a22[k]);*/ \
  /*}*/ \
  \
  /* Check, if gamma is zero. */ \
  if (gamma == 0.0) \
  { \
    /*fprintf(stderr,"gamma == 0!\n");*/ \
    /* Perform multiplication only for second row. */ \
    M2_FUNCTION(norm,b,b,a22,a,a21,length); \
  } \
  else \
  { \
    /*fprintf(stderr,"gamma != 0!\n");*/ \
    /* Perform multiplication for both rows. */ \
    M2_FUNCTION(norm,set->z,b,a22,a,a21,length); \
    M1_FUNCTION(norm*gamma,a,a,a11,b,a12,length); \
    memcpy(b,set->z,length*sizeof(double complex)); \
    /* Compute Chebyshev-coefficients using a DCT-II. */ \
    fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a); \
    /* Compensate for factors introduced by a raw DCT-II. */ \
    a[0] *= 0.5; \
  } \
  \
  /* Compute Chebyshev-coefficients using a DCT-II. */ \
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b); \
  /* Compensate for factors introduced by a raw DCT-II. */ \
  b[0] *= 0.5; \
}

FPT_DO_STEP(fpt_do_step,auvxpwy,auvxpwy)
FPT_DO_STEP(fpt_do_step_symmetric,auvxpwy_symmetric,auvxpwy_symmetric)
FPT_DO_STEP(fpt_do_step_symmetric_u,auvxpwy_symmetric,auvxpwy)
FPT_DO_STEP(fpt_do_step_symmetric_l,auvxpwy,auvxpwy_symmetric)

#define FPT_DO_STEP_TRANSPOSED(NAME,M1_FUNCTION,M2_FUNCTION) \
inline void NAME(double complex  *a, double complex *b, double *a11, \
  double *a12, double *a21, double *a22, double gamma, int tau, fpt_set set) \
{ \
  /** The length of the coefficient arrays. */ \
  int length = 1<<(tau+1); \
  /** Twice the length of the coefficient arrays. */ \
  double norm = 1.0/(length<<1); \
  \
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */ \
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)a,(double*)a); \
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)b,(double*)b); \
  \
  /* Perform matrix multiplication. */ \
  M1_FUNCTION(norm,gamma,set->z,a,a11,b,a21,length); \
  M2_FUNCTION(norm,gamma,b,a,a12,b,a22,length); \
  memcpy(a,set->z,length*sizeof(double complex)); \
  \
  /* Compute Chebyshev-coefficients using a DCT-II. */ \
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a); \
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b); \
}

FPT_DO_STEP_TRANSPOSED(fpt_do_step_transposed,abuvxpwy,abuvxpwy)
FPT_DO_STEP_TRANSPOSED(fpt_do_step_transposed_symmetric,abuvxpwy_symmetric1,abuvxpwy_symmetric2)
FPT_DO_STEP_TRANSPOSED(fpt_do_step_transposed_symmetric_u,abuvxpwy_symmetric1_1,abuvxpwy_symmetric1_2)
FPT_DO_STEP_TRANSPOSED(fpt_do_step_transposed_symmetric_l,abuvxpwy_symmetric2_2,abuvxpwy_symmetric2_1)

void eval_clenshaw(const double *x, double *y, int size, int k, const double *alpha,
  const double *beta, const double *gamma)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  double a,b,x_val_act,a_old;
  const double *x_act;
  double *y_act;
  const double *alpha_act, *beta_act, *gamma_act;

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;

    if (k == 0)
    {
      *y_act = 1.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gamma[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
         b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
    }
    x_act++;
    y_act++;
  }
}

int eval_clenshaw_thresh(const double *x, double *y, int size, int k, 
  const double *alpha, const double *beta, const double *gamma, const 
  double threshold)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  double a,b,x_val_act,a_old;
  const double *x_act;
  double *y_act;
  const double *alpha_act, *beta_act, *gamma_act;

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;

    if (k == 0)
    {
     *y_act = 1.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gamma[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
         b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
      if (fabs(*y_act) > threshold)
      {
        return 1;
      }
    }
    x_act++;
    y_act++;
  }
  return 0;
}

/**
 * Clenshaw algorithm
 *
 * Evaluates a sum of real-valued functions \f$P_k : \mathbb{R} \rightarrow
 * \mathbb{R}\f$
 * \f[
 *   f(x) = \sum_{k=0}^N a_k P_k(x) \quad (N \in \mathbb{N}_0)
 * \f]
 * obeying a three-term recurrence relation
 * \f[
 *   P_{k+1}(x) = (alpha_k * x + beta_k)*P_{k}(x) + gamma_k P_{k-1}(x) \quad
 *   (alpha_k, beta_k, gamma_k \in \mathbb{R},\; k \ge 0)
 * \f]
 * with initial conditions \f$P_{-1}(x) := 0\f$, \f$P_0(x) := \lambda\f$
 * for given double complex coefficients \f$\left(a_k\right)_{k=0}^N \in
 * \mathbb{C}^{N+1}\f$ at given nodes \f$\left(x_j\right)_{j=0}^M \in
 * \mathbb{R}^{M+1}\f$, \f$M \in \mathbb{N}_0\f$.
 */
void eval_sum_clenshaw(int N, int M, double complex* a, double *x, double complex *y,
  double complex *temp, double *alpha, double *beta, double *gamma, double lambda)
{
  int j,k;
  double complex* it1 = temp;
  double complex* it2 = y;
  double complex aux;

  /* Clenshaw's algorithm */
  for (j = 0; j <= M; j++)
  {
    it2[j] = a[N];
  }

  if (N > 0)
  {
    for (j = 0; j <= M; j++)
    {
      it1[j] = a[N-1];
    }

    //fprintf(stdout,"N = %d\n",N);
    for (k = N; k > 1; k--)
    {

      for (j = 0; j <= M; j++)
      {
        aux = a[k-2] + it2[j] * gamma[k-1];
        it2[j] = it1[j] + it2[j] * (alpha[k-1] * x[j] + beta[k-1]);
        it1[j] = aux;
      }
    }


    /* Compute final step. */
    for (j = 0; j <= M; j++)
    {
      it2[j] = it1[j] + it2[j] * (alpha[0] * x[j] + beta[0]);
    }
  }

  /* Compute final result by multiplying with the constant lambda */
  for (j = 0; j <= M; j++)
  {
    y[j] = lambda * it2[j];
  }
}

/**
 * Clenshaw algorithm
 *
 * Evaluates a sum of real-valued functions \f$P_k : \mathbb{R} \rightarrow
 * \mathbb{R}\f$
 * \f[
 *   f(x) = \sum_{k=0}^N a_k P_k(x) \quad (N \in \mathbb{N}_0)
 * \f]
 * obeying a three-term recurrence relation
 * \f[
 *   P_{k+1}(x) = (alpha_k * x + beta_k)*P_{k}(x) + gamma_k P_{k-1}(x) \quad
 *   (alpha_k, beta_k, gamma_k \in \mathbb{R},\; k \ge 0)
 * \f]
 * with initial conditions \f$P_{-1}(x) := 0\f$, \f$P_0(x) := \lambda\f$
 * for given double complex coefficients \f$\left(a_k\right)_{k=0}^N \in
 * \mathbb{C}^{N+1}\f$ at given nodes \f$\left(x_j\right)_{j=0}^M \in
 * \mathbb{R}^{M+1}\f$, \f$M \in \mathbb{N}_0\f$.
 */
void eval_sum_clenshaw_transposed(int N, int M, double complex* a, double *x,
  double complex *y, double complex *temp, double *alpha, double *beta, double *gamma,
  double lambda)
{
  int j,k;
  double complex* it1 = temp;
  double complex* it2 = y;
  double complex aux;

  /* Compute final result by multiplying with the constant lambda */
  a[0] = 0.0;
  for (j = 0; j <= M; j++)
  {
    it2[j] = lambda * y[j];
    a[0] += it2[j];
  }

  if (N > 0)
  {
    /* Compute final step. */
    a[1] = 0.0;
    for (j = 0; j <= M; j++)
    {
      it1[j] = it2[j];
      it2[j] = it2[j] * (alpha[0] * x[j] + beta[0]);
      a[1] += it2[j];
    }

    for (k = 2; k <= N; k++)
    {
      a[k] = 0.0;
      for (j = 0; j <= M; j++)
      {
        aux = it1[j];
        it1[j] = it2[j];
        it2[j] = it2[j]*(alpha[k-1] * x[j] + beta[k-1]) + gamma[k-1] * aux;
        a[k] += it2[j];
      }
    }
  }
}


fpt_set fpt_init(const int M, const int t, const unsigned int flags)
{
  /** Polynomial length */
  int plength;
  /** Cascade level */
  int tau;
  /** Index m */
  int m;
  int k;

  /* Allocate memory for new DPT set. */
  fpt_set_s *set = (fpt_set_s*)malloc(sizeof(fpt_set_s));

  /* Save parameters in structure. */
  set->flags = flags;

  //fprintf(stderr,"\nfpt_init: flags = %d \t %d\n",set->flags,flags);

  set->M = M;
  set->t = t;
  set->N = 1<<t;

  /* Allocate memory for L transforms. */
  set->dpt = (fpt_data*) malloc((M+1)*sizeof(fpt_data));

  /* Initialize with NULL pointer. */
  for (m = 0; m <= set->M; m++)
  {
    set->dpt[m].steps = (fpt_step**) NULL;
  }

  /* Create arrays with Chebyshev nodes. */

  /* Initialize array with Chebyshev coefficients for the polynomial x. This
   * would be trivially an array containing a 1 as second entry with all other
   * coefficients set to zero. In order to compensate for the multiplicative
   * factor 2 introduced by the DCT-III, we set this coefficient to 0.5 here. */

  /* Allocate memory for array of pointers to node arrays. */
  set->xcvecs = (double**) malloc((set->t/*-1*/)*sizeof(double*));
  /* For each polynomial length starting with 4, compute the Chebyshev nodes
   * using a DCT-III. */
  plength = 4;
  for (tau = 1; tau < /*t*/t+1; tau++)
  {
    /* Allocate memory for current array. */
    set->xcvecs[tau-1] = (double*) malloc(plength*sizeof(double));
    for (k = 0; k < plength; k++)
    {
      set->xcvecs[tau-1][k] = cos(((k+0.5)*PI)/plength);
    }
    plength = plength << 1;
  }

  /** Allocate memory for auxilliary arrays. */
  set->work = (double complex*) malloc((2*set->N)*sizeof(double complex));
  set->result = (double complex*) malloc((2*set->N)*sizeof(double complex));

  /* Check if fast transform is activated. */
  if (set->flags & FPT_NO_FAST_ALGORITHM)
  {
  }
  else
  {
    /** Allocate memory for auxilliary arrays. */
    set->vec3 = (double complex*) malloc(set->N*sizeof(double complex));
    set->vec4 = (double complex*) malloc(set->N*sizeof(double complex));
    set->z = (double complex*) malloc(set->N*sizeof(double complex));

    /** Initialize FFTW plans. */
    set->plans_dct3 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(set->t/*-1*/));
    set->plans_dct2 = (fftw_plan*) fftw_malloc(sizeof(fftw_plan)*(set->t/*-1*/));
    set->kinds      = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    set->kinds[0]   = FFTW_REDFT01;
    set->kinds[1]   = FFTW_REDFT01;
    set->kindsr     = (fftw_r2r_kind*) malloc(2*sizeof(fftw_r2r_kind));
    set->kindsr[0]  = FFTW_REDFT10;
    set->kindsr[1]  = FFTW_REDFT10;
    set->lengths    = (int*) malloc((set->t/*-1*/)*sizeof(int));
    for (tau = 0, plength = 4; tau < set->t/*-1*/; tau++, plength<<=1)
    {
      set->lengths[tau] = plength;
      set->plans_dct3[tau] =
        fftw_plan_many_r2r(1, &set->lengths[tau], 2, (double*)set->work, NULL,
                           2, 1, (double*)set->result, NULL, 2, 1, set->kinds,
                           0);
      set->plans_dct2[tau] =
        fftw_plan_many_r2r(1, &set->lengths[tau], 2, (double*)set->work, NULL,
                           2, 1, (double*)set->result, NULL, 2, 1,set->kindsr,
                           0);
    }
    free(set->lengths);
    free(set->kinds);
    free(set->kindsr);
    set->lengths = NULL;
    set->kinds = NULL;
    set->kindsr = NULL;
  }

  if (set->flags & FPT_NO_DIRECT_ALGORITHM)
  {
  }
  else
  {
    set->xc_slow = (double*) malloc((set->N+1)*sizeof(double));
    set->temp = (double complex*) calloc((set->N+1),sizeof(double complex));
  }

  /* Return the newly created DPT set. */
  return set;
}

void fpt_precompute(fpt_set set, const int m, const double *alpha,
                    const double *beta, const double *gamma, int k_start,
                    const double threshold)
{

  int tau;          /**< Cascade level                                       */
  int l;            /**< Level index                                         */
  int plength;      /**< Length of polynomials for the next level in the
                         cascade                                             */
  int degree;       /**< Degree of polynomials for the current level in the
                         cascade                                             */
  int firstl;       /**< First index l for current cascade level             */
  int lastl;        /**< Last index l for current cascade level and current  */
  int plength_stab; /**< Length of polynomials for the next level in the
                         cascade for stabilization                           */
  int degree_stab;  /**< Degree of polynomials for the current level in the
                         cascade for stabilization                           */
  double *a11;      /**< Array containing function values of the
                         (1,1)-component of U_k^n.                           */
  double *a12;      /**< Array containing function values of the
                         (1,2)-component of U_k^n.                           */
  double *a21;      /**< Array containing function values of the
                         (2,1)-component of U_k^n.                           */
  double *a22;      /**< Array containing function values of the
                         (2,2)-component of U_k^n.                           */
  const double *calpha;
  const double *cbeta;
  const double *cgamma;
  int needstab = 0; /**< Used to indicate that stabilization is neccessary.  */
  int k_start_tilde;
  int N_tilde;
  int clength;
  int clength_1;
  int clength_2;
  int t_stab, N_stab;

  fpt_data *data;

  /* Allocate memory for DPT transform data. */
  //set->dpt[m] = (fpt_data*) malloc(sizeof(fpt_data));

  /* Get pointer to DPT data. */
  data = &(set->dpt[m]);

  /* Check, if already precomputed. */
  if (data->steps != NULL)
  {
    return;
  }

  /* Save k_start. */
  data->k_start = k_start;

  /* Check if fast transform is activated. */
  if (set->flags & FPT_NO_FAST_ALGORITHM)
  {
  }
  else
  {
    /* Save recursion coefficients. */
    data->alphaN = (double*) malloc((set->t-1)*sizeof(double complex));
    data->betaN = (double*) malloc((set->t-1)*sizeof(double complex));
    data->gammaN = (double*) malloc((set->t-1)*sizeof(double complex));

    for (tau = 2; tau <= set->t; tau++)
    {

      data->alphaN[tau-2] = alpha[1<<tau];
      data->betaN[tau-2] = beta[1<<tau];
      data->gammaN[tau-2] = gamma[1<<tau];
    }

    data->alpha_0 = alpha[1];
    data->beta_0 = beta[1];
    data->gamma_m1 = gamma[0];

    k_start_tilde = K_START_TILDE(data->k_start,nfft_next_power_of_2(data->k_start)
      /*set->N*/);
    N_tilde = N_TILDE(set->N);

    /* Allocate memory for the cascade with t = log_2(N) many levels. */
    data->steps = (fpt_step**) malloc(sizeof(fpt_step*)*set->t);

    /* For tau = 1,...t compute the matrices U_{n,tau,l}. */
    plength = 4;
    for (tau = 1; tau < set->t; tau++)
    {
      /* Compute auxilliary values. */
      degree = plength>>1;
      /* Compute first l. */
      firstl = FIRST_L(k_start_tilde,plength);
      /* Compute last l. */
      lastl = LAST_L(N_tilde,plength);

      /* Allocate memory for current level. This level will contain 2^{t-tau-1}
       * many matrices. */
      data->steps[tau] = (fpt_step*) fftw_malloc(sizeof(fpt_step)
                         * (lastl+1));

      /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
      for (l = firstl; l <= lastl; l++)
      {
        if (set->flags & FPT_AL_SYMMETRY && IS_SYMMETRIC(l,m,plength))
        {
          //fprintf(stderr,"fpt_precompute(%d): symmetric step\n",m);
          //fflush(stderr);
          clength = plength/2;
        }
        else
        {
          clength = plength;
        }

        /* Allocate memory for the components of U_{n,tau,l}. */
        a11 = (double*) fftw_malloc(sizeof(double)*clength);
        a12 = (double*) fftw_malloc(sizeof(double)*clength);
        a21 = (double*) fftw_malloc(sizeof(double)*clength);
        a22 = (double*) fftw_malloc(sizeof(double)*clength);

        /* Evaluate the associated polynomials at the 2^{tau+1} Chebyshev
         * nodes. */

        /* Get the pointers to the three-term recurrence coeffcients. */
        calpha = &(alpha[plength*l+1+1]);
        cbeta = &(beta[plength*l+1+1]);
        cgamma = &(gamma[plength*l+1+1]);

        if (set->flags & FPT_NO_STABILIZATION)
        {
          /* Evaluate P_{2^{tau}-2}^n(\cdot,2^{tau+1}l+2). */
          eval_clenshaw(set->xcvecs[tau-1], a11, clength, degree-2, calpha, cbeta,
            cgamma);
          eval_clenshaw(set->xcvecs[tau-1], a12, clength, degree-1, calpha, cbeta,
            cgamma);
          calpha--;
          cbeta--;
          cgamma--;
          eval_clenshaw(set->xcvecs[tau-1], a21, clength, degree-1, calpha, cbeta,
            cgamma);
          eval_clenshaw(set->xcvecs[tau-1], a22, clength, degree, calpha, cbeta,
            cgamma);
          needstab = 0;
        }
        else
        {
          needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a11, clength, degree-2,
            calpha, cbeta, cgamma, threshold);
          if (needstab == 0)
          {
            /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+2). */
            needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a12, clength, degree-1,
              calpha, cbeta, cgamma, threshold);
            if (needstab == 0)
            {
              calpha--;
              cbeta--;
              cgamma--;
              /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1). */
              needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a21, clength,
                degree-1, calpha, cbeta, cgamma, threshold);
              if (needstab == 0)
              {
                /* Evaluate P_{2^{tau}}^n(\cdot,2^{tau+1}l+1). */
                needstab = eval_clenshaw_thresh(set->xcvecs[tau-1], a22, clength,
                  degree, calpha, cbeta, cgamma, threshold);
              }
            }
          }
        }

        /* Check if stabilization needed. */
        if (needstab == 0)
        {
          data->steps[tau][l].a11 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].a12 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].a21 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].a22 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].gamma = (double*) fftw_malloc(sizeof(double));
          /* No stabilization needed. */
          data->steps[tau][l].a11[0] = a11;
          data->steps[tau][l].a12[0] = a12;
          data->steps[tau][l].a21[0] = a21;
          data->steps[tau][l].a22[0] = a22;
          data->steps[tau][l].gamma[0] = gamma[plength*l+1+1];
          data->steps[tau][l].stable = true;
        }
        else
        {
          /* Stabilize. */
          degree_stab = degree*(2*l+1);
          nfft_next_power_of_2_exp((l+1)*(1<<(tau+1)),&N_stab,&t_stab);
          /*fprintf(stderr,"(l+1)*(1<<(tau+2)) = %d, N_stab = %d, t_stab = %d\n",
            (l+1)*(1<<(tau+2)),N_stab,t_stab);*/

          /* Old arrays are to small. */
          fftw_free(a11);
          fftw_free(a12);
          fftw_free(a21);
          fftw_free(a22);

          data->steps[tau][l].a11 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].a12 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].a21 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].a22 = (double**) fftw_malloc(sizeof(double*));
          data->steps[tau][l].gamma = (double*) fftw_malloc(sizeof(double));

          plength_stab = N_stab;

          if (m <= 1)
          {
            clength_1 = plength_stab/2;
            clength_2 = plength_stab/2;
          }
          else if (m%2 == 0)
          {
            clength_1 = plength_stab/2;
            clength_2 = plength_stab;
          }
          else
          {
            clength_1 = plength_stab;
            clength_2 = plength_stab/2;
          }


          /* Allocate memory for arrays. */
          a11 = (double*) fftw_malloc(sizeof(double)*clength_1);
          a12 = (double*) fftw_malloc(sizeof(double)*clength_1);
          a21 = (double*) fftw_malloc(sizeof(double)*clength_2);
          a22 = (double*) fftw_malloc(sizeof(double)*clength_2);

          /* Get the pointers to the three-term recurrence coeffcients. */
          calpha = &(alpha[2]);
          cbeta = &(beta[2]);
          cgamma = &(gamma[2]);
          /* Evaluate P_{2^{tau}(2l+1)-2}^n(\cdot,2). */
          eval_clenshaw(set->xcvecs[t_stab-2], a11, clength_1, degree_stab-2,
            calpha, cbeta, cgamma);
          /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,2). */
          eval_clenshaw(set->xcvecs[t_stab-2], a12, clength_1, degree_stab-1,
            calpha, cbeta, cgamma);
          calpha--;
          cbeta--;
          cgamma--;
          /* Evaluate P_{2^{tau}(2l+1)-1}^n(\cdot,1). */
          eval_clenshaw(set->xcvecs[t_stab-2], a21, clength_2, degree_stab-1,
            calpha, cbeta, cgamma);
          /* Evaluate P_{2^{tau}(2l+1)}^n(\cdot,1). */
          eval_clenshaw(set->xcvecs[t_stab-2], a22, clength_2, degree_stab+0,
            calpha, cbeta, cgamma);

          data->steps[tau][l].a11[0] = a11;
          data->steps[tau][l].a12[0] = a12;
          data->steps[tau][l].a21[0] = a21;
          data->steps[tau][l].a22[0] = a22;
          data->steps[tau][l].gamma[0] =  gamma[1+1];
          data->steps[tau][l].stable = false;
          data->steps[tau][l].t_stab = t_stab;
          data->steps[tau][l].N_stab = N_stab;
        }
      }
      /** Increase polynomial degree to next power of two. */
      plength = plength << 1;
    }
  }

  if (set->flags & FPT_NO_DIRECT_ALGORITHM)
  {
  }
  else
  {
    /* Check, if recurrence coefficients must be copied. */
    if (set->flags & FPT_PERSISTENT_DATA)
    {
      data->alpha = (double*) alpha;
      data->beta = (double*) beta;
      data->gamma = (double*) gamma;
    }
    else
    {
      data->alpha = (double*) malloc((set->N+1)*sizeof(double));
      data->beta = (double*) malloc((set->N+1)*sizeof(double));
      data->gamma = (double*) malloc((set->N+1)*sizeof(double));
      memcpy(data->alpha,alpha,(set->N+1)*sizeof(double));
      memcpy(data->beta,beta,(set->N+1)*sizeof(double));
      memcpy(data->gamma,gamma,(set->N+1)*sizeof(double));
    }
  }
}

void dpt_trafo(fpt_set set, const int m, const double complex *x, double complex *y,
  const int k_end, const unsigned int flags)
{
  int j;
  fpt_data *data = &(set->dpt[m]);
  int Nk;
  int tk;
  double norm;

  nfft_next_power_of_2_exp(k_end+1,&Nk,&tk);
  norm = 2.0/(Nk<<1);

  if (set->flags & FPT_NO_DIRECT_ALGORITHM)
  {
    return;
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
    /* Fill array with Chebyshev nodes. */
    for (j = 0; j <= k_end; j++)
    {
      set->xc_slow[j] = cos((PI*(j+0.5))/(k_end+1));
    }

    memset(set->result,0U,data->k_start*sizeof(double complex));
    memcpy(&set->result[data->k_start],x,(k_end-data->k_start+1)*sizeof(double complex));

    eval_sum_clenshaw(k_end, k_end, set->result, set->xc_slow,
      y, set->work, &data->alpha[1], &data->beta[1], &data->gamma[1],
      data->gamma_m1);
  }
  else
  {
    memset(set->temp,0U,data->k_start*sizeof(double complex));
    memcpy(&set->temp[data->k_start],x,(k_end-data->k_start+1)*sizeof(double complex));

    eval_sum_clenshaw(k_end, Nk-1, set->temp, set->xcvecs[tk-2],
      set->result, set->work, &data->alpha[1], &data->beta[1], &data->gamma[1],
      data->gamma_m1);

    fftw_execute_r2r(set->plans_dct2[tk-2],(double*)set->result,
      (double*)set->result);

    set->result[0] *= 0.5;
    for (j = 0; j < Nk; j++)
    {
      set->result[j] *= norm;
    }

    memcpy(y,set->result,(k_end+1)*sizeof(double complex));
  }
}

void fpt_trafo(fpt_set set, const int m, const double complex *x, double complex *y,
  const int k_end, const unsigned int flags)
{
  /* Get transformation data. */
  fpt_data *data = &(set->dpt[m]);
  /** */
  int Nk;
  /** */
  int tk;
  /** */
  int k_start_tilde;
  /** */
  int k_end_tilde;

  /** Level index \f$tau\f$ */
  int tau;
  /** Index of first block at current level */
  int firstl;
  /** Index of last block at current level */
  int lastl;
  /** Block index \f$l\f$ */
  int l;
  /** Length of polynomial coefficient arrays at next level */
  int plength;
  /** Polynomial array length for stabilization */
  int plength_stab;
  int t_stab;
  /** Current matrix \f$U_{n,tau,l}\f$ */
  fpt_step *step;
  /** */
  fftw_plan plan;
  int length = k_end+1;
  fftw_r2r_kind kinds[2] = {FFTW_REDFT01,FFTW_REDFT01};

  /** Loop counter */
  int k;

  double complex *work_ptr;
  const double complex *x_ptr;
  
  /* Check, if slow transformation should be used due to small bandwidth. */
  if (k_end < FPT_BREAK_EVEN)
  {
    /* Use NDSFT. */
    dpt_trafo(set, m, x, y, k_end, flags);
  }
  
  nfft_next_power_of_2_exp(k_end,&Nk,&tk);
  k_start_tilde = K_START_TILDE(data->k_start,Nk);
  k_end_tilde = K_END_TILDE(k_end,Nk);

  /* Check if fast transform is activated. */
  if (set->flags & FPT_NO_FAST_ALGORITHM)
  {
    return;
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
    plan = fftw_plan_many_r2r(1, &length, 2, (double*)set->work, NULL, 2, 1,
      (double*)set->work, NULL, 2, 1, kinds, 0U);
  }

  /* Initialize working arrays. */
  memset(set->result,0U,2*Nk*sizeof(double complex));

  /* The first step. */

  /* Set the first 2*data->k_start coefficients to zero. */
  memset(set->work,0U,2*data->k_start*sizeof(double complex));

  work_ptr = &set->work[2*data->k_start];
  x_ptr = x;

  for (k = 0; k < k_end_tilde-data->k_start+1; k++)
  {
    *work_ptr++ = *x_ptr++;
    *work_ptr++ = 0;
  }

  /* Set the last 2*(set->N-1-k_end_tilde) coefficients to zero. */
  memset(&set->work[2*(k_end_tilde+1)],0U,2*(Nk-1-k_end_tilde)*sizeof(double complex));

  /* If k_end == Nk, use three-term recurrence to map last coefficient x_{Nk} to
   * x_{Nk-1} and x_{Nk-2}. */
  if (k_end == Nk)
  {
    set->work[2*(Nk-2)]   += data->gammaN[tk-2]*x[Nk-data->k_start];
    set->work[2*(Nk-1)]   += data->betaN[tk-2]*x[Nk-data->k_start];
    set->work[2*(Nk-1)+1]  = data->alphaN[tk-2]*x[Nk-data->k_start];
  }

  /*--------*/
  /*for (k = 0; k < 2*Nk; k++)
  {
    fprintf(stderr,"work[%2d] = %le + I*%le\tresult[%2d] = %le + I*%le\n",
      k,creal(set->work[k]),cimag(set->work[k]),k,creal(set->result[k]),
      cimag(set->result[k]));
  }*/
  /*--------*/

  /* Compute the remaining steps. */
  plength = 4;
  for (tau = 1; tau < tk; tau++)
  {
    /* Compute first l. */
    firstl = FIRST_L(k_start_tilde,plength);
    /* Compute last l. */
    lastl = LAST_L(k_end_tilde,plength);

    /* Compute the multiplication steps. */
    for (l = firstl; l <= lastl; l++)
    {
      /* Copy vectors to multiply into working arrays zero-padded to twice the length. */
      memcpy(set->vec3,&(set->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(double complex));
      memcpy(set->vec4,&(set->work[(plength/2)*(4*l+3)]),(plength/2)*sizeof(double complex));
      memset(&set->vec3[plength/2],0U,(plength/2)*sizeof(double complex));
      memset(&set->vec4[plength/2],0U,(plength/2)*sizeof(double complex));

      /* Copy coefficients into first half. */
      memcpy(&(set->work[(plength/2)*(4*l+2)]),&(set->work[(plength/2)*(4*l+1)]),(plength/2)*sizeof(double complex));
      memset(&(set->work[(plength/2)*(4*l+1)]),0U,(plength/2)*sizeof(double complex));
      memset(&(set->work[(plength/2)*(4*l+3)]),0U,(plength/2)*sizeof(double complex));

      /* Get matrix U_{n,tau,l} */
      step = &(data->steps[tau][l]);

      /* Check if step is stable. */
      if (step->stable)
      {
        /* Check, if we should do a symmetrizised step. */
        if (set->flags & FPT_AL_SYMMETRY && IS_SYMMETRIC(l,m,plength))
        {
          /*for (k = 0; k < plength; k++)
          {
            fprintf(stderr,"fpt_trafo: a11 = %le, a12 = %le, a21 = %le, a22 = %le\n",
              step->a11[0][k],step->a12[0][k],step->a21[0][k],step->a22[0][k]);
          }*/
          /* Multiply third and fourth polynomial with matrix U. */
          //fprintf(stderr,"\nhallo\n");
          fpt_do_step_symmetric(set->vec3, set->vec4, step->a11[0],
            step->a12[0], step->a21[0], step->a22[0], step->gamma[0], tau, set);
        }
        else
        {
          /* Multiply third and fourth polynomial with matrix U. */
          fpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], tau, set);
        }

        if (step->gamma[0] != 0.0)
        {
          for (k = 0; k < plength; k++)
          {
            set->work[plength*2*l+k] += set->vec3[k];
          }
        }
        for (k = 0; k < plength; k++)
        {
          set->work[plength*(2*l+1)+k] += set->vec4[k];
        }
      }
      else
      {
        /* Stabilize. */

        /* The lengh of the polynomials */
        plength_stab = step->N_stab;
        t_stab = step->t_stab;

        /*---------*/
        /*fprintf(stderr,"\nfpt_trafo: stabilizing at tau = %d, l = %d.\n",tau,l);
        fprintf(stderr,"\nfpt_trafo: plength_stab = %d.\n",plength_stab);
        fprintf(stderr,"\nfpt_trafo: tk = %d.\n",tk);
        fprintf(stderr,"\nfpt_trafo: index = %d.\n",tk-tau-1);*/
        /*---------*/

        /* Set rest of vectors explicitely to zero */
        /*fprintf(stderr,"fpt_trafo: stabilizing: plength = %d, plength_stab = %d\n",
          plength, plength_stab);*/
        memset(&set->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(double complex));
        memset(&set->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(double complex));

        /* Multiply third and fourth polynomial with matrix U. */
        if (m <= 1)
        {
          fpt_do_step_symmetric(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);
        }
        else if (m%2 == 0)
        {
          fpt_do_step_symmetric_u(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);
        }
        else
        {
          fpt_do_step_symmetric_l(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);
        }

        if (step->gamma[0] != 0.0)
        {
          for (k = 0; k < plength_stab; k++)
          {
            set->result[k] += set->vec3[k];
          }
        }

        for (k = 0; k < plength_stab; k++)
        {
          set->result[Nk+k] += set->vec4[k];
        }
      }
    }
    /* Double length of polynomials. */
    plength = plength<<1;

    /*--------*/
    /*for (k = 0; k < 2*Nk; k++)
    {
      fprintf(stderr,"work[%2d] = %le + I*%le\tresult[%2d] = %le + I*%le\n",
        k,creal(set->work[k]),cimag(set->work[k]),k,creal(set->result[k]),
        cimag(set->result[k]));
    }*/
    /*--------*/
  }

  /* Add the resulting cascade coeffcients to the coeffcients accumulated from
   * the stabilization steps. */
  for (k = 0; k < 2*Nk; k++)
  {
    set->result[k] += set->work[k];
  }

  /* The last step. Compute the Chebyshev coeffcients c_k^n from the
   * polynomials in front of P_0^n and P_1^n. */
  y[0] = data->gamma_m1*(set->result[0] + data->beta_0*set->result[Nk] +
    data->alpha_0*set->result[Nk+1]*0.5);
  y[1] = data->gamma_m1*(set->result[1] + data->beta_0*set->result[Nk+1]+
    data->alpha_0*(set->result[Nk]+set->result[Nk+2]*0.5));
  y[k_end-1] = data->gamma_m1*(set->result[k_end-1] +
    data->beta_0*set->result[Nk+k_end-1] +
    data->alpha_0*set->result[Nk+k_end-2]*0.5);
  y[k_end] = data->gamma_m1*(0.5*data->alpha_0*set->result[Nk+k_end-1]);
  for (k = 2; k <= k_end-2; k++)
  {
    y[k] = data->gamma_m1*(set->result[k] + data->beta_0*set->result[Nk+k] +
      data->alpha_0*0.5*(set->result[Nk+k-1]+set->result[Nk+k+1]));
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
    y[0] *= 2.0;
    fftw_execute_r2r(plan,(double*)y,(double*)y);
    fftw_destroy_plan(plan);
    for (k = 0; k <= k_end; k++)
    {
      y[k] *= 0.5;
    }
  }
}

void dpt_transposed(fpt_set set, const int m, double complex *x, 
  double complex *y, const int k_end, const unsigned int flags)
{
  int j;
  fpt_data *data = &(set->dpt[m]);
  int Nk;
  int tk;
  double norm;

  nfft_next_power_of_2_exp(k_end+1,&Nk,&tk);
  norm = 2.0/(Nk<<1);

  if (set->flags & FPT_NO_DIRECT_ALGORITHM)
  {
    return;
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
    for (j = 0; j <= k_end; j++)
    {
      set->xc_slow[j] = cos((PI*(j+0.5))/(k_end+1));
    }

    eval_sum_clenshaw_transposed(k_end, k_end, set->result, set->xc_slow,
      y, set->work, &data->alpha[1], &data->beta[1], &data->gamma[1],
      data->gamma_m1);

    memcpy(x,&set->result[data->k_start],(k_end-data->k_start+1)*
      sizeof(double complex));
  }
  else
  {
    memcpy(set->result,y,(k_end+1)*sizeof(double complex));
    memset(&set->result[k_end+1],0U,(Nk-k_end-1)*sizeof(double complex));

    for (j = 0; j < Nk; j++)
    {
      set->result[j] *= norm;
    }

    fftw_execute_r2r(set->plans_dct3[tk-2],(double*)set->result,
      (double*)set->result);

    eval_sum_clenshaw_transposed(k_end, Nk-1, set->temp, set->xcvecs[tk-2],
      set->result, set->work, &data->alpha[1], &data->beta[1], &data->gamma[1],
      data->gamma_m1);

    memcpy(x,&set->temp[data->k_start],(k_end-data->k_start+1)*sizeof(double complex));
  }
}

void fpt_transposed(fpt_set set, const int m, double complex *x, 
  const double complex *y, const int k_end, const unsigned int flags)
{
  /* Get transformation data. */
  fpt_data *data = &(set->dpt[m]);
  /** */
  int Nk;
  /** */
  int tk;
  /** */
  int k_start_tilde;
  /** */
  int k_end_tilde;

  /** Level index \f$tau\f$ */
  int tau;
  /** Index of first block at current level */
  int firstl;
  /** Index of last block at current level */
  int lastl;
  /** Block index \f$l\f$ */
  int l;
  /** Length of polynomial coefficient arrays at next level */
  int plength;
  /** Polynomial array length for stabilization */
  int plength_stab;
  /** Current matrix \f$U_{n,tau,l}\f$ */
  fpt_step *step;
  /** */
  fftw_plan plan;
  int length = k_end+1;
  fftw_r2r_kind kinds[2] = {FFTW_REDFT10,FFTW_REDFT10};
  /** Loop counter */
  int k;
  int t_stab;

  /* Check, if slow transformation should be used due to small bandwidth. */
  if (k_end < FPT_BREAK_EVEN)
  {
    /* Use NDSFT. */
    dpt_transposed(set, m, x, y, k_end, flags);
  }

  nfft_next_power_of_2_exp(k_end,&Nk,&tk);
  k_start_tilde = K_START_TILDE(data->k_start,Nk);
  k_end_tilde = K_END_TILDE(k_end,Nk);

  /* Check if fast transform is activated. */
  if (set->flags & FPT_NO_FAST_ALGORITHM)
  {
    return;
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
    plan = fftw_plan_many_r2r(1, &length, 2, (double*)set->work, NULL, 2, 1,
      (double*)set->work, NULL, 2, 1, kinds, 0U);
    fftw_execute_r2r(plan,(double*)y,(double*)set->result);
    fftw_destroy_plan(plan);
    for (k = 0; k <= k_end; k++)
    {
      set->result[k] *= 0.5;
    }
  }
  else
  {
    memcpy(set->result,y,(k_end+1)*sizeof(double complex));
  }

  /* Initialize working arrays. */
  memset(set->work,0U,2*Nk*sizeof(double complex));

  /* The last step is now the first step. */
  for (k = 0; k <= k_end; k++)
  {
    set->work[k] = data->gamma_m1*set->result[k];
  }
  //memset(&set->work[k_end+1],0U,(Nk+1-k_end)*sizeof(double complex));

  set->work[Nk] = data->gamma_m1*(data->beta_0*set->result[0] +
    data->alpha_0*set->result[1]);
  for (k = 1; k < k_end; k++)
  {
    set->work[Nk+k] = data->gamma_m1*(data->beta_0*set->result[k] +
      data->alpha_0*0.5*(set->result[k-1]+set->result[k+1]));
  }
  if (k_end<Nk)
  {
    memset(&set->work[k_end],0U,(Nk-k_end)*sizeof(double complex));
  }

  /** Save copy of inpute data for stabilization steps. */
  memcpy(set->result,set->work,2*Nk*sizeof(double complex));

  /* Compute the remaining steps. */
  plength = Nk;
  for (tau = tk-1; tau >= 1; tau--)
  {
    /* Compute first l. */
    firstl = FIRST_L(k_start_tilde,plength);
    /* Compute last l. */
    lastl = LAST_L(k_end_tilde,plength);

    /* Compute the multiplication steps. */
    for (l = firstl; l <= lastl; l++)
    {
      /* Initialize second half of coefficient arrays with zeros. */
      memcpy(set->vec3,&(set->work[(plength/2)*(4*l+0)]),plength*sizeof(double complex));
      memcpy(set->vec4,&(set->work[(plength/2)*(4*l+2)]),plength*sizeof(double complex));

      memcpy(&set->work[(plength/2)*(4*l+1)],&(set->work[(plength/2)*(4*l+2)]),
        (plength/2)*sizeof(double complex));

      /* Get matrix U_{n,tau,l} */
      step = &(data->steps[tau][l]);

      /* Check if step is stable. */
      if (step->stable)
      {
        if (set->flags & FPT_AL_SYMMETRY && IS_SYMMETRIC(l,m,plength))
        {
          /* Multiply third and fourth polynomial with matrix U. */
          fpt_do_step_transposed_symmetric(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], tau, set);
        }
        else
        {
          /* Multiply third and fourth polynomial with matrix U. */
          fpt_do_step_transposed(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], tau, set);
        }
        memcpy(&(set->vec3[plength/2]), set->vec4,(plength/2)*sizeof(double complex));

        for (k = 0; k < plength; k++)
        {
          set->work[plength*(4*l+2)/2+k] = set->vec3[k];
        }
      }
      else
      {
        /* Stabilize. */
        plength_stab = step->N_stab;
        t_stab = step->t_stab;

        memcpy(set->vec3,set->result,plength_stab*sizeof(double complex));
        memcpy(set->vec4,&(set->result[Nk]),plength_stab*sizeof(double complex));

        /* Multiply third and fourth polynomial with matrix U. */
        if (m <= 1)
        {
          fpt_do_step_transposed_symmetric(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);
        }
        else if (m%2 == 0)
        {
          fpt_do_step_transposed_symmetric_u(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);
        }
        else
        {
          fpt_do_step_transposed_symmetric_l(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);
        }

        memcpy(&(set->vec3[plength/2]),set->vec4,(plength/2)*sizeof(double complex));

        for (k = 0; k < plength; k++)
        {
          set->work[(plength/2)*(4*l+2)+k] = set->vec3[k];
        }
       }
    }
    /* Half the length of polynomial arrays. */
    plength = plength>>1;
  }

  /* First step */
  for (k = 0; k <= k_end_tilde-data->k_start; k++)
  {
    x[k] = set->work[2*(data->k_start+k)];
  }
  if (k_end == Nk)
  {
    x[Nk-data->k_start] =
        data->gammaN[tk-2]*set->work[2*(Nk-2)]
      + data->betaN[tk-2] *set->work[2*(Nk-1)]
      + data->alphaN[tk-2]*set->work[2*(Nk-1)+1];
  }
}

void fpt_finalize(fpt_set set)
{
  int tau;
  int l;
  int m;
  fpt_data *data;
  int k_start_tilde;
  int N_tilde;
  int firstl, lastl;
  int plength;

  /* TODO Clean up DPT transform data structures. */
  for (m = 0; m <= set->M; m++)
  {
    /* Check if precomputed. */
    data = &set->dpt[m];
    if (data->steps != (fpt_step**)NULL)
    {
      free(data->alphaN);
      free(data->betaN);
      free(data->gammaN);
      data->alphaN = NULL;
      data->betaN = NULL;
      data->gammaN = NULL;

      /* Free precomputed data. */
      k_start_tilde = K_START_TILDE(data->k_start,nfft_next_power_of_2(data->k_start)
        /*set->N*/);
      N_tilde = N_TILDE(set->N);
      plength = 4;
      for (tau = 1; tau < set->t; tau++)
      {
        /* Compute first l. */
        firstl = FIRST_L(k_start_tilde,plength);
        /* Compute last l. */
        lastl = LAST_L(N_tilde,plength);

        /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
        for (l = firstl; l <= lastl; l++)
        {
          /* Free components. */
          free(data->steps[tau][l].a11[0]);
          free(data->steps[tau][l].a12[0]);
          free(data->steps[tau][l].a21[0]);
          free(data->steps[tau][l].a22[0]);
          data->steps[tau][l].a11[0] = NULL;
          data->steps[tau][l].a12[0] = NULL;
          data->steps[tau][l].a21[0] = NULL;
          data->steps[tau][l].a22[0] = NULL;
          /* Free components. */
          free(data->steps[tau][l].a11);
          free(data->steps[tau][l].a12);
          free(data->steps[tau][l].a21);
          free(data->steps[tau][l].a22);
          free(data->steps[tau][l].gamma);
          data->steps[tau][l].a11 = NULL;
          data->steps[tau][l].a12 = NULL;
          data->steps[tau][l].a21 = NULL;
          data->steps[tau][l].a22 = NULL;
          data->steps[tau][l].gamma = NULL;
        }
        /* Free pointers for current level. */
        free(data->steps[tau]);
        data->steps[tau] = NULL;
        /* Double length of polynomials. */
        plength = plength<<1;
      }
      /* Free steps. */
      free(data->steps);
      data->steps = NULL;
    }

    if (set->flags & FPT_NO_DIRECT_ALGORITHM)
    {
    }
    else
    {
      /* Check, if recurrence coefficients must be copied. */
      //fprintf(stderr,"\nfpt_finalize: %d\n",set->flags & FPT_PERSISTENT_DATA);
      if (set->flags & FPT_PERSISTENT_DATA)
      {
      }
      else
      {
        free(data->alpha);
        free(data->beta);
        free(data->gamma);
      }
      data->alpha = NULL;
      data->beta = NULL;
      data->gamma = NULL;
    }
  }

  /* Delete array of DPT transform data. */
  free(set->dpt);
  set->dpt = NULL;

  for (tau = 1; tau < /*set->t*/set->t+1; tau++)
  {
    free(set->xcvecs[tau-1]);
    set->xcvecs[tau-1] = NULL;
  }
  free(set->xcvecs);
  set->xcvecs = NULL;

  /* Free auxilliary arrays. */
  free(set->work);
  free(set->result);

  /* Check if fast transform is activated. */
  if (set->flags & FPT_NO_FAST_ALGORITHM)
  {
  }
  else
  {
    /* Free auxilliary arrays. */
    free(set->vec3);
    free(set->vec4);
    free(set->z);
    set->work = NULL;
    set->result = NULL;
    set->vec3 = NULL;
    set->vec4 = NULL;
    set->z = NULL;

    /* Free FFTW plans. */
    for(tau = 0; tau < set->t/*-1*/; tau++)
    {
      fftw_destroy_plan(set->plans_dct3[tau]);
      fftw_destroy_plan(set->plans_dct2[tau]);
      set->plans_dct3[tau] = NULL;
      set->plans_dct2[tau] = NULL;
    }

    free(set->plans_dct3);
    free(set->plans_dct2);
    set->plans_dct3 = NULL;
    set->plans_dct2 = NULL;
  }

  //fprintf(stderr,"fpt_finalize: flags = %d\n",set->flags);

  if (set->flags & FPT_NO_DIRECT_ALGORITHM)
  {
  }
  else
  {
    /* Delete arrays of Chebyshev nodes. */
    free(set->xc_slow);
    set->xc_slow = NULL;
    free(set->temp);
    set->temp = NULL;
  }

  /* Free DPT set structure. */
  free(set);
}
