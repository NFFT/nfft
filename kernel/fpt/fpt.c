/*
 * Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
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

/**
 * \file fpt.c
 * \brief Implementation file for the FPT module
 * \author Jens Keiner
 */

/* configure header */
#include "config.h"

/* complex datatype (maybe) */
#ifdef HAVE_COMPLEX_H
#include<complex.h>
#endif

/* NFFT headers */
#include "nfft3.h"
#include "infft.h"

#include <stdbool.h>

#undef X
#define X(name) FPT(name)

/* Macros for index calculation. */

/** Minimum degree at top of a cascade */
#define K_START_TILDE(x,y) (MAX(MIN(x,y-2),0))

/** Maximum degree at top of a cascade */
#define K_END_TILDE(x,y) MIN(x,y-1)

/** Index of first block of four functions at level */
#define FIRST_L(x,y) (LRINT(FLOOR((x)/(R)y)))

/** Index of last block of four functions at level */
#define LAST_L(x,y) (LRINT(CEIL(((x)+1)/(R)y))-1)

#define N_TILDE(y) (y-1)

#define IS_SYMMETRIC(x,y,z) (x >= ((y-1.0)/z))

#define FPT_BREAK_EVEN 4

/**
 * Holds data for a single multiplication step in the cascade summation.
 */
typedef struct fpt_step_
{
  /** Whether this instance stores data for a stable fast step */
  bool stable;

  int Ns;                                 /**< TODO Add comment here.         */
  int ts;                                 /**< TODO Add comment here.         */
  R **a11,**a12,**a21,**a22;         /**< The matrix components          */
  R *g;                              /**<                                */
} fpt_step;

/**
 * Holds data for a single cascade summation.
 */
typedef struct fpt_data_
{
  fpt_step **steps;                       /**< The cascade summation steps    */
  int k_start;                            /**< TODO Add comment here.         */
  R *alphaN;                         /**< TODO Add comment here.         */
  R *betaN;                          /**< TODO Add comment here.         */
  R *gammaN;                         /**< TODO Add comment here.         */
  R alpha_0;                         /**< TODO Add comment here.         */
  R beta_0;                          /**< TODO Add comment here.         */
  R gamma_m1;                        /**< TODO Add comment here.         */
  /* Data for direct transform. */        /**< TODO Add comment here.         */
  R *_alpha;                          /**< TODO Add comment here.         */
  R *_beta;                           /**< TODO Add comment here.         */
  R *_gamma;                          /**< TODO Add comment here.         */
} fpt_data;

/**
 * Holds data for a set of cascade summations.
 */
typedef struct fpt_set_s_
{
  unsigned int flags;                              /**< The flags                     */
  int M;                                  /**< The number of DPT transforms  */
  int N;                                  /**< The transform length. Must be
                                               a power of two.               */
  int t;                                  /**< The exponent of N             */
  fpt_data *dpt;                          /**< The DPT transform data        */
  R **xcvecs;                        /**< Array of pointers to arrays
                                               containing the Chebyshev
                                               nodes                         */
  R *xc;                             /**< Array for Chebychev-nodes.    */
  C *temp;                          /**< */
  C *work;                          /**< */
  C *result;                        /**< */
  C *vec3;
  C *vec4;
  C *z;
  FFTW(plan) *plans_dct3;                  /**< Transform plans for the fftw
                                               library                       */
  FFTW(plan) *plans_dct2;                  /**< Transform plans for the fftws
                                               library                       */
  FFTW(r2r_kind) *kinds;                   /**< Transform kinds for fftw
                                               library                       */
  FFTW(r2r_kind) *kindsr;                  /**< Transform kinds for fftw
                                               library                       */

  int *lengths; /**< Transform lengths for fftw library */

  /* Data for slow transforms. */
  R *xc_slow;
} fpt_set_s;

static inline void abuvxpwy(R a, R b, C* u, C* x,
  R* v, C* y, R* w, int n)
{
  int l; C *u_ptr = u, *x_ptr = x, *y_ptr = y;
  R *v_ptr = v, *w_ptr = w;
  for (l = 0; l < n; l++)
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++));
}

#define ABUVXPWY_SYMMETRIC(NAME,S1,S2) \
static inline void NAME(R a, R b, C* u, C* x, \
  R* v, C* y, R* w, int n) \
{ \
  const int n2 = n>>1; \
  int l; C *u_ptr = u, *x_ptr = x, *y_ptr = y; \
  R *v_ptr = v, *w_ptr = w; \
  for (l = 0; l < n2; l++) \
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++)); \
  v_ptr--; w_ptr--; \
  for (l = 0; l < n2; l++) \
    *u_ptr++ = a * (b * S1 * (*v_ptr--) * (*x_ptr++) + S2 * (*w_ptr--) * (*y_ptr++)); \
}

ABUVXPWY_SYMMETRIC(abuvxpwy_symmetric1,1.0,-1.0)
ABUVXPWY_SYMMETRIC(abuvxpwy_symmetric2,-1.0,1.0)

#define ABUVXPWY_SYMMETRIC_1(NAME,S1) \
static inline void NAME(R a, R b, C* u, C* x, \
  R* v, C* y, int n, R *xx) \
{ \
  const int n2 = n>>1; \
  int l; C *u_ptr = u, *x_ptr = x, *y_ptr = y; \
  R *v_ptr = v, *xx_ptr = xx; \
  for (l = 0; l < n2; l++, v_ptr++) \
    *u_ptr++ = a * (b * (*v_ptr) * (*x_ptr++) + ((*v_ptr)*(1.0+*xx_ptr++)) * (*y_ptr++)); \
  v_ptr--; \
  for (l = 0; l < n2; l++, v_ptr--) \
    *u_ptr++ = a * (b * S1 * (*v_ptr) * (*x_ptr++) + (S1 * (*v_ptr) * (1.0+*xx_ptr++)) * (*y_ptr++)); \
}

ABUVXPWY_SYMMETRIC_1(abuvxpwy_symmetric1_1,1.0)
ABUVXPWY_SYMMETRIC_1(abuvxpwy_symmetric1_2,-1.0)

#define ABUVXPWY_SYMMETRIC_2(NAME,S1) \
static inline void NAME(R a, R b, C* u, C* x, \
  C* y, R* w, int n, R *xx) \
{ \
  const int n2 = n>>1; \
  int l; C *u_ptr = u, *x_ptr = x, *y_ptr = y; \
  R *w_ptr = w, *xx_ptr = xx; \
  for (l = 0; l < n2; l++, w_ptr++) \
    *u_ptr++ = a * (b * (*w_ptr/(1.0+*xx_ptr++)) * (*x_ptr++) + (*w_ptr) * (*y_ptr++)); \
  w_ptr--; \
  for (l = 0; l < n2; l++, w_ptr--) \
    *u_ptr++ = a * (b * (S1 * (*w_ptr)/(1.0+*xx_ptr++) ) * (*x_ptr++) + S1 * (*w_ptr) * (*y_ptr++)); \
}

ABUVXPWY_SYMMETRIC_2(abuvxpwy_symmetric2_1,1.0)
ABUVXPWY_SYMMETRIC_2(abuvxpwy_symmetric2_2,-1.0)

static inline void auvxpwy(R a, C* u, C* x, R* v,
  C* y, R* w, int n)
{
  int l;
  C *u_ptr = u, *x_ptr = x, *y_ptr = y;
  R *v_ptr = v, *w_ptr = w;
  for (l = n; l > 0; l--)
    *u_ptr++ = a * ((*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++));
}

static inline void auvxpwy_symmetric(R a, C* u, C* x,
  R* v, C* y, R* w, int n)
{
  const int n2 = n>>1; \
  int l;
  C *u_ptr = u, *x_ptr = x, *y_ptr = y;
  R *v_ptr = v, *w_ptr = w;
  for (l = n2; l > 0; l--)
    *u_ptr++ = a * ((*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++));
  v_ptr--; w_ptr--;
  for (l = n2; l > 0; l--)
    *u_ptr++ = a * ((*v_ptr--) * (*x_ptr++) - (*w_ptr--) * (*y_ptr++));
}

static inline void auvxpwy_symmetric_1(R a, C* u, C* x,
  R* v, C* y, R* w, int n, R *xx)
{
  const int n2 = n>>1; \
  int l;
  C *u_ptr = u, *x_ptr = x, *y_ptr = y;
  R *v_ptr = v, *w_ptr = w, *xx_ptr = xx;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (((*v_ptr++)*(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr++)*(1.0+*xx_ptr)) * (*y_ptr++));
  v_ptr--; w_ptr--;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (-((*v_ptr--)*(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr--)*(1.0+*xx_ptr)) * (*y_ptr++));
}

static inline void auvxpwy_symmetric_2(R a, C* u, C* x,
  R* v, C* y, R* w, int n, R *xx)
{
  const int n2 = n>>1; \
  int l;
  C *u_ptr = u, *x_ptr = x, *y_ptr = y;
  R *v_ptr = v, *w_ptr = w, *xx_ptr = xx;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (((*v_ptr++)/(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr++)/(1.0+*xx_ptr)) * (*y_ptr++));
  v_ptr--; w_ptr--;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (-((*v_ptr--)/(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr--)/(1.0+*xx_ptr)) * (*y_ptr++));
}

#define FPT_DO_STEP(NAME,M1_FUNCTION,M2_FUNCTION) \
static inline void NAME(C  *a, C *b, R *a11, R *a12, \
  R *a21, R *a22, R g, int tau, fpt_set set) \
{ \
  /** The length of the coefficient arrays. */ \
  int length = 1<<(tau+1); \
  /** Twice the length of the coefficient arrays. */ \
  R norm = 1.0/(length<<1); \
  \
  /* Compensate for factors introduced by a raw DCT-III. */ \
  a[0] *= 2.0; \
  b[0] *= 2.0; \
  \
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */ \
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)a,(R*)a); \
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)b,(R*)b); \
  \
  /*for (k = 0; k < length; k++)*/ \
  /*{*/ \
    /*fprintf(stderr,"fpt_do_step: a11 = %le, a12 = %le, a21 = %le, a22 = %le\n",*/ \
    /*  a11[k],a12[k],a21[k],a22[k]);*/ \
  /*}*/ \
  \
  /* Check, if gamma is zero. */ \
  if (g == 0.0) \
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
    M1_FUNCTION(norm*g,a,a,a11,b,a12,length); \
    memcpy(b,set->z, (size_t)(length) * sizeof(C)); \
    /* Compute Chebyshev-coefficients using a DCT-II. */ \
    fftw_execute_r2r(set->plans_dct2[tau-1],(R*)a,(R*)a); \
    /* Compensate for factors introduced by a raw DCT-II. */ \
    a[0] *= 0.5; \
  } \
  \
  /* Compute Chebyshev-coefficients using a DCT-II. */ \
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)b,(R*)b); \
  /* Compensate for factors introduced by a raw DCT-II. */ \
  b[0] *= 0.5; \
}

FPT_DO_STEP(fpt_do_step,auvxpwy,auvxpwy)
FPT_DO_STEP(fpt_do_step_symmetric,auvxpwy_symmetric,auvxpwy_symmetric)
/*FPT_DO_STEP(fpt_do_step_symmetric_u,auvxpwy_symmetric,auvxpwy)
FPT_DO_STEP(fpt_do_step_symmetric_l,auvxpwy,auvxpwy_symmetric)*/

static inline void fpt_do_step_symmetric_u(C *a, C *b,
  R *a11, R *a12, R *a21, R *a22, R *x,
  R gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  R norm = 1.0/(length<<1);

  UNUSED(a21); UNUSED(a22);

  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)a,(R*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)b,(R*)b);

  /*for (k = 0; k < length; k++)*/
  /*{*/
    /*fprintf(stderr,"fpt_do_step: a11 = %le, a12 = %le, a21 = %le, a22 = %le\n",*/
    /*  a11[k],a12[k],a21[k],a22[k]);*/
  /*}*/

  /* Check, if gamma is zero. */
  if (gam == 0.0)
  {
    /*fprintf(stderr,"gamma == 0!\n");*/
    /* Perform multiplication only for second row. */
    auvxpwy_symmetric_1(norm,b,b,a12,a,a11,length,x);
  }
  else
  {
    /*fprintf(stderr,"gamma != 0!\n");*/
    /* Perform multiplication for both rows. */
    auvxpwy_symmetric_1(norm,set->z,b,a12,a,a11,length,x);
    auvxpwy_symmetric(norm*gam,a,a,a11,b,a12,length);
    memcpy(b,set->z,(size_t)(length) * sizeof(C));
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(set->plans_dct2[tau-1],(R*)a,(R*)a);
    /* Compensate for factors introduced by a raw DCT-II. */
    a[0] *= 0.5;
  }

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)b,(R*)b);
  /* Compensate for factors introduced by a raw DCT-II. */
  b[0] *= 0.5;
}

static inline void fpt_do_step_symmetric_l(C  *a, C *b,
  R *a11, R *a12, R *a21, R *a22, R *x, R gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  R norm = 1.0/(length<<1);

  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;

  UNUSED(a11); UNUSED(a12);

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)a,(R*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)b,(R*)b);

  /*for (k = 0; k < length; k++)*/
  /*{*/
    /*fprintf(stderr,"fpt_do_step: a11 = %le, a12 = %le, a21 = %le, a22 = %le\n",*/
    /*  a11[k],a12[k],a21[k],a22[k]);*/
  /*}*/

  /* Check, if gamma is zero. */
  if (gam == 0.0)
  {
    /*fprintf(stderr,"gamma == 0!\n");*/
    /* Perform multiplication only for second row. */
    auvxpwy_symmetric(norm,b,b,a22,a,a21,length);
  }
  else
  {
    /*fprintf(stderr,"gamma != 0!\n");*/
    /* Perform multiplication for both rows. */
    auvxpwy_symmetric(norm,set->z,b,a22,a,a21,length);
    auvxpwy_symmetric_2(norm*gam,a,a,a21,b,a22,length,x);
    memcpy(b,set->z, (size_t)(length) * sizeof(C));
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(set->plans_dct2[tau-1],(R*)a,(R*)a);
    /* Compensate for factors introduced by a raw DCT-II. */
    a[0] *= 0.5;
  }

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)b,(R*)b);
  /* Compensate for factors introduced by a raw DCT-II. */
  b[0] *= 0.5;
}

#define FPT_DO_STEP_TRANSPOSED(NAME,M1_FUNCTION,M2_FUNCTION) \
static inline void NAME(C  *a, C *b, R *a11, \
  R *a12, R *a21, R *a22, R g, int tau, fpt_set set) \
{ \
  /** The length of the coefficient arrays. */ \
  int length = 1<<(tau+1); \
  /** Twice the length of the coefficient arrays. */ \
  R norm = 1.0/(length<<1); \
  \
  /* Compute function values from Chebyshev-coefficients using a DCT-III. */ \
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)a,(R*)a); \
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)b,(R*)b); \
  \
  /* Perform matrix multiplication. */ \
  M1_FUNCTION(norm,g,set->z,a,a11,b,a21,length); \
  M2_FUNCTION(norm,g,b,a,a12,b,a22,length); \
  memcpy(a,set->z, (size_t)(length) * sizeof(C)); \
  \
  /* Compute Chebyshev-coefficients using a DCT-II. */ \
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)a,(R*)a); \
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)b,(R*)b); \
}

FPT_DO_STEP_TRANSPOSED(fpt_do_step_t,abuvxpwy,abuvxpwy)
FPT_DO_STEP_TRANSPOSED(fpt_do_step_t_symmetric,abuvxpwy_symmetric1,abuvxpwy_symmetric2)
/*FPT_DO_STEP_TRANSPOSED(fpt_do_step_t_symmetric_u,abuvxpwy_symmetric1_1,abuvxpwy_symmetric1_2)*/
/*FPT_DO_STEP_TRANSPOSED(fpt_do_step_t_symmetric_l,abuvxpwy_symmetric2_2,abuvxpwy_symmetric2_1)*/


static inline void fpt_do_step_t_symmetric_u(C  *a,
  C *b, R *a11, R *a12, R *x,
  R gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  R norm = 1.0/(length<<1);

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)a,(R*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)b,(R*)b);

  /* Perform matrix multiplication. */
  abuvxpwy_symmetric1_1(norm,gam,set->z,a,a11,b,length,x);
  abuvxpwy_symmetric1_2(norm,gam,b,a,a12,b,length,x);
  memcpy(a,set->z, (size_t)(length) * sizeof(C));

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)a,(R*)a);
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)b,(R*)b);
}

static inline void fpt_do_step_t_symmetric_l(C  *a,
  C *b, R *a21, R *a22,
  R *x, R gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  R norm = 1.0/(length<<1);

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)a,(R*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(R*)b,(R*)b);

  /* Perform matrix multiplication. */
  abuvxpwy_symmetric2_2(norm,gam,set->z,a,b,a21,length,x);
  abuvxpwy_symmetric2_1(norm,gam,b,a,b,a22,length,x);
  memcpy(a,set->z, (size_t)(length) * sizeof(C));

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)a,(R*)a);
  fftw_execute_r2r(set->plans_dct2[tau-1],(R*)b,(R*)b);
}


static void eval_clenshaw(const R *x, R *y, int size, int k, const R *alpha,
  const R *beta, const R *gam)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  R a,b,x_val_act,a_old;
  const R *x_act;
  R *y_act;
  const R *alpha_act, *beta_act, *gamma_act;

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
      gamma_act = &(gam[k]);
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

static void eval_clenshaw2(const R *x, R *z, R *y, int size1, int size, int k, const R *alpha,
  const R *beta, const R *gam)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  R a,b,x_val_act,a_old;
  const R *x_act;
  R *y_act, *z_act;
  const R *alpha_act, *beta_act, *gamma_act;

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  z_act = z;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;

    if (k == 0)
    {
      *y_act = 1.0;
      *z_act = 0.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gam[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
        b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      if (i < size1)
        *z_act = a;
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
    }

    x_act++;
    y_act++;
    z_act++;
  }
  /*gamma_act++;
  FILE *f = fopen("/Users/keiner/Desktop/nfsft_debug.txt","a");
  fprintf(f,"size1: %10d, size: %10d\n",size1,size);
  fclose(f);*/
}

static int eval_clenshaw_thresh2(const R *x, R *z, R *y, int size, int k,
  const R *alpha, const R *beta, const R *gam, const
  R threshold)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  R a,b,x_val_act,a_old;
  const R *x_act;
  R *y_act, *z_act;
  const R *alpha_act, *beta_act, *gamma_act;
  R max = -nfft_float_property(NFFT_R_MAX);
  const R t = LOG10(FABS(threshold));

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  z_act = z;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;

    if (k == 0)
    {
     *y_act = 1.0;
     *z_act = 0.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gam[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
         b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *z_act = a;
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
      max = FMAX(max,LOG10(FABS(*y_act)));
      if (max > t)
        return 1;
    }
    x_act++;
    y_act++;
    z_act++;
  }
  return 0;
}

static inline void eval_sum_clenshaw_fast(const int N, const int M,
  const C *a, const R *x, C *y,
  const R *alpha, const R *beta, const R *gam,
  const R lambda)
{
  int j,k;
  C tmp1, tmp2, tmp3;
  R xc;
  
  /*fprintf(stderr, "Executing eval_sum_clenshaw_fast.\n");  
  fprintf(stderr, "Before transform:\n");  
  for (j = 0; j < N; j++)
    fprintf(stderr, "a[%4d] = %e.\n", j, a[j]);  
  for (j = 0; j <= M; j++)
    fprintf(stderr, "x[%4d] = %e, y[%4d] = %e.\n", j, x[j], j, y[j]);*/

  if (N == 0)
    for (j = 0; j <= M; j++)
      y[j] = a[0];
  else
  {
    for (j = 0; j <= M; j++)
    {
#if 0
      xc = x[j];
      tmp2 = a[N];
      tmp1 = a[N-1] + (alpha[N-1] * xc + beta[N-1])*tmp2;
      for (k = N-1; k > 0; k--)
      {
        tmp3 =   a[k-1]
               + (alpha[k-1] * xc + beta[k-1]) * tmp1
               + gam[k] * tmp2;
        tmp2 = tmp1;
        tmp1 = tmp3;
      }
      y[j] = lambda * tmp1;
#else
      xc = x[j];
      tmp1 = a[N-1];
      tmp2 = a[N];
      for (k = N-1; k > 0; k--)
      {
        tmp3 = a[k-1] + tmp2 * gam[k];
        tmp2 *= (alpha[k] * xc + beta[k]);
        tmp2 += tmp1;
        tmp1 = tmp3;
        /*if (j == 1515) 
        {
          fprintf(stderr, "k = %d, tmp1 = %e, tmp2 = %e.\n", k, tmp1, tmp2);  
        }*/
      }
      tmp2 *= (alpha[0] * xc + beta[0]);
        //fprintf(stderr, "alpha[0] = %e, beta[0] = %e.\n", alpha[0], beta[0]);  
      y[j] = lambda * (tmp2 + tmp1);
        //fprintf(stderr, "lambda = %e.\n", lambda);  
#endif
    }
  }
  /*fprintf(stderr, "Before transform:\n");  
  for (j = 0; j < N; j++)
    fprintf(stderr, "a[%4d] = %e.\n", j, a[j]);  
  for (j = 0; j <= M; j++)
    fprintf(stderr, "x[%4d] = %e, y[%4d] = %e.\n", j, x[j], j, y[j]);  */
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
 * for given complex coefficients \f$\left(a_k\right)_{k=0}^N \in
 * \mathbb{C}^{N+1}\f$ at given nodes \f$\left(x_j\right)_{j=0}^M \in
 * \mathbb{R}^{M+1}\f$, \f$M \in \mathbb{N}_0\f$.
 */
static void eval_sum_clenshaw_transposed(int N, int M, C* a, R *x,
  C *y, C *temp, R *alpha, R *beta, R *gam,
  R lambda)
{
  int j,k;
  C* it1 = temp;
  C* it2 = y;
  C aux;

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
        it2[j] = it2[j]*(alpha[k-1] * x[j] + beta[k-1]) + gam[k-1] * aux;
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
#ifdef _OPENMP
  int nthreads = X(get_num_threads)();
#endif

  /* Allocate memory for new DPT set. */
  fpt_set_s *set = (fpt_set_s*)Y(malloc)(sizeof(fpt_set_s));

  /* Save parameters in structure. */
  set->flags = flags;

  set->M = M;
  set->t = t;
  set->N = 1<<t;

  /* Allocate memory for M transforms. */
  set->dpt = (fpt_data*) Y(malloc)((size_t)(M) * sizeof(fpt_data));

  /* Initialize with NULL pointer. */
  for (m = 0; m < set->M; m++)
    set->dpt[m].steps = 0;

  /* Create arrays with Chebyshev nodes. */

  /* Initialize array with Chebyshev coefficients for the polynomial x. This
   * would be trivially an array containing a 1 as second entry with all other
   * coefficients set to zero. In order to compensate for the multiplicative
   * factor 2 introduced by the DCT-III, we set this coefficient to 0.5 here. */

  /* Allocate memory for array of pointers to node arrays. */
  set->xcvecs = (R**) Y(malloc)((size_t)(set->t) * sizeof(R*));
  /* For each polynomial length starting with 4, compute the Chebyshev nodes
   * using a DCT-III. */
  plength = 4;
  for (tau = 1; tau < t+1; tau++)
  {
    /* Allocate memory for current array. */
    set->xcvecs[tau-1] = (R*) Y(malloc)((size_t)(plength) * sizeof(R));
    for (k = 0; k < plength; k++)
    {
      set->xcvecs[tau-1][k] = cos(((k+0.5)*KPI)/plength);
    }
    plength = plength << 1;
  }

  /** Allocate memory for auxilliary arrays. */
  set->work = (C*) Y(malloc)((size_t)(2 * set->N) * sizeof(C));
  set->result = (C*) Y(malloc)((size_t)(2 * set->N) * sizeof(C));

  set->lengths = (int*) Y(malloc)((size_t)(set->t) * sizeof(int));
  set->plans_dct2 = (fftw_plan*) Y(malloc)((size_t)(set->t) * sizeof(fftw_plan));
  set->kindsr     = (fftw_r2r_kind*) Y(malloc)(2 * sizeof(fftw_r2r_kind));
  set->kindsr[0]  = FFTW_REDFT10;
  set->kindsr[1]  = FFTW_REDFT10;
  for (tau = 0, plength = 4; tau < set->t/*-1*/; tau++, plength<<=1)
  {
    set->lengths[tau] = plength;
#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    fftw_plan_with_nthreads(nthreads);
#endif
    set->plans_dct2[tau] =
      fftw_plan_many_r2r(1, &set->lengths[tau], 2, (R*)set->work, NULL,
                         2, 1, (R*)set->result, NULL, 2, 1,set->kindsr,
                         0);
#ifdef _OPENMP
}
#endif
  }

  /* Check if fast transform is activated. */
  if (!(set->flags & FPT_NO_FAST_ALGORITHM))
  {
    /** Allocate memory for auxilliary arrays. */
    set->vec3 = (C*) Y(malloc)((size_t)(set->N) * sizeof(C));
    set->vec4 = (C*) Y(malloc)((size_t)(set->N) * sizeof(C));
    set->z = (C*) Y(malloc)((size_t)(set->N) * sizeof(C));

    /** Initialize FFTW plans. */
    set->plans_dct3 = (fftw_plan*) Y(malloc)((size_t)(set->t) * sizeof(fftw_plan));
    set->kinds      = (fftw_r2r_kind*) Y(malloc)(2 * sizeof(fftw_r2r_kind));
    set->kinds[0]   = FFTW_REDFT01;
    set->kinds[1]   = FFTW_REDFT01;
    for (tau = 0, plength = 4; tau < set->t/*-1*/; tau++, plength<<=1)
    {
      set->lengths[tau] = plength;
#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    fftw_plan_with_nthreads(nthreads);
#endif
      set->plans_dct3[tau] =
        fftw_plan_many_r2r(1, &set->lengths[tau], 2, (R*)set->work, NULL,
                           2, 1, (R*)set->result, NULL, 2, 1, set->kinds,
                           0);
#ifdef _OPENMP
}
#endif
    }
    nfft_free(set->lengths);
    nfft_free(set->kinds);
    nfft_free(set->kindsr);
    set->lengths = NULL;
    set->kinds = NULL;
    set->kindsr = NULL;
  }

  if (!(set->flags & FPT_NO_DIRECT_ALGORITHM))
  {
    set->xc_slow = (R*) Y(malloc)((size_t)(set->N + 1) * sizeof(R));
    set->temp = (C*) Y(malloc)((size_t)(set->N + 1) * sizeof(C));
  }

  /* Return the newly created DPT set. */
  return set;
}

void fpt_precompute(fpt_set set, const int m, R *alpha, R *beta,
  R *gam, int k_start, const R threshold)
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
  R *a11;      /**< Array containing function values of the
                         (1,1)-component of U_k^n.                           */
  R *a12;      /**< Array containing function values of the
                         (1,2)-component of U_k^n.                           */
  R *a21;      /**< Array containing function values of the
                         (2,1)-component of U_k^n.                           */
  R *a22;      /**< Array containing function values of the
                         (2,2)-component of U_k^n.                           */
  const R *calpha;
  const R *cbeta;
  const R *cgamma;
  int needstab = 0; /**< Used to indicate that stabilization is neccessary.  */
  int k_start_tilde;
  int N_tilde;
  int clength;
  int clength_1;
  int clength_2;
  int t_stab, N_stab;
  fpt_data *data;

  /* Get pointer to DPT data. */
  data = &(set->dpt[m]);

  /* Check, if already precomputed. */
  if (data->steps != NULL)
    return;

  /* Save k_start. */
  data->k_start = k_start;

  data->gamma_m1 = gam[0];

  /* Check if fast transform is activated. */
  if (!(set->flags & FPT_NO_FAST_ALGORITHM))
  {
    /* Save recursion coefficients. */
    data->alphaN = (R*) Y(malloc)((size_t)(set->t - 1) * sizeof(C));
    data->betaN = (R*) Y(malloc)((size_t)(set->t - 1) * sizeof(C));
    data->gammaN = (R*) Y(malloc)((size_t)(set->t - 1) * sizeof(C));

    for (tau = 2; tau <= set->t; tau++)
    {

      data->alphaN[tau-2] = alpha[1<<tau];
      data->betaN[tau-2] = beta[1<<tau];
      data->gammaN[tau-2] = gam[1<<tau];
    }

    data->alpha_0 = alpha[1];
    data->beta_0 = beta[1];

    k_start_tilde = K_START_TILDE(data->k_start,Y(next_power_of_2)(data->k_start)
      /*set->N*/);
    N_tilde = N_TILDE(set->N);

    /* Allocate memory for the cascade with t = log_2(N) many levels. */
    data->steps = (fpt_step**) Y(malloc)((size_t)(set->t) * sizeof(fpt_step*));

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
      data->steps[tau] = (fpt_step*) Y(malloc)(sizeof(fpt_step)
                         * (size_t)(lastl + 1));

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
        a11 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);
        a12 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);
        a21 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);
        a22 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);

        /* Evaluate the associated polynomials at the 2^{tau+1} Chebyshev
         * nodes. */

        /* Get the pointers to the three-term recurrence coeffcients. */
        calpha = &(alpha[plength*l+1+1]);
        cbeta = &(beta[plength*l+1+1]);
        cgamma = &(gam[plength*l+1+1]);

        if (set->flags & FPT_NO_STABILIZATION)
        {
          /* Evaluate P_{2^{tau}-2}^n(\cdot,2^{tau+1}l+2). */
          calpha--;
          cbeta--;
          cgamma--;
          eval_clenshaw2(set->xcvecs[tau-1], a11, a21, clength, clength, degree-1, calpha, cbeta,
            cgamma);
          eval_clenshaw2(set->xcvecs[tau-1], a12, a22, clength, clength, degree, calpha, cbeta,
            cgamma);
          needstab = 0;
        }
        else
        {
          calpha--;
          cbeta--;
          cgamma--;
          /* Evaluate P_{2^{tau}-1}^n(\cdot,2^{tau+1}l+1). */
          needstab = eval_clenshaw_thresh2(set->xcvecs[tau-1], a11, a21, clength,
            degree-1, calpha, cbeta, cgamma, threshold);
          if (needstab == 0)
          {
            /* Evaluate P_{2^{tau}}^n(\cdot,2^{tau+1}l+1). */
            needstab = eval_clenshaw_thresh2(set->xcvecs[tau-1], a12, a22, clength,
              degree, calpha, cbeta, cgamma, threshold);
          }
        }
        

        /* Check if stabilization needed. */
        if (needstab == 0)
        {
          data->steps[tau][l].a11 = (R**) Y(malloc)(sizeof(R*));
          data->steps[tau][l].a12 = (R**) Y(malloc)(sizeof(R*));
          data->steps[tau][l].a21 = (R**) Y(malloc)(sizeof(R*));
          data->steps[tau][l].a22 = (R**) Y(malloc)(sizeof(R*));
          data->steps[tau][l].g = (R*) Y(malloc)(sizeof(R));
          /* No stabilization needed. */
          data->steps[tau][l].a11[0] = a11;
          data->steps[tau][l].a12[0] = a12;
          data->steps[tau][l].a21[0] = a21;
          data->steps[tau][l].a22[0] = a22;
          data->steps[tau][l].g[0] = gam[plength*l+1+1];
          data->steps[tau][l].stable = true;
        }
        else
        {
          /* Stabilize. */
          degree_stab = degree*(2*l+1);
          Y(next_power_of_2_exp)((l+1)*(1<<(tau+1)),&N_stab,&t_stab);

          /* Old arrays are to small. */
          nfft_free(a11);
          nfft_free(a12);
          nfft_free(a21);
          nfft_free(a22);

          data->steps[tau][l].a11 = (R**) Y(malloc)(sizeof(R*));
          data->steps[tau][l].a12 = (R**)Y(malloc)(sizeof(R*));
          data->steps[tau][l].a21 = (R**) Y(malloc)(sizeof(R*));
          data->steps[tau][l].a22 = (R**) Y(malloc)(sizeof(R*));
          data->steps[tau][l].g = (R*) Y(malloc)(sizeof(R));

          plength_stab = N_stab;

          if (set->flags & FPT_AL_SYMMETRY)
          {
            if (m <= 1)
            {
              /* This should never be executed */
              clength_1 = plength_stab;
              clength_2 = plength_stab;
              /* Allocate memory for arrays. */
              a11 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_1);
              a12 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_1);
              a21 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_2);
              a22 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_2);
              /* Get the pointers to the three-term recurrence coeffcients. */
              calpha = &(alpha[1]); cbeta = &(beta[1]); cgamma = &(gam[1]);
              eval_clenshaw2(set->xcvecs[t_stab-2], a11, a21, clength_1,
                clength_2, degree_stab-1, calpha, cbeta, cgamma);
              eval_clenshaw2(set->xcvecs[t_stab-2], a12, a22, clength_1,
                clength_2, degree_stab+0, calpha, cbeta, cgamma);
            }
            else
            {
              clength = plength_stab/2;
              if (m%2 == 0)
              {
                a11 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);
                a12 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);
                a21 = 0;
                a22 = 0;
                calpha = &(alpha[2]); cbeta = &(beta[2]); cgamma = &(gam[2]);
                eval_clenshaw(set->xcvecs[t_stab-2], a11, clength,
                  degree_stab-2, calpha, cbeta, cgamma);
                eval_clenshaw(set->xcvecs[t_stab-2], a12, clength,
                  degree_stab-1, calpha, cbeta, cgamma);
              }
              else
              {
                a11 = 0;
                a12 = 0;
                a21 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);
                a22 = (R*) Y(malloc)(sizeof(R) * (size_t)clength);
                calpha = &(alpha[1]); cbeta = &(beta[1]); cgamma = &(gam[1]);
                eval_clenshaw(set->xcvecs[t_stab-2], a21, clength,
                  degree_stab-1,calpha, cbeta, cgamma);
                eval_clenshaw(set->xcvecs[t_stab-2], a22, clength,
                  degree_stab+0, calpha, cbeta, cgamma);
              }
            }
          }
          else
          {
            clength_1 = plength_stab;
            clength_2 = plength_stab;
            a11 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_1);
            a12 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_1);
            a21 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_2);
            a22 = (R*) Y(malloc)(sizeof(R) * (size_t)clength_2);
            calpha = &(alpha[2]);
            cbeta = &(beta[2]);
            cgamma = &(gam[2]);
            calpha--;
            cbeta--;
            cgamma--;
            eval_clenshaw2(set->xcvecs[t_stab-2], a11, a21, clength_1, clength_2, degree_stab-1,
              calpha, cbeta, cgamma);
            eval_clenshaw2(set->xcvecs[t_stab-2], a12, a22, clength_1, clength_2, degree_stab+0,
              calpha, cbeta, cgamma);

          }
          data->steps[tau][l].a11[0] = a11;
          data->steps[tau][l].a12[0] = a12;
          data->steps[tau][l].a21[0] = a21;
          data->steps[tau][l].a22[0] = a22;

          data->steps[tau][l].g[0] =  gam[1+1];
          data->steps[tau][l].stable = false;
          data->steps[tau][l].ts = t_stab;
          data->steps[tau][l].Ns = N_stab;
        }
      }
      /** Increase polynomial degree to next power of two. */
      plength = plength << 1;
    }
  }

  if (!(set->flags & FPT_NO_DIRECT_ALGORITHM))
  {
    /* Check, if recurrence coefficients must be copied. */
    if (set->flags & FPT_PERSISTENT_DATA)
    {
      data->_alpha = (R*) alpha;
      data->_beta = (R*) beta;
      data->_gamma = (R*) gam;
    }
    else
    {
      data->_alpha = (R*) Y(malloc)((size_t)(set->N + 1) * sizeof(R));
      data->_beta = (R*) Y(malloc)((size_t)(set->N + 1) * sizeof(R));
      data->_gamma = (R*) Y(malloc)((size_t)(set->N + 1) * sizeof(R));
      memcpy(data->_alpha,alpha, (size_t)((set->N + 1)) * sizeof(R));
      memcpy(data->_beta,beta, (size_t)((set->N + 1)) * sizeof(R));
      memcpy(data->_gamma,gam, (size_t)((set->N + 1)) * sizeof(R));
    }
  }
}

void fpt_trafo_direct(fpt_set set, const int m, const C *x, C *y,
  const int k_end, const unsigned int flags)
{
  int j;
  fpt_data *data = &(set->dpt[m]);
  int Nk;
  int tk;
  R norm;
  
    //fprintf(stderr, "Executing dpt.\n");  

  Y(next_power_of_2_exp)(k_end+1,&Nk,&tk);
  norm = 2.0/(Nk<<1);

    //fprintf(stderr, "Norm = %e.\n", norm);  

  if (set->flags & FPT_NO_DIRECT_ALGORITHM)
  {
    return;
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
    /* Fill array with Chebyshev nodes. */
    for (j = 0; j <= k_end; j++)
    {
      set->xc_slow[j] = cos((KPI*(j+0.5))/(k_end+1));
        //fprintf(stderr, "x[%4d] = %e.\n", j, set->xc_slow[j]);  
    }

    memset(set->result,0U,data->k_start*sizeof(C));
    memcpy(&set->result[data->k_start],x, (size_t)((k_end - data->k_start + 1)) * sizeof(C));

    /*eval_sum_clenshaw(k_end, k_end, set->result, set->xc_slow,
      y, set->work, &data->alpha[1], &data->beta[1], &data->gamma[1],
      data->gamma_m1);*/
    eval_sum_clenshaw_fast(k_end, k_end, set->result, set->xc_slow,
      y, &data->_alpha[1], &data->_beta[1], &data->_gamma[1], data->gamma_m1);
  }
  else
  {
    memset(set->temp,0U,data->k_start*sizeof(C));
    memcpy(&set->temp[data->k_start],x, (size_t)((k_end - data->k_start + 1)) * sizeof(C));

    eval_sum_clenshaw_fast(k_end, Nk-1, set->temp, set->xcvecs[tk-2],
      set->result, &data->_alpha[1], &data->_beta[1], &data->_gamma[1],
      data->gamma_m1);

    fftw_execute_r2r(set->plans_dct2[tk-2],(R*)set->result,
      (R*)set->result);

    set->result[0] *= 0.5;
    for (j = 0; j < Nk; j++)
    {
      set->result[j] *= norm;
    }

    memcpy(y,set->result, (size_t)((k_end + 1)) * sizeof(C));
  }
}

void fpt_trafo(fpt_set set, const int m, const C *x, C *y,
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
  fftw_plan plan = 0;
  int length = k_end+1;
  fftw_r2r_kind kinds[2] = {FFTW_REDFT01,FFTW_REDFT01};

  /** Loop counter */
  int k;

  C *work_ptr;
  const C *x_ptr;

  /* Check, if slow transformation should be used due to small bandwidth. */
  if (k_end < FPT_BREAK_EVEN)
  {
    /* Use NDSFT. */
    fpt_trafo_direct(set, m, x, y, k_end, flags);
    return;
  }

  Y(next_power_of_2_exp)(k_end,&Nk,&tk);
  k_start_tilde = K_START_TILDE(data->k_start,Nk);
  k_end_tilde = K_END_TILDE(k_end,Nk);

  /* Check if fast transform is activated. */
  if (set->flags & FPT_NO_FAST_ALGORITHM)
    return;

  if (flags & FPT_FUNCTION_VALUES)
  {
#ifdef _OPENMP
    int nthreads = X(get_num_threads)();
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    fftw_plan_with_nthreads(nthreads);
#endif
    plan = fftw_plan_many_r2r(1, &length, 2, (R*)set->work, NULL, 2, 1,
      (R*)set->work, NULL, 2, 1, kinds, 0U);
#ifdef _OPENMP
}
#endif
  }

  /* Initialize working arrays. */
  memset(set->result,0U,2*Nk*sizeof(C));

  /* The first step. */

  /* Set the first 2*data->k_start coefficients to zero. */
  memset(set->work,0U,2*data->k_start*sizeof(C));

  work_ptr = &set->work[2*data->k_start];
  x_ptr = x;

  for (k = 0; k <= k_end_tilde-data->k_start; k++)
  {
    *work_ptr++ = *x_ptr++;
    *work_ptr++ = K(0.0);
  }

  /* Set the last 2*(set->N-1-k_end_tilde) coefficients to zero. */
  memset(&set->work[2*(k_end_tilde+1)],0U,2*(Nk-1-k_end_tilde)*sizeof(C));

  /* If k_end == Nk, use three-term recurrence to map last coefficient x_{Nk} to
   * x_{Nk-1} and x_{Nk-2}. */
  if (k_end == Nk)
  {
    set->work[2*(Nk-2)] += data->gammaN[tk-2]*x[Nk-data->k_start];
    set->work[2*(Nk-1)] += data->betaN[tk-2]*x[Nk-data->k_start];
    set->work[2*(Nk-1)+1] = data->alphaN[tk-2]*x[Nk-data->k_start];
  }

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
      memcpy(set->vec3,&(set->work[(plength/2)*(4*l+2)]), (size_t)((plength / 2)) * sizeof(C));
      memcpy(set->vec4,&(set->work[(plength/2)*(4*l+3)]), (size_t)((plength / 2)) * sizeof(C));
      memset(&set->vec3[plength/2],0U,(plength/2)*sizeof(C));
      memset(&set->vec4[plength/2],0U,(plength/2)*sizeof(C));

      /* Copy coefficients into first half. */
      memcpy(&(set->work[(plength/2)*(4*l+2)]),&(set->work[(plength/2)*(4*l+1)]), (size_t)((plength / 2)) * sizeof(C));
      memset(&(set->work[(plength/2)*(4*l+1)]),0U,(plength/2)*sizeof(C));
      memset(&(set->work[(plength/2)*(4*l+3)]),0U,(plength/2)*sizeof(C));

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
            step->a12[0], step->a21[0], step->a22[0], step->g[0], tau, set);
        }
        else
        {
          /* Multiply third and fourth polynomial with matrix U. */
          fpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->g[0], tau, set);
        }

        if (step->g[0] != 0.0)
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
        plength_stab = step->Ns;
        t_stab = step->ts;

        /*---------*/
        /*fprintf(stderr,"\nfpt_trafo: stabilizing at tau = %d, l = %d.\n",tau,l);
        fprintf(stderr,"\nfpt_trafo: plength_stab = %d.\n",plength_stab);
        fprintf(stderr,"\nfpt_trafo: tk = %d.\n",tk);
        fprintf(stderr,"\nfpt_trafo: index = %d.\n",tk-tau-1);*/
        /*---------*/

        /* Set rest of vectors explicitely to zero */
        /*fprintf(stderr,"fpt_trafo: stabilizing: plength = %d, plength_stab = %d\n",
          plength, plength_stab);*/
        memset(&set->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(C));
        memset(&set->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(C));

        /* Multiply third and fourth polynomial with matrix U. */
        /* Check for symmetry. */
        if (set->flags & FPT_AL_SYMMETRY)
        {
          if (m <= 1)
          {
            fpt_do_step_symmetric(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->g[0], t_stab-1, set);
          }
          else if (m%2 == 0)
          {
            /*fpt_do_step_symmetric_u(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);*/
            fpt_do_step_symmetric_u(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0],
              set->xcvecs[t_stab-2], step->g[0], t_stab-1, set);
            /*fpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);*/
          }
          else
          {
            /*fpt_do_step_symmetric_l(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);*/
            fpt_do_step_symmetric_l(set->vec3, set->vec4,
              step->a11[0], step->a12[0],
              step->a21[0],
              step->a22[0], set->xcvecs[t_stab-2], step->g[0], t_stab-1, set);
            /*fpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->gamma[0], t_stab-1, set);*/
          }
        }
        else
        {
            fpt_do_step(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->g[0], t_stab-1, set);
        }

        if (step->g[0] != 0.0)
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
    fftw_execute_r2r(plan,(R*)y,(R*)y);
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    fftw_destroy_plan(plan);
    for (k = 0; k <= k_end; k++)
    {
      y[k] *= 0.5;
    }
  }
}

void fpt_transposed_direct(fpt_set set, const int m, C *x,
  C *y, const int k_end, const unsigned int flags)
{
  int j;
  fpt_data *data = &(set->dpt[m]);
  int Nk;
  int tk;
  R norm;

  Y(next_power_of_2_exp)(k_end+1,&Nk,&tk);
  norm = 2.0/(Nk<<1);

  if (set->flags & FPT_NO_DIRECT_ALGORITHM)
  {
    return;
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
    for (j = 0; j <= k_end; j++)
    {
      set->xc_slow[j] = cos((KPI*(j+0.5))/(k_end+1));
    }

    eval_sum_clenshaw_transposed(k_end, k_end, set->result, set->xc_slow,
      y, set->work, &data->_alpha[1], &data->_beta[1], &data->_gamma[1],
      data->gamma_m1);

    memcpy(x,&set->result[data->k_start], (size_t)((k_end - data->k_start + 1)) *
      sizeof(C));
  }
  else
  {
    memcpy(set->result,y, (size_t)((k_end + 1)) * sizeof(C));
    memset(&set->result[k_end+1],0U,(Nk-k_end-1)*sizeof(C));

    for (j = 0; j < Nk; j++)
    {
      set->result[j] *= norm;
    }

    fftw_execute_r2r(set->plans_dct3[tk-2],(R*)set->result,
      (R*)set->result);

    eval_sum_clenshaw_transposed(k_end, Nk-1, set->temp, set->xcvecs[tk-2],
      set->result, set->work, &data->_alpha[1], &data->_beta[1], &data->_gamma[1],
      data->gamma_m1);

    memcpy(x,&set->temp[data->k_start], (size_t)((k_end - data->k_start + 1)) * sizeof(C));
  }
}

void fpt_transposed(fpt_set set, const int m, C *x,
  C *y, const int k_end, const unsigned int flags)
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
    fpt_transposed_direct(set, m, x, y, k_end, flags);
    return;
  }

  Y(next_power_of_2_exp)(k_end,&Nk,&tk);
  k_start_tilde = K_START_TILDE(data->k_start,Nk);
  k_end_tilde = K_END_TILDE(k_end,Nk);

  /* Check if fast transform is activated. */
  if (set->flags & FPT_NO_FAST_ALGORITHM)
  {
    return;
  }

  if (flags & FPT_FUNCTION_VALUES)
  {
#ifdef _OPENMP
    int nthreads = X(get_num_threads)();
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    fftw_plan_with_nthreads(nthreads);
#endif
    plan = fftw_plan_many_r2r(1, &length, 2, (R*)set->work, NULL, 2, 1,
      (R*)set->work, NULL, 2, 1, kinds, 0U);
#ifdef _OPENMP
}
#endif
    fftw_execute_r2r(plan,(R*)y,(R*)set->result);
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    fftw_destroy_plan(plan);
    for (k = 0; k <= k_end; k++)
    {
      set->result[k] *= 0.5;
    }
  }
  else
  {
    memcpy(set->result,y, (size_t)((k_end + 1)) * sizeof(C));
  }

  /* Initialize working arrays. */
  memset(set->work,0U,2*Nk*sizeof(C));

  /* The last step is now the first step. */
  for (k = 0; k <= k_end; k++)
  {
    set->work[k] = data->gamma_m1*set->result[k];
  }
  //memset(&set->work[k_end+1],0U,(Nk+1-k_end)*sizeof(C));

  set->work[Nk] = data->gamma_m1*(data->beta_0*set->result[0] +
    data->alpha_0*set->result[1]);
  for (k = 1; k < k_end; k++)
  {
    set->work[Nk+k] = data->gamma_m1*(data->beta_0*set->result[k] +
      data->alpha_0*0.5*(set->result[k-1]+set->result[k+1]));
  }
  if (k_end<Nk)
  {
    memset(&set->work[k_end],0U,(Nk-k_end)*sizeof(C));
  }

  /** Save copy of inpute data for stabilization steps. */
  memcpy(set->result,set->work, (size_t)(2 * Nk) * sizeof(C));

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
      memcpy(set->vec3,&(set->work[(plength/2)*(4*l+0)]), (size_t)(plength) * sizeof(C));
      memcpy(set->vec4,&(set->work[(plength/2)*(4*l+2)]), (size_t)(plength) * sizeof(C));

      memcpy(&set->work[(plength/2)*(4*l+1)],&(set->work[(plength/2)*(4*l+2)]),
        (size_t)((plength / 2)) * sizeof(C));

      /* Get matrix U_{n,tau,l} */
      step = &(data->steps[tau][l]);

      /* Check if step is stable. */
      if (step->stable)
      {
        if (set->flags & FPT_AL_SYMMETRY && IS_SYMMETRIC(l,m,plength))
        {
          /* Multiply third and fourth polynomial with matrix U. */
          fpt_do_step_t_symmetric(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->g[0], tau, set);
        }
        else
        {
          /* Multiply third and fourth polynomial with matrix U. */
          fpt_do_step_t(set->vec3, set->vec4, step->a11[0], step->a12[0],
            step->a21[0], step->a22[0], step->g[0], tau, set);
        }
        memcpy(&(set->vec3[plength/2]), set->vec4, (size_t)((plength / 2)) * sizeof(C));

        for (k = 0; k < plength; k++)
        {
          set->work[plength*(4*l+2)/2+k] = set->vec3[k];
        }
      }
      else
      {
        /* Stabilize. */
        plength_stab = step->Ns;
        t_stab = step->ts;

        memcpy(set->vec3,set->result, (size_t)(plength_stab) * sizeof(C));
        memcpy(set->vec4,&(set->result[Nk]), (size_t)(plength_stab) * sizeof(C));

        /* Multiply third and fourth polynomial with matrix U. */
        if (set->flags & FPT_AL_SYMMETRY)
        {
          if (m <= 1)
          {
            fpt_do_step_t_symmetric(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->g[0], t_stab-1, set);
          }
          else if (m%2 == 0)
          {
            fpt_do_step_t_symmetric_u(set->vec3, set->vec4, step->a11[0], step->a12[0],
              set->xcvecs[t_stab-2], step->g[0], t_stab-1, set);
          }
          else
          {
            fpt_do_step_t_symmetric_l(set->vec3, set->vec4,
              step->a21[0], step->a22[0], set->xcvecs[t_stab-2], step->g[0], t_stab-1, set);
          }
        }
        else
        {
            fpt_do_step_t(set->vec3, set->vec4, step->a11[0], step->a12[0],
              step->a21[0], step->a22[0], step->g[0], t_stab-1, set);
        }

        memcpy(&(set->vec3[plength/2]),set->vec4, (size_t)((plength / 2)) * sizeof(C));

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
  const int M = set->M;

  /* TODO Clean up DPT transform data structures. */
  for (m = 0; m < M; m++)
  {
    /* Check if precomputed. */
    data = &set->dpt[m];
    if (data->steps != (fpt_step**)NULL)
    {
      nfft_free(data->alphaN);
      nfft_free(data->betaN);
      nfft_free(data->gammaN);
      data->alphaN = NULL;
      data->betaN = NULL;
      data->gammaN = NULL;

      /* Free precomputed data. */
      k_start_tilde = K_START_TILDE(data->k_start,Y(next_power_of_2)(data->k_start)
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
          nfft_free(data->steps[tau][l].a11[0]);
          nfft_free(data->steps[tau][l].a12[0]);
          nfft_free(data->steps[tau][l].a21[0]);
          nfft_free(data->steps[tau][l].a22[0]);
          data->steps[tau][l].a11[0] = NULL;
          data->steps[tau][l].a12[0] = NULL;
          data->steps[tau][l].a21[0] = NULL;
          data->steps[tau][l].a22[0] = NULL;
          /* Free components. */
          nfft_free(data->steps[tau][l].a11);
          nfft_free(data->steps[tau][l].a12);
          nfft_free(data->steps[tau][l].a21);
          nfft_free(data->steps[tau][l].a22);
          nfft_free(data->steps[tau][l].g);
          data->steps[tau][l].a11 = NULL;
          data->steps[tau][l].a12 = NULL;
          data->steps[tau][l].a21 = NULL;
          data->steps[tau][l].a22 = NULL;
          data->steps[tau][l].g = NULL;
        }
        /* Free pointers for current level. */
        nfft_free(data->steps[tau]);
        data->steps[tau] = NULL;
        /* Double length of polynomials. */
        plength = plength<<1;
      }
      /* Free steps. */
      nfft_free(data->steps);
      data->steps = NULL;
    }

    if (!(set->flags & FPT_NO_DIRECT_ALGORITHM))
    {
      /* Check, if recurrence coefficients must be copied. */
      //fprintf(stderr,"\nfpt_finalize: %d\n",set->flags & FPT_PERSISTENT_DATA);
      if (!(set->flags & FPT_PERSISTENT_DATA))
      {
        nfft_free(data->_alpha);
        nfft_free(data->_beta);
        nfft_free(data->_gamma);
      }
      data->_alpha = NULL;
      data->_beta = NULL;
      data->_gamma = NULL;
    }
  }

  /* Delete array of DPT transform data. */
  nfft_free(set->dpt);
  set->dpt = NULL;

  for (tau = 1; tau < set->t+1; tau++)
  {
    nfft_free(set->xcvecs[tau-1]);
    set->xcvecs[tau-1] = NULL;
  }
  nfft_free(set->xcvecs);
  set->xcvecs = NULL;

  /* Free auxilliary arrays. */
  nfft_free(set->work);
  nfft_free(set->result);

  /* Check if fast transform is activated. */
  if (!(set->flags & FPT_NO_FAST_ALGORITHM))
  {
    /* Free auxilliary arrays. */
    nfft_free(set->vec3);
    nfft_free(set->vec4);
    nfft_free(set->z);
    set->work = NULL;
    set->result = NULL;
    set->vec3 = NULL;
    set->vec4 = NULL;
    set->z = NULL;

    /* Free FFTW plans. */
    for(tau = 0; tau < set->t/*-1*/; tau++)
    {
#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
{
      fftw_destroy_plan(set->plans_dct3[tau]);
      fftw_destroy_plan(set->plans_dct2[tau]);
}
      set->plans_dct3[tau] = NULL;
      set->plans_dct2[tau] = NULL;
    }

    nfft_free(set->plans_dct3);
    nfft_free(set->plans_dct2);
    set->plans_dct3 = NULL;
    set->plans_dct2 = NULL;
  }

  if (!(set->flags & FPT_NO_DIRECT_ALGORITHM))
  {
    /* Delete arrays of Chebyshev nodes. */
    nfft_free(set->xc_slow);
    set->xc_slow = NULL;
    nfft_free(set->temp);
    set->temp = NULL;
  }

  /* Free DPT set structure. */
  nfft_free(set);
}
