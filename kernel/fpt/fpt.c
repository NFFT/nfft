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

/**
 * \file fpt.c
 * \brief Implementation file for the FPT module
 * \author Jens Keiner
 */

#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

/**
 * If defined, perform critical computations in three-term recurrence
 * always in long double instead of using adaptive version that starts
 * in double and switches to long double if required.
 */
#undef  FPT_CLENSHAW_USE_ONLY_LONG_DOUBLE

/* Macros for index calculation. */

/** Minimum degree at top of a cascade */
#define K_START_TILDE(x,y) (MAX(MIN(x,y-2),0))

/** Maximum degree at top of a cascade */
#define K_END_TILDE(x,y) MIN(x,y-1)

/** Index of first block of four functions at level */
#define FIRST_L(x,y) (LRINT(floor((x)/(double)y)))

/** Index of last block of four functions at level */
#define LAST_L(x,y) (LRINT(ceil(((x)+1)/(double)y))-1)

#define N_TILDE(y) (y-1)

#define IS_SYMMETRIC(x,y,z) (x >= ((y-1.0)/z))

#define FPT_BREAK_EVEN 4

#include "fpt.h"


static inline void abuvxpwy(double a, double b, double _Complex* u, double _Complex* x,
  double* v, double _Complex* y, double* w, int n)
{
  int l; double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y;
  double *v_ptr = v, *w_ptr = w;
  for (l = 0; l < n; l++)
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++));
}

#define ABUVXPWY_SYMMETRIC(NAME,S1,S2) \
static inline void NAME(double a, double b, double _Complex* u, double _Complex* x, \
  double* v, double _Complex* y, double* w, int n) \
{ \
  const int n2 = n>>1; \
  int l; double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y; \
  double *v_ptr = v, *w_ptr = w; \
  for (l = 0; l < n2; l++) \
    *u_ptr++ = a * (b * (*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++)); \
  v_ptr--; w_ptr--; \
  for (l = 0; l < n2; l++) \
    *u_ptr++ = a * (b * S1 * (*v_ptr--) * (*x_ptr++) + S2 * (*w_ptr--) * (*y_ptr++)); \
}

ABUVXPWY_SYMMETRIC(abuvxpwy_symmetric1,1.0,-1.0)
ABUVXPWY_SYMMETRIC(abuvxpwy_symmetric2,-1.0,1.0)

#define ABUVXPWY_SYMMETRIC_1(NAME,S1) \
static inline void NAME(double a, double b, double _Complex* u, double _Complex* x, \
  double* v, double _Complex* y, int n, double *xx) \
{ \
  const int n2 = n>>1; \
  int l; double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y; \
  double *v_ptr = v, *xx_ptr = xx; \
  for (l = 0; l < n2; l++, v_ptr++) \
    *u_ptr++ = a * (b * (*v_ptr) * (*x_ptr++) + ((*v_ptr)*(1.0+*xx_ptr++)) * (*y_ptr++)); \
  v_ptr--; \
  for (l = 0; l < n2; l++, v_ptr--) \
    *u_ptr++ = a * (b * S1 * (*v_ptr) * (*x_ptr++) + (S1 * (*v_ptr) * (1.0+*xx_ptr++)) * (*y_ptr++)); \
}

ABUVXPWY_SYMMETRIC_1(abuvxpwy_symmetric1_1,1.0)
ABUVXPWY_SYMMETRIC_1(abuvxpwy_symmetric1_2,-1.0)

#define ABUVXPWY_SYMMETRIC_2(NAME,S1) \
static inline void NAME(double a, double b, double _Complex* u, double _Complex* x, \
  double _Complex* y, double* w, int n, double *xx) \
{ \
  const int n2 = n>>1; \
  int l; double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y; \
  double *w_ptr = w, *xx_ptr = xx; \
  for (l = 0; l < n2; l++, w_ptr++) \
    *u_ptr++ = a * (b * (*w_ptr/(1.0+*xx_ptr++)) * (*x_ptr++) + (*w_ptr) * (*y_ptr++)); \
  w_ptr--; \
  for (l = 0; l < n2; l++, w_ptr--) \
    *u_ptr++ = a * (b * (S1 * (*w_ptr)/(1.0+*xx_ptr++) ) * (*x_ptr++) + S1 * (*w_ptr) * (*y_ptr++)); \
}

ABUVXPWY_SYMMETRIC_2(abuvxpwy_symmetric2_1,1.0)
ABUVXPWY_SYMMETRIC_2(abuvxpwy_symmetric2_2,-1.0)

static inline void auvxpwy(double a, double _Complex* u, double _Complex* x, double* v,
  double _Complex* y, double* w, int n)
{
  int l;
  double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y;
  double *v_ptr = v, *w_ptr = w;
  for (l = n; l > 0; l--)
    *u_ptr++ = a * ((*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++));
}

static inline void auvxpwy_symmetric(double a, double _Complex* u, double _Complex* x,
  double* v, double _Complex* y, double* w, int n)
{
  const int n2 = n>>1; \
  int l;
  double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y;
  double *v_ptr = v, *w_ptr = w;
  for (l = n2; l > 0; l--)
    *u_ptr++ = a * ((*v_ptr++) * (*x_ptr++) + (*w_ptr++) * (*y_ptr++));
  v_ptr--; w_ptr--;
  for (l = n2; l > 0; l--)
    *u_ptr++ = a * ((*v_ptr--) * (*x_ptr++) - (*w_ptr--) * (*y_ptr++));
}

static inline void auvxpwy_symmetric_1(double a, double _Complex* u, double _Complex* x,
  double* v, double _Complex* y, double* w, int n, double *xx)
{
  const int n2 = n>>1; \
  int l;
  double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y;
  double *v_ptr = v, *w_ptr = w, *xx_ptr = xx;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (((*v_ptr++)*(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr++)*(1.0+*xx_ptr)) * (*y_ptr++));
  v_ptr--; w_ptr--;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (-((*v_ptr--)*(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr--)*(1.0+*xx_ptr)) * (*y_ptr++));
}

static inline void auvxpwy_symmetric_2(double a, double _Complex* u, double _Complex* x,
  double* v, double _Complex* y, double* w, int n, double *xx)
{
  const int n2 = n>>1; \
  int l;
  double _Complex *u_ptr = u, *x_ptr = x, *y_ptr = y;
  double *v_ptr = v, *w_ptr = w, *xx_ptr = xx;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (((*v_ptr++)/(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr++)/(1.0+*xx_ptr)) * (*y_ptr++));
  v_ptr--; w_ptr--;
  for (l = n2; l > 0; l--, xx_ptr++)
    *u_ptr++ = a * (-((*v_ptr--)/(1.0+*xx_ptr)) * (*x_ptr++) + ((*w_ptr--)/(1.0+*xx_ptr)) * (*y_ptr++));
}

#define FPT_DO_STEP(NAME,M1_FUNCTION,M2_FUNCTION) \
static inline void NAME(double _Complex  *a, double _Complex *b, double *a11, double *a12, \
  double *a21, double *a22, double g, int tau, fpt_set set) \
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
    memcpy(b,set->z,length*sizeof(double _Complex)); \
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
/*FPT_DO_STEP(fpt_do_step_symmetric_u,auvxpwy_symmetric,auvxpwy)
FPT_DO_STEP(fpt_do_step_symmetric_l,auvxpwy,auvxpwy_symmetric)*/

static inline void fpt_do_step_symmetric_u(double _Complex *a, double _Complex *b,
  double *a11, double *a12, double *a21, double *a22, double *x,
  double gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  double norm = 1.0/(length<<1);

  UNUSED(a21); UNUSED(a22);

  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)b,(double*)b);

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
    memcpy(b,set->z,length*sizeof(double _Complex));
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a);
    /* Compensate for factors introduced by a raw DCT-II. */
    a[0] *= 0.5;
  }

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b);
  /* Compensate for factors introduced by a raw DCT-II. */
  b[0] *= 0.5;
}

static inline void fpt_do_step_symmetric_l(double _Complex  *a, double _Complex *b,
  double *a11, double *a12, double *a21, double *a22, double *x, double gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  double norm = 1.0/(length<<1);

  /* Compensate for factors introduced by a raw DCT-III. */
  a[0] *= 2.0;
  b[0] *= 2.0;

  UNUSED(a11); UNUSED(a12);

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)b,(double*)b);

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
    memcpy(b,set->z,length*sizeof(double _Complex));
    /* Compute Chebyshev-coefficients using a DCT-II. */
    fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a);
    /* Compensate for factors introduced by a raw DCT-II. */
    a[0] *= 0.5;
  }

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b);
  /* Compensate for factors introduced by a raw DCT-II. */
  b[0] *= 0.5;
}

#define FPT_DO_STEP_TRANSPOSED(NAME,M1_FUNCTION,M2_FUNCTION) \
static inline void NAME(double _Complex  *a, double _Complex *b, double *a11, \
  double *a12, double *a21, double *a22, double g, int tau, fpt_set set) \
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
  M1_FUNCTION(norm,g,set->z,a,a11,b,a21,length); \
  M2_FUNCTION(norm,g,b,a,a12,b,a22,length); \
  memcpy(a,set->z,length*sizeof(double _Complex)); \
  \
  /* Compute Chebyshev-coefficients using a DCT-II. */ \
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a); \
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b); \
}

FPT_DO_STEP_TRANSPOSED(fpt_do_step_t,abuvxpwy,abuvxpwy)
FPT_DO_STEP_TRANSPOSED(fpt_do_step_t_symmetric,abuvxpwy_symmetric1,abuvxpwy_symmetric2)
/*FPT_DO_STEP_TRANSPOSED(fpt_do_step_t_symmetric_u,abuvxpwy_symmetric1_1,abuvxpwy_symmetric1_2)*/
/*FPT_DO_STEP_TRANSPOSED(fpt_do_step_t_symmetric_l,abuvxpwy_symmetric2_2,abuvxpwy_symmetric2_1)*/


static inline void fpt_do_step_t_symmetric_u(double _Complex  *a,
  double _Complex *b, double *a11, double *a12, double *x,
  double gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  double norm = 1.0/(length<<1);

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)b,(double*)b);

  /* Perform matrix multiplication. */
  abuvxpwy_symmetric1_1(norm,gam,set->z,a,a11,b,length,x);
  abuvxpwy_symmetric1_2(norm,gam,b,a,a12,b,length,x);
  memcpy(a,set->z,length*sizeof(double _Complex));

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b);
}

static inline void fpt_do_step_t_symmetric_l(double _Complex  *a,
  double _Complex *b, double *a21, double *a22,
  double *x, double gam, int tau, fpt_set set)
{
  /** The length of the coefficient arrays. */
  int length = 1<<(tau+1);
  /** Twice the length of the coefficient arrays. */
  double norm = 1.0/(length<<1);

  /* Compute function values from Chebyshev-coefficients using a DCT-III. */
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(set->plans_dct3[tau-1],(double*)b,(double*)b);

  /* Perform matrix multiplication. */
  abuvxpwy_symmetric2_2(norm,gam,set->z,a,b,a21,length,x);
  abuvxpwy_symmetric2_1(norm,gam,b,a,b,a22,length,x);
  memcpy(a,set->z,length*sizeof(double _Complex));

  /* Compute Chebyshev-coefficients using a DCT-II. */
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)a,(double*)a);
  fftw_execute_r2r(set->plans_dct2[tau-1],(double*)b,(double*)b);
}


static void eval_clenshaw(const double *x, double *y, int size, int k, const double *alpha,
  const double *beta, const double *gam)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  const double *x_act;
  double *y_act;
  const double *alpha_act, *beta_act, *gamma_act;

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  for (i = 0; i < size; i++)
  {
    double x_val_act = *x_act;

    if (k == 0)
    {
      *y_act = 1.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gam[k]);

#ifdef FPT_CLENSHAW_USE_ONLY_LONG_DOUBLE
      long double a = 1.0;
      long double b = 0.0;
      for (j = k; j > 1; j--)
      {
        long double a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
        b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
#else
      double a = 1.0;
      double b = 0.0;
      /* 1e247 should not trigger for NFSFT with N <= 1024 */
      for (j = k; j > 1 && fabs(a) < 1e247; j--)
      {
        double a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
        b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      if (j <= 1)
        *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
      else /* fabs(a) >= 1e247, continue in long double */
      {
        long double a_ld = a;
        long double b_ld = b;
        for (; j > 1; j--)
        {
          long double a_old = a_ld;
          a_ld = b_ld + a_old*((*alpha_act)*x_val_act+(*beta_act));
          b_ld = a_old*(*gamma_act);
          alpha_act--;
          beta_act--;
          gamma_act--;
        }
        *y_act = (a_ld*((*alpha_act)*x_val_act+(*beta_act))+b_ld);
      }
#endif

    }
    x_act++;
    y_act++;
  }
}

static void eval_clenshaw2(const double *x, double *z, double *y, int size1, int size, int k, const double *alpha,
  const double *beta, const double *gam)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  double a,b,x_val_act,a_old;
  const double *x_act;
  double *y_act, *z_act;
  const double *alpha_act, *beta_act, *gamma_act;

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

static int eval_clenshaw_thresh2(const double *x, double *z, double *y, int size, int k,
  const double *alpha, const double *beta, const double *gam, const
  double threshold)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  double x_val_act;
  const double *x_act;
  double *y_act, *z_act;
  const double *alpha_act, *beta_act, *gamma_act;
  const R threshold_abs = FABS(threshold);

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  z_act = z;
  for (i = 0; i < size; i++)
  {
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

#ifdef FPT_CLENSHAW_USE_ONLY_LONG_DOUBLE
      long double a = 1.0;
      long double b = 0.0;
      for (j = k; j > 1; j--)
      {
        long double a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
        b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *z_act = a;
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
#else
      double a = 1.0;
      double b = 0.0;
      for (j = k; j > 1 && fabs(a) < 1e247; j--)
      {
        double a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
        b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      if (j <= 1)
      {
        *z_act = a;
        *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
      }
      else /* fabs(a) >= 1e247, continue in long double */
      {
        long double a_ld = a;
        long double b_ld = b;
        for (; j > 1; j--)
        {
          long double a_old = a_ld;
          a_ld = b_ld + a_old*((*alpha_act)*x_val_act+(*beta_act));
          b_ld = a_old*(*gamma_act);
          alpha_act--;
          beta_act--;
          gamma_act--;
        }
        *z_act = a_ld;
        *y_act = (a_ld*((*alpha_act)*x_val_act+(*beta_act))+b_ld);
      }
#endif
      if (FABS(*y_act) > threshold_abs)
        return 1;
    }
    x_act++;
    y_act++;
    z_act++;
  }
  return 0;
}

static inline void eval_sum_clenshaw_fast(const int N, const int M,
  const double _Complex *a, const double *x, double _Complex *y,
  const double *alpha, const double *beta, const double *gam,
  const double lambda)
{
  int j,k;
  
  if (N == 0)
    for (j = 0; j <= M; j++)
      y[j] = a[0];
  else
  {
    for (j = 0; j <= M; j++)
    {
      double xc = x[j];
#ifdef FPT_CLENSHAW_USE_ONLY_LONG_DOUBLE
      long double _Complex tmp1 = a[N-1];
      long double _Complex tmp2 = a[N];
      for (k = N-1; k > 0; k--)
      {
        long double _Complex tmp3 = a[k-1] + tmp2 * gam[k];
        tmp2 *= (alpha[k] * xc + beta[k]);
        tmp2 += tmp1;
        tmp1 = tmp3;
      }
      tmp2 *= (alpha[0] * xc + beta[0]);
      y[j] = lambda * (tmp2 + tmp1);
#else
      double _Complex tmp1 = a[N-1];
      double _Complex tmp2 = a[N];
      /* 1e247 should not trigger for NFSFT with N <= 1024 */
      for (k = N-1; k > 0 && fabs(creal(tmp2)) < 1e247 && fabs(cimag(tmp2)) < 1e247; k--)
      {
        double _Complex tmp3 = a[k-1] + tmp2 * gam[k];
        tmp2 *= (alpha[k] * xc + beta[k]);
        tmp2 += tmp1;
        tmp1 = tmp3;
      }
      if (k <= 0)
      {
        tmp2 *= (alpha[0] * xc + beta[0]);
        y[j] = lambda * (tmp2 + tmp1);
      }
      else /* fabs(tmp2) >= 1e247 */
      {
        long double _Complex tmp1_ld = tmp1;
        long double _Complex tmp2_ld = tmp2;
        for (; k > 0; k--)
        {
          long double _Complex tmp3_ld = a[k-1] + tmp2_ld * gam[k];
          tmp2_ld *= (alpha[k] * xc + beta[k]);
          tmp2_ld += tmp1_ld;
          tmp1_ld = tmp3_ld;
        }
        tmp2_ld *= (alpha[0] * xc + beta[0]);
        y[j] = lambda * (tmp2_ld + tmp1_ld);
      } /* end fabs(tmp2) >= 1e247 */
#endif
    } /* for j */
  } /* N > 0 */
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
 * for given double _Complex coefficients \f$\left(a_k\right)_{k=0}^N \in
 * \mathbb{C}^{N+1}\f$ at given nodes \f$\left(x_j\right)_{j=0}^M \in
 * \mathbb{R}^{M+1}\f$, \f$M \in \mathbb{N}_0\f$.
 */
static void eval_sum_clenshaw_transposed(int N, int M, double _Complex* a, double *x,
  double _Complex *y, double _Complex *temp, double *alpha, double *beta, double *gam,
  double lambda)
{
  int j,k;
  double _Complex* it1 = temp;
  double _Complex* it2 = y;
  double _Complex aux;

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

static void eval_sum_clenshaw_transposed_ld(int N, int M, double _Complex* a, double *x,
  double _Complex *y, double _Complex *temp, double *alpha, double *beta, double *gam,
  double lambda)
{
  int j,k;

  for (k = 0; k <= N; k++)
    a[k] = 0.0;

  if (N == 0)
    for (j = 0; j <= M; j++)
      a[0] += lambda * y[j];
  else
  {
    for (j = 0; j <= M; j++)
    {
      /* Compute final result by multiplying with the constant lambda */
      long double _Complex it2 = lambda * y[j];
      a[0] += it2;

      /* Compute final step. */
      long double _Complex it1 = it2;
      it2 = it2 * (alpha[0] * x[j] + beta[0]);
      a[1] += it2;

      for (k = 2; k <= N; k++)
      {
        long double _Complex aux = it1;
        it1 = it2;
        it2 = it2*(alpha[k-1] * x[j] + beta[k-1]) + gam[k-1] * aux;
        a[k] += it2;
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
  fpt_set_s *set = (fpt_set_s*)nfft_malloc(sizeof(fpt_set_s));

  /* Save parameters in structure. */
  set->flags = flags;

  set->M = M;
  set->t = t;
  set->N = 1<<t;

  if (!(flags & FPT_NO_INIT_FPT_DATA))
  {
    /* Allocate memory for M transforms. */
    set->dpt = (fpt_data*) nfft_malloc(M*sizeof(fpt_data));

    /* Initialize with NULL pointer. */
    for (m = 0; m < set->M; m++)
      {
      set->dpt[m].steps = NULL;
      set->dpt[m].precomputed = false;
      }
  }
  else
    set->dpt = NULL;

  /* Create arrays with Chebyshev nodes. */

  /* Initialize array with Chebyshev coefficients for the polynomial x. This
   * would be trivially an array containing a 1 as second entry with all other
   * coefficients set to zero. In order to compensate for the multiplicative
   * factor 2 introduced by the DCT-III, we set this coefficient to 0.5 here. */

  /* Allocate memory for array of pointers to node arrays. */
  set->xcvecs = (double**) nfft_malloc((set->t)*sizeof(double*));
  /* For each polynomial length starting with 4, compute the Chebyshev nodes
   * using a DCT-III. */
  plength = 4;
  for (tau = 1; tau < t+1; tau++)
  {
    /* Allocate memory for current array. */
    set->xcvecs[tau-1] = (double*) nfft_malloc(plength*sizeof(double));
    for (k = 0; k < plength; k++)
    {
      set->xcvecs[tau-1][k] = cos(((k+0.5)*KPI)/plength);
    }
    plength = plength << 1;
  }

  /** Allocate memory for auxilliary arrays. */
  set->work = (double _Complex*) nfft_malloc((2*set->N)*sizeof(double _Complex));
  set->result = (double _Complex*) nfft_malloc((2*set->N)*sizeof(double _Complex));

  set->plans_dct2 = (fftw_plan*) nfft_malloc(sizeof(fftw_plan)*(set->t/*-1*/));
  set->kindsr     = (fftw_r2r_kind*) nfft_malloc(2*sizeof(fftw_r2r_kind));
  set->kindsr[0]  = FFTW_REDFT10;
  set->kindsr[1]  = FFTW_REDFT10;
  for (tau = 0, plength = 4; tau < set->t/*-1*/; tau++, plength<<=1)
  {
#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    fftw_plan_with_nthreads(nthreads);
#endif
    set->plans_dct2[tau] =
      fftw_plan_many_r2r(1, &plength, 2, (double*)set->work, NULL,
                         2, 1, (double*)set->result, NULL, 2, 1,set->kindsr,
                         0);
#ifdef _OPENMP
}
#endif
  }

  /** Initialize FFTW plans. */
  set->plans_dct3 = (fftw_plan*) nfft_malloc(sizeof(fftw_plan)*(set->t/*-1*/));
  set->kinds      = (fftw_r2r_kind*) nfft_malloc(2*sizeof(fftw_r2r_kind));
  set->kinds[0]   = FFTW_REDFT01;
  set->kinds[1]   = FFTW_REDFT01;
  for (tau = 0, plength = 4; tau < set->t/*-1*/; tau++, plength<<=1)
  {
#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
  fftw_plan_with_nthreads(nthreads);
#endif
    set->plans_dct3[tau] =
      fftw_plan_many_r2r(1, &plength, 2, (double*)set->work, NULL,
                         2, 1, (double*)set->result, NULL, 2, 1, set->kinds,
                         0);
#ifdef _OPENMP
}
#endif
  }
  nfft_free(set->kinds);
  nfft_free(set->kindsr);
  set->kinds = NULL;
  set->kindsr = NULL;

  set->vec3 = NULL;
  set->vec4 = NULL;
  set->z = NULL;

  set->xc_slow = NULL;
  set->temp = NULL;

  /* Check if fast transform is activated. */
  if (!(set->flags & FPT_NO_FAST_ALGORITHM))
  {
    /** Allocate memory for auxilliary arrays. */
    set->vec3 = (double _Complex*) nfft_malloc(set->N*sizeof(double _Complex));
    set->vec4 = (double _Complex*) nfft_malloc(set->N*sizeof(double _Complex));
    set->z = (double _Complex*) nfft_malloc(set->N*sizeof(double _Complex));
  }

  if (!(set->flags & FPT_NO_DIRECT_ALGORITHM))
  {
    set->xc_slow = (double*) nfft_malloc((set->N+1)*sizeof(double));
    set->temp = (double _Complex*) nfft_malloc((set->N+1)*sizeof(double _Complex));

    if (!(flags & FPT_NO_INIT_FPT_DATA))
    {
      for (m = 0; m < set->M; m++)
      {
        fpt_data *data = &(set->dpt[m]);
        data->_alpha = NULL;
        data->_beta = NULL;
        data->_gamma = NULL;
      }
    }
  }

  /* Return the newly created DPT set. */
  return set;
}

void fpt_precompute_1(fpt_set set, const int m, int k_start)
{
  int tau;          /**< Cascade level                                       */
  int l;            /**< Level index                                         */
  int plength;      /**< Length of polynomials for the next level in the
                         cascade                                             */
  int degree;       /**< Degree of polynomials for the current level in the
                         cascade                                             */
  int firstl;       /**< First index l for current cascade level             */
  int lastl;        /**< Last index l for current cascade level and current  */
  int k_start_tilde;
  int N_tilde;
  int clength;
  fpt_data *data;

  /* Get pointer to DPT data. */
  data = &(set->dpt[m]);

  /* Check, if already precomputed. */
  if (data->steps != NULL)
    return;

  /* Save k_start. */
  data->k_start = k_start;

  data->alphaN = NULL;
  data->betaN = NULL;
  data->gammaN = NULL;

  if (!(set->flags & FPT_NO_FAST_ALGORITHM))
  {
    /* Save recursion coefficients. */
    data->alphaN = (double*) nfft_malloc(3*(set->t-1)*sizeof(double));
    data->betaN = data->alphaN + (set->t-1);
    data->gammaN = data->betaN + (set->t-1);

    k_start_tilde = K_START_TILDE(data->k_start,X(next_power_of_2)(data->k_start)
          /*set->N*/);
    N_tilde = N_TILDE(set->N);
  
    /* Allocate memory for the cascade with t = log_2(N) many levels. */
    data->steps = (fpt_step**) nfft_malloc(sizeof(fpt_step*)*set->t);

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
      data->steps[tau] = (fpt_step*) nfft_malloc(sizeof(fpt_step)
                         * (lastl+1));

      /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
      for (l = firstl; l <= lastl; l++)
      {
        if ((set->flags & FPT_AL_SYMMETRY) && IS_SYMMETRIC(l,m,plength))
        {
          clength = plength/2;
        }
        else
        {
          clength = plength;
        }

        /* Allocate memory for the components of U_{n,tau,l}. */
        data->steps[tau][l].a = (double*) nfft_malloc(sizeof(double)*clength*4);
      }
      /** Increase polynomial degree to next power of two. */
      plength = plength << 1;
    }
  }

  if (!(set->flags & FPT_NO_DIRECT_ALGORITHM) && !(set->flags & FPT_PERSISTENT_DATA) && (data->_alpha == NULL))
  {
    data->_alpha = (double*) nfft_malloc(3*(set->N+1)*sizeof(double));
    data->_beta = data->_alpha + (set->N+1);
    data->_gamma = data->_beta + (set->N+1);
  }
}

void fpt_precompute_2(fpt_set set, const int m, double *alpha, double *beta,
  double *gam, int k_start, const double threshold)
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
  /*double *a11;*/      /**< Array containing function values of the
                         (1,1)-component of U_k^n.                           */
  /*double *a12;*/      /**< Array containing function values of the
                         (1,2)-component of U_k^n.                           */
  /*double *a21;*/      /**< Array containing function values of the
                         (2,1)-component of U_k^n.                           */
  /*double *a22;*/      /**< Array containing function values of the
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

  /* Get pointer to DPT data. */
  data = &(set->dpt[m]);

  /* Check, if already precomputed. */
  if ((data->steps != NULL) && (data->precomputed))
    return;

  /* Save k_start. */
  data->k_start = k_start;

  data->gamma_m1 = gam[0];
/* moved to fpt_precompute_1
  data->alphaN = NULL;
  data->betaN = NULL;
  data->gammaN = NULL;*/

  /* Check if fast transform is activated. */
  if (!(set->flags & FPT_NO_FAST_ALGORITHM))
  {
    /* Save recursion coefficients. moved to fpt_precompute_1
    data->alphaN = (double*) nfft_malloc((set->t-1)*sizeof(double _Complex));
    data->betaN = (double*) nfft_malloc((set->t-1)*sizeof(double _Complex));
    data->gammaN = (double*) nfft_malloc((set->t-1)*sizeof(double _Complex)); */

    for (tau = 2; tau <= set->t; tau++)
    {

      data->alphaN[tau-2] = alpha[1<<tau];
      data->betaN[tau-2] = beta[1<<tau];
      data->gammaN[tau-2] = gam[1<<tau];
    }

    data->alpha_0 = alpha[1];
    data->beta_0 = beta[1];

    k_start_tilde = K_START_TILDE(data->k_start,X(next_power_of_2)(data->k_start)
      /*set->N*/);
    N_tilde = N_TILDE(set->N);

    /* Allocate memory for the cascade with t = log_2(N) many levels. moved to fpt_precompute_1
    data->steps = (fpt_step**) nfft_malloc(sizeof(fpt_step*)*set->t); */

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
       * many matrices. moved to fpt_precompute_1
      data->steps[tau] = (fpt_step*) nfft_malloc(sizeof(fpt_step)
                         * (lastl+1)); */

      /* For l = 0,...2^{t-tau-1}-1 compute the matrices U_{n,tau,l}. */
      for (l = firstl; l <= lastl; l++)
      {
        if ((set->flags & FPT_AL_SYMMETRY) && IS_SYMMETRIC(l,m,plength))
        {
          //fprintf(stderr,"fpt_precompute(%d): symmetric step\n",m);
          //fflush(stderr);
          clength = plength/2;
        }
        else
        {
          clength = plength;
        }

        /* Allocate memory for the components of U_{n,tau,l}. moved to fpt_precompute_1
        data->steps[tau][l].a11 = (double*) nfft_malloc(sizeof(double)*clength);
        data->steps[tau][l].a12 = (double*) nfft_malloc(sizeof(double)*clength);
        data->steps[tau][l].a21 = (double*) nfft_malloc(sizeof(double)*clength);
        data->steps[tau][l].a22 = (double*) nfft_malloc(sizeof(double)*clength); */

        /* Evaluate the associated polynomials at the 2^{tau+1} Chebyshev
         * nodes. */

        /* Get the pointers to the three-term recurrence coeffcients. */
        calpha = &(alpha[plength*l+1+1]);
        cbeta = &(beta[plength*l+1+1]);
        cgamma = &(gam[plength*l+1+1]);

        double *a11 = data->steps[tau][l].a;
        double *a12 = a11+clength;
        double *a21 = a12+clength;
        double *a22 = a21+clength;

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
          /* No stabilization needed. */
          data->steps[tau][l].g = gam[plength*l+1+1];
          data->steps[tau][l].stable = true;
        }
        else
        {
          /* Stabilize. */
          degree_stab = degree*(2*l+1);
          X(next_power_of_2_exp_int)((l+1)*(1<<(tau+1)),&N_stab,&t_stab);

          /* Old arrays are to small. */
          nfft_free(data->steps[tau][l].a);
          a11 = NULL;
          a12 = NULL;
          a21 = NULL;
          a22 = NULL;

          plength_stab = N_stab;

          if (set->flags & FPT_AL_SYMMETRY)
          {
            if (m <= 1)
            {
              /* This should never be executed */
              clength_1 = plength_stab;
              clength_2 = plength_stab;
              /* Allocate memory for arrays. */
              data->steps[tau][l].a = (double*) nfft_malloc(sizeof(double)*(clength_1*2+clength_2*2));
              a11 = data->steps[tau][l].a;
              a12 = a11+clength_1;
              a21 = a12+clength_1;
              a22 = a21+clength_2;
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
              data->steps[tau][l].a = (double*) nfft_malloc(sizeof(double)*clength*2);
              if (m%2 == 0)
              {
                a11 = data->steps[tau][l].a;
                a12 = a11+clength;
                calpha = &(alpha[2]); cbeta = &(beta[2]); cgamma = &(gam[2]);
                eval_clenshaw(set->xcvecs[t_stab-2], a11, clength,
                  degree_stab-2, calpha, cbeta, cgamma);
                eval_clenshaw(set->xcvecs[t_stab-2], a12, clength,
                  degree_stab-1, calpha, cbeta, cgamma);
              }
              else
              {
                a21 = data->steps[tau][l].a;
                a22 = a21+clength;
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
            data->steps[tau][l].a = (double*) nfft_malloc(sizeof(double)*(clength_1*2+clength_2*2));
            a11 = data->steps[tau][l].a;
            a12 = a11+clength_1;
            a21 = a12+clength_1;
            a22 = a21+clength_2;
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

          data->steps[tau][l].g =  gam[1+1];
          data->steps[tau][l].stable = false;
          data->steps[tau][l].ts = t_stab;
          data->steps[tau][l].Ns = N_stab;
        }
      }
      /** Increase polynomial degree to next power of two. */
      plength = plength << 1;
    }
    data->precomputed = true;
  }

  if (!(set->flags & FPT_NO_DIRECT_ALGORITHM))
  {
    /* Check, if recurrence coefficients must be copied. */
    if (set->flags & FPT_PERSISTENT_DATA)
    {
      data->_alpha = (double*) alpha;
      data->_beta = (double*) beta;
      data->_gamma = (double*) gam;
    }
    else
    {/* moved to fpt_precompute_1
      data->_alpha = (double*) nfft_malloc((set->N+1)*sizeof(double));
      data->_beta = (double*) nfft_malloc((set->N+1)*sizeof(double));
      data->_gamma = (double*) nfft_malloc((set->N+1)*sizeof(double));*/
      memcpy(data->_alpha,alpha,(set->N+1)*sizeof(double));
      memcpy(data->_beta,beta,(set->N+1)*sizeof(double));
      memcpy(data->_gamma,gam,(set->N+1)*sizeof(double));
    }
  }
}

void fpt_precompute(fpt_set set, const int m, double *alpha, double *beta,
  double *gam, int k_start, const double threshold)
{
  fpt_precompute_1(set, m, k_start);
  fpt_precompute_2(set, m, alpha, beta, gam, k_start, threshold);
}

void fpt_trafo_direct(fpt_set set, const int m, const double _Complex *x, double _Complex *y,
  const int k_end, const unsigned int flags)
{
  int j;
  fpt_data *data = &(set->dpt[m]);
  int Nk;
  int tk;
  double norm;
  
    //fprintf(stderr, "Executing dpt.\n");  

  X(next_power_of_2_exp_int)(k_end+1,&Nk,&tk);
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

    memset(set->result,0U,data->k_start*sizeof(double _Complex));
    memcpy(&set->result[data->k_start],x,(k_end-data->k_start+1)*sizeof(double _Complex));

    /*eval_sum_clenshaw(k_end, k_end, set->result, set->xc_slow,
      y, set->work, &data->alpha[1], &data->beta[1], &data->gamma[1],
      data->gamma_m1);*/
    eval_sum_clenshaw_fast(k_end, k_end, set->result, set->xc_slow,
      y, &data->_alpha[1], &data->_beta[1], &data->_gamma[1], data->gamma_m1);
  }
  else
  {
    memset(set->temp,0U,data->k_start*sizeof(double _Complex));
    memcpy(&set->temp[data->k_start],x,(k_end-data->k_start+1)*sizeof(double _Complex));

    eval_sum_clenshaw_fast(k_end, Nk-1, set->temp, set->xcvecs[tk-2],
      set->result, &data->_alpha[1], &data->_beta[1], &data->_gamma[1],
      data->gamma_m1);

    fftw_execute_r2r(set->plans_dct2[tk-2],(double*)set->result,
      (double*)set->result);

    set->result[0] *= 0.5;
    for (j = 0; j < Nk; j++)
    {
      set->result[j] *= norm;
    }

    memcpy(y,set->result,(k_end+1)*sizeof(double _Complex));
  }
}

void fpt_trafo(fpt_set set, const int m, const double _Complex *x, double _Complex *y,
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

  double _Complex *work_ptr;
  const double _Complex *x_ptr;

  /* Check, if slow transformation should be used due to small bandwidth. */
  if (k_end < FPT_BREAK_EVEN)
  {
    /* Use NDSFT. */
    fpt_trafo_direct(set, m, x, y, k_end, flags);
    return;
  }

  X(next_power_of_2_exp_int)(k_end,&Nk,&tk);
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
    plan = fftw_plan_many_r2r(1, &length, 2, (double*)set->work, NULL, 2, 1,
      (double*)set->work, NULL, 2, 1, kinds, 0U);
#ifdef _OPENMP
}
#endif
  }

  /* Initialize working arrays. */
  memset(set->result,0U,2*Nk*sizeof(double _Complex));

  /* The first step. */

  /* Set the first 2*data->k_start coefficients to zero. */
  memset(set->work,0U,2*data->k_start*sizeof(double _Complex));

  work_ptr = &set->work[2*data->k_start];
  x_ptr = x;

  for (k = 0; k <= k_end_tilde-data->k_start; k++)
  {
    *work_ptr++ = *x_ptr++;
    *work_ptr++ = K(0.0);
  }

  /* Set the last 2*(set->N-1-k_end_tilde) coefficients to zero. */
  memset(&set->work[2*(k_end_tilde+1)],0U,2*(Nk-1-k_end_tilde)*sizeof(double _Complex));

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
      memcpy(set->vec3,&(set->work[(plength/2)*(4*l+2)]),(plength/2)*sizeof(double _Complex));
      memcpy(set->vec4,&(set->work[(plength/2)*(4*l+3)]),(plength/2)*sizeof(double _Complex));
      memset(&set->vec3[plength/2],0U,(plength/2)*sizeof(double _Complex));
      memset(&set->vec4[plength/2],0U,(plength/2)*sizeof(double _Complex));

      /* Copy coefficients into first half. */
      memcpy(&(set->work[(plength/2)*(4*l+2)]),&(set->work[(plength/2)*(4*l+1)]),(plength/2)*sizeof(double _Complex));
      memset(&(set->work[(plength/2)*(4*l+1)]),0U,(plength/2)*sizeof(double _Complex));
      memset(&(set->work[(plength/2)*(4*l+3)]),0U,(plength/2)*sizeof(double _Complex));

      /* Get matrix U_{n,tau,l} */
      step = &(data->steps[tau][l]);

      /* Check if step is stable. */
      if (step->stable)
      {
        /* Check, if we should do a symmetrizised step. */
        if ((set->flags & FPT_AL_SYMMETRY) && IS_SYMMETRIC(l,m,plength))
        {
          int clength = 1<<(tau);
          double *a11 = step->a;
          double *a12 = a11+clength;
          double *a21 = a12+clength;
          double *a22 = a21+clength;
          /*for (k = 0; k < plength; k++)
          {
            fprintf(stderr,"fpt_trafo: a11 = %le, a12 = %le, a21 = %le, a22 = %le\n",
              a11[k],a12[k],a21[k],step->a22[k]);
          }*/
          /* Multiply third and fourth polynomial with matrix U. */
          //fprintf(stderr,"\nhallo\n");
          fpt_do_step_symmetric(set->vec3, set->vec4, a11,
            a12, a21, a22, step->g, tau, set);
        }
        else
        {
          int clength = 1<<(tau+1);
          double *a11 = step->a;
          double *a12 = a11+clength;
          double *a21 = a12+clength;
          double *a22 = a21+clength;
          /* Multiply third and fourth polynomial with matrix U. */
            fpt_do_step(set->vec3, set->vec4, a11, a12,
              a21, a22, step->g, tau, set);
        }

        if (step->g != 0.0)
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
        memset(&set->vec3[plength/2],0U,(plength_stab-plength/2)*sizeof(double _Complex));
        memset(&set->vec4[plength/2],0U,(plength_stab-plength/2)*sizeof(double _Complex));

        /* Multiply third and fourth polynomial with matrix U. */
        /* Check for symmetry. */
        if (set->flags & FPT_AL_SYMMETRY)
        {
          if (m <= 1)
          {
            int clength_1 = plength_stab;
            int clength_2 = plength_stab;
            double *a11 = step->a;
            double *a12 = a11+clength_1;
            double *a21 = a12+clength_1;
            double *a22 = a21+clength_2;
            fpt_do_step_symmetric(set->vec3, set->vec4, a11, a12,
              a21, a22, step->g, t_stab-1, set);
          }
          else if (m%2 == 0)
          {
            int clength = plength_stab/2;
            double *a11 = step->a;
            double *a12 = a11+clength;
            double *a21 = NULL;
            double *a22 = NULL;
            fpt_do_step_symmetric_u(set->vec3, set->vec4, a11, a12,
              a21, a22,
              set->xcvecs[t_stab-2], step->g, t_stab-1, set);
          }
          else
          {
              int clength = plength_stab/2;
              double *a11 = NULL;
              double *a12 = NULL;
              double *a21 = step->a;
              double *a22 = a21+clength;
              fpt_do_step_symmetric_l(set->vec3, set->vec4,
                a11, a12,
                a21,
                a22, set->xcvecs[t_stab-2], step->g, t_stab-1, set);
          }
        }
        else
        {
          int clength_1 = plength_stab;
          int clength_2 = plength_stab;
          double *a11 = step->a;
          double *a12 = a11+clength_1;
          double *a21 = a12+clength_1;
          double *a22 = a21+clength_2;
          fpt_do_step(set->vec3, set->vec4, a11, a12,
            a21, a22, step->g, t_stab-1, set);
        }

        if (step->g != 0.0)
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

void fpt_transposed_direct(fpt_set set, const int m, double _Complex *x,
  double _Complex *y, const int k_end, const unsigned int flags)
{
  int j;
  fpt_data *data = &(set->dpt[m]);
  int Nk;
  int tk;
  double norm;

  X(next_power_of_2_exp_int)(k_end+1,&Nk,&tk);
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

    memcpy(x,&set->result[data->k_start],(k_end-data->k_start+1)*
      sizeof(double _Complex));
  }
  else
  {
    memcpy(set->result,y,(k_end+1)*sizeof(double _Complex));
    memset(&set->result[k_end+1],0U,(Nk-k_end-1)*sizeof(double _Complex));

    for (j = 0; j < Nk; j++)
    {
      set->result[j] *= norm;
    }

    fftw_execute_r2r(set->plans_dct3[tk-2],(double*)set->result,
      (double*)set->result);

    if (set->N > 1024)
      eval_sum_clenshaw_transposed_ld(k_end, Nk-1, set->temp, set->xcvecs[tk-2],
        set->result, set->work, &data->_alpha[1], &data->_beta[1], &data->_gamma[1],
        data->gamma_m1);
    else
      eval_sum_clenshaw_transposed(k_end, Nk-1, set->temp, set->xcvecs[tk-2],
        set->result, set->work, &data->_alpha[1], &data->_beta[1], &data->_gamma[1],
        data->gamma_m1);

    memcpy(x,&set->temp[data->k_start],(k_end-data->k_start+1)*sizeof(double _Complex));
  }
}

void fpt_transposed(fpt_set set, const int m, double _Complex *x,
  double _Complex *y, const int k_end, const unsigned int flags)
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

  X(next_power_of_2_exp_int)(k_end,&Nk,&tk);
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
    plan = fftw_plan_many_r2r(1, &length, 2, (double*)set->work, NULL, 2, 1,
      (double*)set->work, NULL, 2, 1, kinds, 0U);
#ifdef _OPENMP
}
#endif
    fftw_execute_r2r(plan,(double*)y,(double*)set->result);
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
    memcpy(set->result,y,(k_end+1)*sizeof(double _Complex));
  }

  /* Initialize working arrays. */
  memset(set->work,0U,2*Nk*sizeof(double _Complex));

  /* The last step is now the first step. */
  for (k = 0; k <= k_end; k++)
  {
    set->work[k] = data->gamma_m1*set->result[k];
  }
  //memset(&set->work[k_end+1],0U,(Nk+1-k_end)*sizeof(double _Complex));

  set->work[Nk] = data->gamma_m1*(data->beta_0*set->result[0] +
    data->alpha_0*set->result[1]);
  for (k = 1; k < k_end; k++)
  {
    set->work[Nk+k] = data->gamma_m1*(data->beta_0*set->result[k] +
      data->alpha_0*0.5*(set->result[k-1]+set->result[k+1]));
  }
  if (k_end<Nk)
  {
    memset(&set->work[k_end],0U,(Nk-k_end)*sizeof(double _Complex));
  }

  /** Save copy of inpute data for stabilization steps. */
  memcpy(set->result,set->work,2*Nk*sizeof(double _Complex));

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
      memcpy(set->vec3,&(set->work[(plength/2)*(4*l+0)]),plength*sizeof(double _Complex));
      memcpy(set->vec4,&(set->work[(plength/2)*(4*l+2)]),plength*sizeof(double _Complex));

      memcpy(&set->work[(plength/2)*(4*l+1)],&(set->work[(plength/2)*(4*l+2)]),
        (plength/2)*sizeof(double _Complex));

      /* Get matrix U_{n,tau,l} */
      step = &(data->steps[tau][l]);

      /* Check if step is stable. */
      if (step->stable)
      {
        if ((set->flags & FPT_AL_SYMMETRY) && IS_SYMMETRIC(l,m,plength))
        {
          /* Multiply third and fourth polynomial with matrix U. */
          int clength = 1<<(tau);
          double *a11 = step->a;
          double *a12 = a11+clength;
          double *a21 = a12+clength;
          double *a22 = a21+clength;
          fpt_do_step_t_symmetric(set->vec3, set->vec4, a11, a12,
            a21, a22, step->g, tau, set);
        }
        else
        {
          /* Multiply third and fourth polynomial with matrix U. */
          int clength = 1<<(tau+1);
          double *a11 = step->a;
          double *a12 = a11+clength;
          double *a21 = a12+clength;
          double *a22 = a21+clength;
          fpt_do_step_t(set->vec3, set->vec4, a11, a12,
            a21, a22, step->g, tau, set);
        }
        memcpy(&(set->vec3[plength/2]), set->vec4,(plength/2)*sizeof(double _Complex));

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

        memcpy(set->vec3,set->result,plength_stab*sizeof(double _Complex));
        memcpy(set->vec4,&(set->result[Nk]),plength_stab*sizeof(double _Complex));

        /* Multiply third and fourth polynomial with matrix U. */
        if (set->flags & FPT_AL_SYMMETRY)
        {
          if (m <= 1)
          {
            int clength_1 = plength_stab;
            int clength_2 = plength_stab;
            double *a11 = step->a;
            double *a12 = a11+clength_1;
            double *a21 = a12+clength_1;
            double *a22 = a21+clength_2;
            fpt_do_step_t_symmetric(set->vec3, set->vec4, a11, a12,
              a21, a22, step->g, t_stab-1, set);
          }
          else if (m%2 == 0)
          {
            int clength = plength_stab/2;
            double *a11 = step->a;
            double *a12 = a11+clength;
            fpt_do_step_t_symmetric_u(set->vec3, set->vec4, a11, a12,
              set->xcvecs[t_stab-2], step->g, t_stab-1, set);
          }
          else
          {
            int clength = plength_stab/2;
            double *a21 = step->a;
            double *a22 = a21+clength;
            fpt_do_step_t_symmetric_l(set->vec3, set->vec4,
              a21, a22, set->xcvecs[t_stab-2], step->g, t_stab-1, set);
          }
        }
        else
        {
          int clength_1 = plength_stab;
          int clength_2 = plength_stab;
          double *a11 = step->a;
          double *a12 = a11+clength_1;
          double *a21 = a12+clength_1;
          double *a22 = a21+clength_2;
          fpt_do_step_t(set->vec3, set->vec4, a11, a12,
            a21, a22, step->g, t_stab-1, set);
        }

        memcpy(&(set->vec3[plength/2]),set->vec4,(plength/2)*sizeof(double _Complex));

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
  int k_start_tilde;
  int N_tilde;
  int firstl, lastl;
  int plength;
  const int M = set->M;

  /* TODO Clean up DPT transform data structures. */
  if (!(set->flags & FPT_NO_INIT_FPT_DATA))
  {
    for (m = 0; m < M; m++)
    {
      /* Check if precomputed. */
      fpt_data *data = &set->dpt[m];
      if (data->steps != (fpt_step**)NULL)
      {
        if (!(set->flags & FPT_NO_FAST_ALGORITHM))
        {
          nfft_free(data->alphaN);
          data->alphaN = NULL;
          data->betaN = NULL;
          data->gammaN = NULL;
        }

        /* Free precomputed data. */
        k_start_tilde = K_START_TILDE(data->k_start,X(next_power_of_2)(data->k_start)
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
            if (data->steps[tau][l].a != NULL)
            {
              nfft_free(data->steps[tau][l].a);
              data->steps[tau][l].a = NULL;
            }
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
        if (!(set->flags & FPT_PERSISTENT_DATA))
        {
          if (data->_alpha != NULL)
            nfft_free(data->_alpha);
        }
        data->_alpha = NULL;
        data->_beta = NULL;
        data->_gamma = NULL;
      }
    }

    /* Delete array of DPT transform data. */
    nfft_free(set->dpt);
    set->dpt = NULL;
  }

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
  set->work = NULL;
  set->result = NULL;

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

  /* Check if fast transform is activated. */
  if (!(set->flags & FPT_NO_FAST_ALGORITHM))
  {
    /* Free auxilliary arrays. */
    nfft_free(set->vec3);
    nfft_free(set->vec4);
    nfft_free(set->z);
    set->vec3 = NULL;
    set->vec4 = NULL;
    set->z = NULL;
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
