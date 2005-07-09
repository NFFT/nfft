/**
 * Header file with internal API of the NFSFT library
 */
#ifndef API_H
#define API_H

#include "nfsft.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#else
#  error Need config.h
#endif

#ifdef STDC_HEADERS
#  include <stdio.h>
#  include <math.h>
#  include <string.h>
#else
#  error Need ANSI-C headers.
#endif

#ifdef HAVE_STDBOOL_H
#  include <stdbool.h>
#else
#  warning ISO C99 bool type not available. Defining own bool type.
typedef enum {false = 0,true = 1} bool;
#endif

#ifdef HAVE_FFTW3_H
# include <fftw3.h>
#else
# error Need fftw3.h
#endif

#include "../nfft/nfft.h"

#define ROW(k) (k*(BW_MAX+1))
#define ROWK(k) (k*(BW_MAX+1)+k)

#define FIRST_L (n/plength)/2
#define LAST_L (plength*(int)ceil(((double)M)/plength)-1)/plength

/** \defgroup nfsft_internal NFSFT: Internal API and functions */

/** 
 * Datatype for a set of real 2x2 matrices used in FLT. 
 *
 * \ingroup nfsft_internal
 */
struct U_type
{
  /** 
   * Indicates if the values contained represent a fast or a slow stabilized 
   * step.
   */
  bool stable;
  /** The components */
  double *m1,*m2,*m3,*m4;
};

/** 
 * Structure for a transform plan. 
 *
 * \ingroup nfsft_internal
 */
struct nfsft_plan_s
{
  /** The flags */
  int flags;
  /** The number of nodes. */
  int D;                      
  /** The bandwidth */
  int M;
  /** Next greater power of two with respect to M */
  int N;
  int t;
  /** The angles phi of the nodes */
  double *angles;
  /** The fourier coefficients. */
  complex **f_hat;
  /** The function values. */
  complex *f;
  /** NFFT plan */
  nfft_plan plan_nfft;  
};

/** 
* Structure for an inverse transform plan. 
*
* \ingroup nfsft_internal
*/
struct infsft_plan_s
{
  nfsft_plan direct_plan;
  int infsft_flags;
  complex *given_f;
  complex *r_iter;
  complex *v_iter;
  complex **f_hat_iter;
  complex **f_hat_iter_2nd;
  complex **p_hat_iter;
  complex **z_hat_iter;
  double *w;
  double **w_hat;
  double dot_r_iter;
  double dot_r_iter_old;
  double dot_v_iter;
  double alpha_iter;
  double alpha_iter_2nd;
  double beta_iter;
  double gamma_iter;
  double gamma_iter_old;
  double dot_alpha_iter;
  double dot_z_hat_iter;
  double dot_z_hat_iter_old;
  double dot_p_hat_iter;
  
  /** The number of nodes. */
//  int D;                      
  /** The bandwidth */
//  int M;
  /** Next greater power of two with respect to M */
//  int N;
  /** The angles phi of the nodes */
//  double *angles;
  /** The fourier coefficients. */
//  complex **f_hat;
  /** The function values. */
//  complex *f;
  /** NFFT plan */
//  nfft_plan plan_nfft;  
};


/** 
 * Toplevel wisdom structure. 
 *
 * \ingroup nfsft_internal
 */
struct nfsft_wisdom
{
  /** Indicates wether the structure has been initialized */ 
  bool initialized;
  nfsft_precompute_flags flags;
  /** The logarithm of the bandwidth */
  int t;
  /** Maximum bandwidth */
  int N;
  /** Precomputed recursion coefficients for associated Legendre-functions */
  double *alpha;
  /** Precomputed recursion coefficients for associated Legendre-functions */
  double *beta;
  /** Precomputed recursion coefficients for associated Legendre-functions */
  double *gamma;
  /** The threshold */
  double threshold;
  /* Structure for matrices U */
  struct U_type ****U;    
  /** For FLFT */
  complex *work,*old,*vec1,*vec2,*vec3,*vec4, *a2, *b2;
  /** Transform plans for the fftw library */
  fftw_plan *plans_dct3;
  /** Transform plans for the fftw library */
  fftw_plan *plans_dct2;  
  /** Transform kinds for fftw library */
  fftw_r2r_kind *kinds;
  /** Transform kinds for fftw library */
  fftw_r2r_kind *kindsr;
  /** Transform lengths for fftw library */
  int *lengths;
  complex *ergeb;
  
  
  double *flft_alpha;
  double *flft_beta;
  double *flft_gamma;

  complex *z;
};
#endif
