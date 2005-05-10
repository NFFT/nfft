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

#define FIRST_L (n/plength)/2;
#define LAST_L (plength*(int)ceil(((double)M)/plength)-1)/plength; 

/** \defgroup nfsft_internal NFSFT: Internal API */

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
  /** The number of nodes. */
  int D;                      
  /** The bandwidth */
  int M;                               
  /** The threshold */
  double threshold;
  /** The angles phi of the nodes */
  double *angles;
  /** The fourier coefficients. */
  complex **f_hat;
  /** The function values. */
  complex *f;
};

/** 
 * Structure for wisdom for a specific bandwidth. 
 *
 * \ingroup nfsft_internal
 */
struct nfsft_transform_wisdom
{
  /** The bandwidth */
  int N;
  /** The logarithm of the bandwidth */
  int t;
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
  /* Structure for matrices U */
  struct U_type ***U;    
  /** For FLFT */
  complex *work,*old,*vec1,*vec2,*vec3,*vec4, *a2, *b2;
  complex *ergeb;
  /** NFFT plan */
  nfft_plan plan;  
};

/** 
 * Toplevel wisdom structure. 
 *
 * \ingroup nfsft_internal
 */
struct nfsft_wisdom
{
  /** Indicates wether the structure has been initialized. */ 
  bool initialized;
  /** Array of wisdoms for bandwidth values in powers of 2 */
  struct nfsft_transform_wisdom **transform_wisdoms;  
  /** Precomputed recursion coefficients for associated Legendre-functions */
  double alpha[(BW_MAX+1)*(BW_MAX+1)];
  /** Precomputed recursion coefficients for associated Legendre-functions */
  double beta[(BW_MAX+1)*(BW_MAX+1)];
  /** Precomputed recursion coefficients for associated Legendre-functions */
  double gamma[(BW_MAX+1)*(BW_MAX+1)];
  /** Precomputed recursion coefficients for associated Legendre-functions */
  double gamma_m1[BW_MAX+1];
};
#endif
