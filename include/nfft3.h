/*! \file nfft3.h
 *  \brief Header file for the nfft3 library.
 */
#ifndef NFFT3_H
#define NFFT3_H

/** Include header for C99 complex datatype. */
#include <complex.h>
/** Include header for FFTW3 library. */
#include <fftw3.h>

/** Macros for public members inherited by all plan structures. */
#define MACRO_MV_PLAN(float_type)                                           \
int N_total;                          /**< total number of Fourier coeffs.*/\
int M_total;                          /**< total number of samples        */\
float_type *f_hat;                    /**< Fourier coefficients           */\
float_type *f;                        /**< samples                        */\

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nfft NFFT
 * Direct and fast computation of the
 * discrete Fourier transform at nonequispaced knots
 * @{
 */

/**
 * Constant symbols for precomputation and memory usage
 */
#define PRE_PHI_HUT      (1U<< 0)
#define FG_PSI           (1U<< 1)
#define PRE_LIN_PSI      (1U<< 2)
#define PRE_FG_PSI       (1U<< 3)
#define PRE_PSI          (1U<< 4)
#define PRE_FULL_PSI     (1U<< 5)
#define MALLOC_X         (1U<< 6)
#define MALLOC_F_HAT     (1U<< 7)
#define MALLOC_F         (1U<< 8)
#define FFT_OUT_OF_PLACE (1U<< 9)
#define FFTW_INIT        (1U<< 10)

#define MALLOC_V         (1U<< 11)

#define SNDFT            (1U<< 12)

#define PRE_ONE_PSI (PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)

typedef struct nfft_plan_
{
  /** api */
  MACRO_MV_PLAN(complex);

  int d;                                /**< dimension, rank                 */
  int *N;                               /**< multi bandwidth                 */
  double *sigma;                        /**< oversampling-factor             */
  int *n;                               /**< fftw-length = sigma*N           */
  int n_total;                          /**< total size of fftw              */
  int m;                                /**< cut-off, window function        */
  double *b;                            /**< shape parameters                */
  int K;                                /**< number of precomp. uniform psi  */

  unsigned nfft_flags;                  /**< flags for precomputation, malloc*/
  unsigned fftw_flags;                  /**< flags for the fftw              */

  double *x;                            /**< nodes (in time/spatial domain)  */

  double MEASURE_TIME_t[3];             /**< measured time for each step     */

  /** internal*/
  fftw_plan  my_fftw_plan1;             /**< fftw_plan forward               */
  fftw_plan  my_fftw_plan2;             /**< fftw_plan backward              */

  double **c_phi_inv;                   /**< precomputed data, matrix D      */
  double *psi;                          /**< precomputed data, matrix B      */
  int *psi_index_g;                     /**< only for PRE_FULL_PSI           */
  int *psi_index_f;                     /**< only for PRE_FULL_PSI           */

  complex *g;
  complex *g_hat;
  complex *g1;                          /**< input of fftw                   */
  complex *g2;                          /**< output of fftw                  */

  double *spline_coeffs;            /**< input for de Boor algorithm, if
               B_SPLINE or SINC_2m is defined  */
} nfft_plan;


/**
 * Executes a NDFT, see equation (1.1) in [Guide], computes
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void ndft_trafo(nfft_plan *ths);

/**
 * Executes a NDFT, see equation (1.2) in [Guide], computes
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void ndft_adjoint(nfft_plan *ths);

/**
 * Executes a NFFT, see equation (1.1) in [Guide], computes fast and
 * approximate
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_trafo(nfft_plan *ths);

/**
 * Executes an adjoint  NFFT, see equation (1.2) in [Guide], computes fast and
 * approximate
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_adjoint(nfft_plan *ths);

/**
 * Initialisation of a transform plan, wrapper d=1.
 *
 * \arg ths The pointer to a nfft plan
 * \arg N1 bandwidth
 * \arg M_total The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_1d(nfft_plan *ths, int N1, int M);

/**
 * Initialisation of a transform plan, wrapper d=2.
 *
 * \arg ths The pointer to a nfft plan
 * \arg N1 bandwidth
 * \arg N2 bandwidth
 * \arg M_total The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M);

/**
 * Initialisation of a transform plan, wrapper d=3.
 *
 * \arg ths The pointer to a nfft plan
 * \arg N1 bandwidth
 * \arg N2 bandwidth
 * \arg N3 bandwidth
 * \arg M_total The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_3d(nfft_plan *ths, int N1, int N2, int N3, int M);

/**
 * Initialisation of a transform plan, simple.
 *
 * \arg ths The pointer to a nfft plan
 * \arg d The dimension
 * \arg N The multi bandwidth
 * \arg M_total The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init(nfft_plan *ths, int d, int *N, int M);

/**
 * Initialisation of a transform plan, advanced.
 *
 * \arg ths The pointer to a nfft plan
 * \arg d The dimension
 * \arg N The multi bandwidth
 * \arg M_total The number of nodes
 * \arg nfft_flags_on NFFT flags to switch on
 * \arg nfft_flags_off NFFT flags to switch off
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_advanced(nfft_plan *ths, int d, int *N, int M,
      unsigned nfft_flags_on, unsigned nfft_flags_off);

/**
 * Initialisation of a transform plan, guru.
 *
 * \arg ths The pointer to a nfft plan
 * \arg d The dimension
 * \arg N The multi bandwidth
 * \arg M_total The number of nodes
 * \arg n The oversampled multi bandwidth
 * \arg m The spatial cut-off
 * \arg nfft_flags_on NFFT flags to switch on
 * \arg nfft_flags_off NFFT flags to switch off
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_guru(nfft_plan *ths, int d, int *N, int M, int *n,
              int m, unsigned nfft_flags, unsigned fftw_flags);

void nfft_precompute_full_psi(nfft_plan *ths);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis
 *
 * wrapper for precompute*_psi
 *
 * if PRE_*_PSI is set the application program has to call this routine
 * (after) setting the nodes x
 */
void nfft_precompute_one_psi(nfft_plan *ths);

/**
 * Destroys a transform plan.
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_finalize(nfft_plan *ths);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nfct NFCT
 * direct and fast computation of the
 * discrete cosine transform at nonequispaced knots in time/spatial domain
 * @{
 */

/** Structure for a transform plan */
typedef struct nfct_plan_
{
  /** api */
  MACRO_MV_PLAN(double);

  int d;                                /**< dimension, rank                  */
  int *N;                               /**< cut-off-frequencies (kernel)     */
  int *n;                               /**< length of dct-i                  */
  double *sigma;                        /**< oversampling-factor              */
  int m;                                /**< cut-off parameter in time-domain */

  double nfct_full_psi_eps;
  double *b;                            /**< shape parameters                 */

  unsigned nfct_flags;                  /**< flags for precomputation, malloc */
  unsigned fftw_flags;                  /**< flags for the fftw               */

  double *x;                            /**< nodes (in time/spatial domain)   */

  double MEASURE_TIME_t[3];             /**< measured time for each step     */

  /** internal */
  fftw_plan  my_fftw_r2r_plan;          /**< fftw_plan                        */
  fftw_r2r_kind *r2r_kind;              /**< r2r transform type (dct-i)       */

  double **c_phi_inv;                   /**< precomputed data, matrix D       */
  double *psi;                          /**< precomputed data, matrix B       */
  int size_psi;                         /**< only for thin B                  */
  int *psi_index_g;                     /**< only for thin B                  */
  int *psi_index_f;                     /**< only for thin B                  */

  double *g;
  double *g_hat;
  double *g1;                           /**< input of fftw                    */
  double *g2;                           /**< output of fftw                   */

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined   */
} nfct_plan;


/**
 * Creates a 1-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth \f$N\f$
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init_1d( nfct_plan *ths_plan, int N0, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init_2d( nfct_plan *ths_plan, int N0, int N1, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg N2 The bandwidth of dimension 3
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init_3d( nfct_plan *ths_plan, int N0, int N1, int N2, int M_total);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init( nfct_plan *ths_plan, int d, int *N, int M_total);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 * \arg n The oversampled bandwidths
 * \arg m The cut-off parameter
 * \arg nfct_flags The flags known to nfct
 * \arg fftw_flags The flags known to fftw
 *
 * \author Steffen Klatt
 */
void nfct_init_guru( nfct_plan *ths_plan, int d, int *N, int M_total, int *n,
                         int m, unsigned nfct_flags, unsigned fftw_flags);

/**
 * precomputes the values psi
 * if the PRE_PSI is set the application program has to call this routine
 * after setting the nodes this_plan->x
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_precompute_psi( nfct_plan *ths_plan);


/**
 * executes a NFCT (approximate,fast), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^C(x_j) = sum_{k \in I_0^{N,d}} \hat{f}_k^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_trafo( nfct_plan *ths_plan);

/**
 * executes a NDCT (exact,slow), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^C(x_j) = sum_{k \in I_0^{N,d}} \hat{f}_k^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndct_trafo( nfct_plan *ths_plan);

/**
 * executes a transposed NFCT (approximate,fast), computes for \f$k \in I_0^{N,d}\f$
 * \f$h^C(k) = sum_{j \in I_0^{(M\_total,1)}} f_j^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_adjoint( nfct_plan *ths_plan);

/**
 * executes a direct transposed NDCT (exact,slow), computes for \f$k \in I_0^{N,d}\f$
 * \f$h^C(k) = sum_{j \in I_0^{(M\_total,1)}} f_j^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndct_adjoint( nfct_plan *ths_plan);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_finalize( nfct_plan *ths_plan);

/**
 * do some adjustments (N,n) then compute PHI_HUT
 *
 * \arg ths_plan the plan for the transform
 * \arg k        index of c_phi
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfct_phi_hut( nfct_plan *ths_plan, int k, int d);

/**
 * do some adjustments (N,n) then compute PHI
 *
 * \arg ths_plan the plan for the transform
 * \arg x        node \f$x\f$
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfct_phi ( nfct_plan *ths_plan, double x, int d);

/**
 * returns 2(n-1),  fftw related issue
 *
 * \arg n       i.e. length of dct-1
 *
 * \author Steffen Klatt
 */
int nfct_fftw_2N( int n);

/**
 * returns 0.5n+1,  fftw related issue
 *
 * \arg n       i.e. length of dct-1
 *
 * \author Steffen Klatt
 */
int nfct_fftw_2N_rev( int n);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nfst NFST
 * direct and fast computation of the
 * discrete sine transform at nonequispaced knots in time/spatial domain
 * @{
 */

/** Structure for a transform plan */
typedef struct nfst_plan_
{
  /** api */
  MACRO_MV_PLAN(double);

  int d;                                /**< dimension, rank                  */
  int *N;                               /**< bandwidth                        */
  int *n;                               /**< length of dst-1                  */
  double *sigma;                        /**< oversampling-factor              */
  int m;                                /**< cut-off parameter in time-domain */

  double nfst_full_psi_eps;
  double *b;                            /**< shape parameters                 */

  unsigned nfst_flags;                  /**< flags for precomputation, malloc */
  unsigned fftw_flags;                  /**< flags for the fftw               */

  double *x;                            /**< nodes (in time/spatial domain)   */

  double MEASURE_TIME_t[3];             /**< measured time for each step     */

  /** internal */
  fftw_plan  my_fftw_r2r_plan;         /**< fftw_plan forward                */
  fftw_r2r_kind *r2r_kind;              /**< r2r transform type (dct-i)       */

  double **c_phi_inv;                   /**< precomputed data, matrix D       */
  double *psi;                          /**< precomputed data, matrix B       */
  int size_psi;                         /**< only for thin B                  */
  int *psi_index_g;                     /**< only for thin B                  */
  int *psi_index_f;                     /**< only for thin B                  */


  double *g;
  double *g_hat;
  double *g1;                           /**< input of fftw                    */
  double *g2;                           /**< output of fftw                   */

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined   */
} nfst_plan;


/**
 * Creates a 1-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth \f$N\f$
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init_1d( nfst_plan *ths_plan, int N0, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init_2d( nfst_plan *ths_plan, int N0, int N1, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg N2 The bandwidth of dimension 3
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init_3d( nfst_plan *ths_plan, int N0, int N1, int N2, int M_total);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init( nfst_plan *ths_plan, int d, int *N, int M_total);

/**
 * Creates a d-dimensional transform plan with pcific m.
 * (just for convenience)
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 * \arg m cut-off parameter
 *
 * \author Steffen Klatt
 */
void nfst_init_m( nfst_plan *ths_plan, int d, int *N, int M_total, int m);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 * \arg n The oversampled bandwidths
 * \arg m The cut-off parameter
 * \arg nfst_flags The flags known to nfst
 * \arg fftw_flags The flags known to fftw
 *
 * \author Steffen Klatt
 */
void nfst_init_guru( nfst_plan *ths_plan, int d, int *N, int M_total, int *n,
                     int m, unsigned nfst_flags, unsigned fftw_flags);

/**
 * precomputes the values psi
 * if the PRE_PSI is set the application program has to call this routine
 * after setting the nodes this_plan->x
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_precompute_psi( nfst_plan *ths_plan);

/**
 * executes a NFST (approximate,fast), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^S(x_j) = sum_{k \in I_1^{N,d}} \hat{f}_k^S * sin(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_trafo( nfst_plan *ths_plan);

/**
 * executes a NDST (exact,slow), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^S(x_j) = sum_{k \in I_1^{N,d}} \hat{f}_k^S * sin(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndst_trafo( nfst_plan *ths_plan);



/**
 * executes a transposed NFST (approximate,fast), computes for \f$k \in I_1^{N,d}\f$
 * \f$h^S(k) = sum_{j \in I_0^{M\_total,1}} f_j^S * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_adjoint( nfst_plan *ths_plan);

/**
 * executes a direct transposed NDST (exact,slow), computes for \f$k \in I_1^{N,d}\f$
 * \f$h^S(k) = sum_{j \in I_0^{M\_total,1}} f_j^S * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndst_adjoint( nfst_plan *ths_plan);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_finalize( nfst_plan *ths_plan);

/**
 *
 *
 * \arg ths_plan The plan for the transform
 *
 */
void nfst_full_psi( nfst_plan *ths_plan, double eps);

/**
 * do some adjustments (N,n) then compute PHI_HUT
 *
 * \arg ths_plan the plan for the transform
 * \arg k        index of c_phi
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfst_phi_hut( nfst_plan *ths_plan, int k, int d);

/**
 * do some adjustments (N,n) then compute PHI
 *
 * \arg ths_plan the plan for the transform
 * \arg x        node \f$x\f$
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfst_phi ( nfst_plan *ths_plan, double x, int d);

/**
 * returns 2(n+1),  fftw related issue
 *
 * \arg n       i.e. length of dst-1
 *
 * \author Steffen Klatt
 */
int nfst_fftw_2N( int n);

/**
 * returns 0.5n-1,  fftw related issue
 *
 * \arg n       i.e. length of dct-1
 *
 * \author Steffen Klatt
 */
int nfst_fftw_2N_rev( int n);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nnfft NNFFT
 * Direct and fast computation of the
 * discrete Fourier transform at nonequispaced knots in time and Fourier domain
 * @{
 */

/**
 * Structure for a transform plan
 */
typedef struct
{
  /** api */
  MACRO_MV_PLAN(complex);

  int d;                                /**< dimension, rank                 */
  double *sigma;                        /**< oversampling-factor            */
  double *a;                            /**< 1 + 2*m/N1                     */
  int *N;                               /**< cut-off-frequencies             */
  int *N1;                              /**< sigma*N                         */
  int *aN1;                             /**< sigma*a*N                       */
  int m;                                /**< cut-off parameter in time-domain*/
  double *b;                            /**< shape parameters                */
  int K;                                /**< number of precomp. uniform psi  */

  int aN1_total;                        /**< aN1_total=aN1[0]* ... *aN1[d-1] */

  nfft_plan *direct_plan;               /**< plan for the nfft               */
  unsigned nnfft_flags;                 /**< flags for precomputation, malloc*/
  int *n;                               /**<  n=N1, just for the window function */

  double *x;                            /**< nodes (in time/spatial domain)  */
  double *v;                            /**< nodes (in fourier domain)       */

  double *c_phi_inv;                    /**< precomputed data, matrix D      */
  double *psi;                          /**< precomputed data, matrix B      */
  int size_psi;                         /**< only for thin B                 */
  int *psi_index_g;                     /**< only for thin B                 */
  int *psi_index_f;                     /**< only for thin B                 */
  complex *F;

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined  */
} nnfft_plan;


/**
 * Creates a transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d The dimension
 * \arg N_total The number of nodes \f$v\f$
 * \arg M_total The number of nodes \f$x\f$
 * \arg N The bandwidth \f$N\f$
 *
 * \author Tobias Knopp
 */
void nnfft_init(nnfft_plan *ths_plan, int d, int N_total, int M_total, int *N);

/**
 * Creates a transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d The dimension
 * \arg N_total The number of nodes \f$v\f$
 * \arg M_total The number of nodes \f$x\f$
 * \arg N The bandwidth \f$N\f$
 * \arg N1 The oversampled bandwidth \f$N\f$
 * \arg m The cut-off parameter
 * \arg nnfft_flags The flags
 *
 * \author Tobias Knopp
 */
void nnfft_init_guru(nnfft_plan *ths_plan, int d, int N_total, int M_total,
                     int *N, int *N1, int m, unsigned nnfft_flags);

/**
 * Executes a direct NNDFT, i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k = 0}^{N_{total}-1} \hat{f}(v_k) {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nndft_trafo(nnfft_plan *ths_plan);

/**
 * Executes a direct adjoint NNDFT, i.e. computes for \f$k=0,...,N_{total}-1\f$
 * \f[
 *   \hat{f}(v_k) = \sum_{j = 0}^{M_{total}-1} f(x_j) {\rm e}^{2 \pi \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nndft_adjoint(nnfft_plan *ths_plan);

/**
 * Executes a NNFFT, i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k = 0}^{N_{total}-1} \hat{f}(v_k) {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nnfft_trafo(nnfft_plan *ths_plan);

/**
 * Executes a adjoint NNFFT, i.e. computes for \f$k=0,...,N_{total}-1\f$
 * \f[
 *   \hat{f}(v_k) = \sum_{j = 0}^{M_{tota}l-1} f(x_j) {\rm e}^{2 \pi \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nnfft_adjoint(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes equally spaced values of the window function psi
 *
 * if PRE_LIN_PSI is set the application program has to call this routine
 */
void nnfft_precompute_lin_psi(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes the values of the window function psi in a tensor product form
 *
 * if PRE_PSI is set the application program has to call this routine after
 * setting the nodes x
 */
void nnfft_precompute_psi(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes the values of the window function psi and their indices in
 * non tensor product form
 *
 * if PRE_FULL_PSI is set the application program has to call this routine
 * after setting the nodes x
 */
void nnfft_precompute_full_psi(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes the values of the fourier transformed window function, i.e. phi_hut
 *
 * if PRE_PHI_HUT is set the application program has to call this routine
 * after setting the nodes v
 */
void nnfft_precompute_phi_hut(nnfft_plan *ths_plan);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nnfft_finalize(nnfft_plan *ths_plan);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nsfft NSFFT
 * Direct and fast computation of the
 * hyperbolic discrete Fourier transform at nonequispaced knots
 * @{
 */
typedef struct nsfft_plan_
{
  MACRO_MV_PLAN(complex);

  int d;                                /**< dimension, rank; d=2,3          */
  int J;                                /**< problem size, i.e.,
                                             d=2: N_total=(J+4) 2^(J+1)
                                             d=3: N_total=2^J 6(2^((J+1)/2+1)
                                                          -1)+2^(3(J/2+1))   */
  int sigma;                            /**< oversampling-factor             */

  unsigned flags;                       /**< flags for precomputation, malloc*/

  int *index_sparse_to_full;            /**< index conversation,
                                             overflow for d=3, J=9!          */

  int r_act_nfft_plan;                  /**< index of current nfft block     */
  nfft_plan *act_nfft_plan;             /**< current nfft block              */
  nfft_plan *center_nfft_plan;          /**< central nfft block              */

  fftw_plan* set_fftw_plan1;            /**< fftw plan for the nfft blocks   */
  fftw_plan* set_fftw_plan2;            /**< fftw plan for the nfft blocks   */

  nfft_plan *set_nfft_plan_1d;          /**< nfft plans for short nffts      */
  nfft_plan *set_nfft_plan_2d;          /**< nfft plans for short nffts      */

  double *x_transposed;                 /**< coordinate exchanged nodes, d=2 */
  double *x_102,*x_201,*x_120,*x_021;   /**< coordinate exchanged nodes, d=3 */

} nsfft_plan;

/**
 * Executes a NSDFT, computes
 * for j=0,...,M-1
 *  f[j] = sum_{k in H_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsdft_trafo(nsfft_plan *ths);

/**
 * Executes an ajoint NSDFT, computes
 * for k in H_N^d
 *  f_hat[k] = sum_{j=0,...,M-1} f[j] * exp(+2 (pi) k x[j])
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Stefan Kunis
 */
void nsdft_adjoint(nsfft_plan *ths);

/**
 * Executes a NSDFT, computes fast and approximate
 * for j=0,...,M-1
 *  f[j] = sum_{k in H_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_trafo(nsfft_plan *ths);

/**
 * Executes a NSDFT, computes fast and approximate
 * for k in H_N^d
 *  f_hat[k] = sum_{j=0,...,M-1} f[j] * exp(+2 (pi) k x[j])
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Stefan Kunis
 */
void nsfft_adjoint(nsfft_plan *ths);

/**
 * Copy coefficients from nsfft plan to a nfft plan.
 *
 * \arg ths Pointers to a nsfft plan and to a nfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_cp(nsfft_plan *ths, nfft_plan *ths_nfft);

/**
 * Initialisation of pseudo random nodes and coefficients.
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_init_random_nodes_coeffs(nsfft_plan *ths);

/**
 * Initialisation of a transform plan.
 *
 * \arg ths The pointer to a nsfft plan
 * \arg d The dimension
 * \arg J The problem size
 * \arg M_total The number of nodes
 * \arg m nfft cut-off parameter
 * \arg flags
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_init(nsfft_plan *ths, int d, int J, int M, int m, unsigned flags);

/**
 * Destroys a transform plan.
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_finalize(nsfft_plan *ths);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup mri_inh MRI_INH
 * @{
 */

/**
 * The structure for the transform plan.
 */
typedef struct
{
  /** api */
  MACRO_MV_PLAN(complex);

  nfft_plan plan;

  int N3;
  double sigma3;
  double *t;
  double *w;
} mri_inh_2d1d_plan;

/**
 * The structure for the transform plan.
 */
typedef struct
{
  /** api */
  MACRO_MV_PLAN(complex);

  nfft_plan plan;

  int N3;
  double sigma3;
  double *t;
  double *w;
} mri_inh_3d_plan;


/**
 * Executes a mri transformation considering the field inhomogeneity with the 2d1d method,
 * i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k \in I_N^2} \hat{f}(k) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_trafo(mri_inh_2d1d_plan *ths);

/**
 * Executes an adjoint mri transformation considering the field inhomogeneity with the 2d1d method,
 * i.e. computes for \f$k \in I_N^2\f$
 * \f[
 *   \hat{f}(k) = \sum_{j=0}^{M_{total}-1} f(x_j) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_adjoint(mri_inh_2d1d_plan *ths);

/**
 * Creates a transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N The bandwidth \f$N\f$
 * \arg M_total The number of nodes \f$x\f$
 * \arg n The oversampled bandwidth \f$N\f$
 * \arg m The cut-off parameter
 * \arg sigma The oversampling factor
 * \arg nnfft_flags The flags
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_init_guru(mri_inh_2d1d_plan *ths, int *N, int M, int *n,
                    int m, double sigma, unsigned nfft_flags, unsigned fftw_flags);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_finalize(mri_inh_2d1d_plan *ths);

/**
 * Executes a mri transformation considering the field inhomogeneity with the 3d method,
 * i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k \in I_N^2} \hat{f}(k) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_3d_trafo(mri_inh_3d_plan *ths);

/**
 * Executes an adjoint mri transformation considering the field inhomogeneity with the 3d method,
 * i.e. computes for \f$k \in I_N^2\f$
 * \f[
 *   \hat{f}(k) = \sum_{j=0}^{M_{total}-1} f(x_j) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_3d_adjoint(mri_inh_3d_plan *ths);

void mri_inh_3d_init_guru(mri_inh_3d_plan *ths, int *N, int M, int *n,
                    int m, double sigma, unsigned nfft_flags, unsigned fftw_flags);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_3d_finalize(mri_inh_3d_plan *ths);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup texture Texture
 * This module provides the basic functions for the Texture Transforms.
 *
 * @author Matthias Schmalz
 *
 * @section texture_transforms Texture Transforms
 * In the following we describe the @ref direct_texture_transform and
 * the @ref adjoint_texture_transform.
 * For the definition of the spherical harmonics @f$ Y_l^n @f$ please see
 * @ref sh in @ref nfsft.
 *
 * @subsection direct_texture_transform Direct Texture Transform
 * The <b>Direct Texture Transform</b> is defined as follows:
 * \f[
 * \begin{array}{rcll}
 *
 * \text{\textbf{Input}} & : &
 * \text{frequency coefficients (frequencies): } &
 * \omega_{l, m, n} \in \mathbb{C} \quad \text{for }
 * l \in [0 \ldots N],\ m \in [-l \ldots l],\ n \in [-l \ldots l],\\[1ex]&&
 * \text{pole figures: } &
 * h_i \in \mathbb{S}^2 \quad \text{for }
 * i \in [1 \ldots N_1] \text{ and }\\[1ex]&&
 * \text{nodes: } &
 * r_{i, j} \in \mathbb{S}^2 \quad \text{for }
 * i \in [1 \ldots N_1],\ j \in [1 \ldots N_2].\\[1em]
 *
 * \text{\textbf{Output}} & : &
 * \text{sample values (samples): }&
 * x_{i, j} \in \mathbb{C} \quad \text{for }
 * i \in [1 \ldots N_1],\ j \in [1 \ldots N_2],
 * \text{ where } \\[1ex]&&&
 * x_{i, j} = \sum_{l = 0}^{N} \sum_{m = -l}^{l} \sum_{n = -l}^{l}
 * \omega_{l, m, n} \overline{Y_l^n(h_i)} Y_l^m(r_{i, j}).
 *
 * \end{array}
 * \f]
 *
 * @subsection adjoint_texture_transform Adjoint Texture Transform
 * The <b>Adjoint Texture Transform</b> is defined as follows:
 *
 * \f[
 * \begin{array}{rcll}
 *
 * \text{\textbf{Input}} & : &
 * \text{sample values (samples): }&
 * x_{i, j} \in \mathbb{C} \quad \text{for }
 * i \in [1 \ldots N_1],\ j \in [1 \ldots N_2], \\[1ex]&&
 * \text{pole figures: } &
 * h_i \in \mathbb{S}^2 \quad \text{for }
 * i \in [1 \ldots N_1] \text{ and }\\[1ex]&&
 * \text{nodes: } &
 * r_{i, j} \in \mathbb{S}^2 \quad \text{for }
 * i \in [1 \ldots N_1],\ j \in [1 \ldots N_2].\\[1em]
 *
 * \text{\textbf{Ouput}} & : &
 * \text{frequency coefficients (frequencies): } &
 * \omega_{l, m, n} \in \mathbb{C} \quad \text{for }
 * l \in [0 \ldots N],\ m \in [-l \ldots l],\ n \in [-l \ldots l],
 * \text{ where}\\[1ex]&&&
 * \omega_{l, m, n} = \sum_{i = 1}^{N_1} \sum_{j = 1}^{N_2}
 * x_{i, j} Y_l^n(h_i) \overline{Y_l^m(r_{i, j})}.
 *
 * \end{array}
 * \f]
 *
 * @section texture_states States of the Transformation
 * For reasons of performance this module has some state based behaviour,
 * i.e. certain functions only yield the correct result, if other certain
 * functions have been called before.
 * For ease of notation we denominate
 * - ::texture_trafo, ::texture_adjoint, ::itexture_before_loop and
 *   ::itexture_loop_one_step as <b>transform functions</b>,
 * - ::texture_precompute and ::texture_precompute_advanced as
 *   <b>precomputation functions</b> and
 * - ::texture_init and ::texture_init_advanced as
 *   <b>initialisation functions</b>.
 *
 * You have to bear in mind the following two points:
 * -# Precomputation
 *  - State 1:
 *   - The behaviour of the transform functions is undefined.
 *   - ::texture_forget has no effect.
 *   - The precomputation functions cause a state change to state 2 and
 *     initialise the precomputed data.
 *   - There is no memory allocated for precomputed data.
 *  - State 2:
 *   - The transform functions yield the correct result, if the
 *     bandwidth of the transform plan is in the valid range according to
 *     the precomputed data.
 *     Otherwise their behaviour is undefined.
 *   - The precomputation functions change the precomputed data.
 *   - ::texture_forget causes a state change to state 1 and frees all memory
 *     for precomputed data.
 *   - There is some memory allocated for precomputed data.
 * -# Manipulation of Transform Plans
 *   - Before using the transform functions with a cerain plan, you have to
 *     initialise it with one of the initialisation functions.
 *   - After the initialisation you can apply the transform functions as often
 *     as you want.
 *     But the only way you may read or manipulate elements of the plan is
 *     using the
 *     utility
 *     functions described in @ref texture_util.
 *   - After the usage of a plan you must destroy its temporary data with
 *     ::texture_finalize
 *     to avoid memory leaks.
 *
 * @section texture_data_rep Data Representation
 * Spherical coordinates are represented by two angles @f$\phi@f$ (latitude)
 * and
 * @f$\theta@f$ (longitude) as described in @ref sc.
 * Their normalisation has to be defined in TEXTURE_MAX_ANGLE before compiling
 * the
 * library.
 * Hence @f$\phi@f$ and @f$\theta@f$ have to satisfy the following conditions:
 * \f[
 * \phi \in
 * [- \frac{TEXTURE\_MAX\_ANGLE}{2}, \frac{TEXTURE\_MAX\_ANGLE}{2})
 * \quad and \quad
 * \theta \in
 * [0, \frac{TEXTURE\_MAX\_ANGLE}{2}].
 * \f]
 *
 * In the following we describe, how the input and output data
 * @f$\omega,\ x,\ h \text{ and } r@f$ is stored in the arguments for
 * ::texture_init or ::texture_init_advanced.
 * Formally the following conditions hold:
 * - @f$\omega_{l, m, n} = @f$ omega[::texture_flat_index (l, m, n)],
 * - @f$x_{i, j} = @f$ x[i * N2 + j],
 * - the latitude @f$\phi@f$ of @f$h_{i} = @f$ h_phi[i],
 * - the longitude @f$\theta@f$ of @f$h_{i} = @f$ h_theta[i],
 * - the latitude @f$\phi@f$ of @f$r_{i, j} = @f$ r[2 * (i * N2 + j)] and
 * - the longitude @f$\theta@f$ of @f$r_{i, j} = @f$ r[2 * (i * N2 + j) + 1]
 *
 * for all @f$l \in [0 \ldots N],\ m \in [-l \ldots l],\ n \in [-l \ldots l],\
 * i \in [1 \ldots N_1] \text{ and } j \in [1 \ldots N_2].@f$
 *
 * To get a better feeling what ::texture_flat_index does, see the following
 * fragment of code:
 * @code
 * int l, m, n;
 * for(l = 0; l <= N; l++) {
 *   for(m = -l; m <= l; m++) {
 *     for(n = -l; n <= l; n++) {
 *       printf("%d\n", texture_flat_index(l, m, n));
 *     }
 *   }
 * }
 * @endcode
 * It will print a list of succeeding numbers from 0.
 * @{
 */

/** @defgroup texture_private Texture: Private Functions
 * This module containes the private functions used for the implementation
 * of the texture transforms.
 * Users of the library can skip this section since it has been written for
 * developers.
 *
 * @author Matthias Schmalz
 */

/** @defgroup texture_util Texture: Utility Functions
 * This module provides functions that perform some basic operations on the
 * @ref texture_plan data structur.
 *
 * @author Matthias Schmalz
 */

/**
 * Constant for the period length of sine (default: @f$2 \pi@f$)
 */
#define TEXTURE_MAX_ANGLE (2*3.1415926535897932384)

/**
 * @addtogroup texture_private
 * @{
 */

/** Default value for texture_precompute_flags.
 * @see texture_precompute_advanced
 */
#define TEXTURE_DEF_PRECOMPUTE_FLAGS 0U

/** Default value for nfsft_precompute_flags.
 * @see texture_precompute_advanced
 */
#define TEXTURE_DEF_NFSFT_PRECOMPUTE_FLAGS 0U

/** Default value for nfsft_threshold.
 * @see texture_precompute_advanced
 */
#define TEXTURE_DEF_NFSFT_THRESHOLD 1000.0

/** Default value for texture_init_flags.
 * @see texture_init_advanced
 */
#define TEXTURE_DEF_INIT_FLAGS 0U

/** Default value for nfsft_init_flags.
 * @see texture_init_advanced
 */
#define TEXTURE_DEF_NFSFT_INIT_FLAGS 0U

/** Default value for nfft_cutoff.
 * @see texture_init_advanced
 */
#define TEXTURE_DEF_NFFT_CUTOFF 8

/**
 * @}
 */

/** @typedef texture_plan.
 * @ingroup texture
 * @brief Stores all data for a direct and adjoint transformation.
 */

/** @struct texture_plan_
 * @ingroup texture_private
 * @brief Definition of the texture_plan.
 *
 * @attention Do not access any member directly!
 * The plan could get an inconsistent state.
 */
typedef struct texture_plan_ {

  /** @var f
   * @brief Used to store the sample data x.
   * @see texture_init
   */

  /** @var f_hat
   * @brief Used to store the frequencies omega.
   * @see texture_init
   */

  /** @var N_total
   * @brief The total length of f_hat.
   */

  /** The total length of f.
   * @var M_total
   */
  MACRO_MV_PLAN(complex);

  /** The bandwidth.
   * @see texture_init
   */
  int N;

  /** The number of pole figures.
   * @see texture_init
   */
  int N1;

  /** The number of samples per pole figure.
   * @see texture_init
   */
  int N2;

  /** The latitudes of the pole figures.
   * @see texture_init
   */
  const double *h_phi;

  /** The longitudes of the pole figures.
   * @see texture_init
   */
  const double *h_theta;

  /** The nodes for the samples for each pole figure.
   * @see texture_init
   */
  const double *r;

  /** The flags for the initialisation of the nfsft.
   * @see texture_init_advanced
   */
  unsigned int nfsft_init_flags;

  /** The nfft_cutoff for the initialisation of the nfsft.
   * @see texture_init_advanced
   */
  unsigned int nfft_cutoff;

  /** Stores the cosines of the components of h_theta.
   */
  double *cos_h_theta;

  /** Stores the sines of the components of h_theta.
   */
  double *sin_h_theta;

  /** Stores the frequencies for the nfsft transformation.
   */
  complex **nfsft_f_hat;

  /** Stores the samples for the nfsft transformation.
   */
  complex *nfsft_f;

  /** Stores the nodes for the nfsft transformation.
   */
  double *nfsft_angles;
} texture_plan;

/** Performes precomputations with default values for all parameters.
 * Afterwards ::texture_trafo and ::texture_adjoint will work with any plans
 * having a bandwidth equal or less than N.
 *
 * @attention To free allocated memory ::texture_forget has to be called.
 *
 * @param N - the maximum bandwidth
 *
 * @see TEXTURE_DEF_PRECOMPUTE_FLAGS
 * @see TEXTURE_DEF_NFSFT_PRECOMPUTE_FLAGS
 * @see TEXTURE_DEF_NFSFT_THRESHOLD
 */
void texture_precompute(int N);

/** Performes precomputations.
 * Afterwards ::texture_trafo and ::texture_adjoint will work with any plans
 * having a bandwidth equal or less than N.
 *
 * @attention To free allocated memory ::texture_forget has to be called.
 * @remark Use ::texture_precompute instead if you do not know, what you are
 * doing.
 *
 * @param N - the maximum bandwidth
 * @param texture_precompute_flags - does not have any effect
 * @param nfsft_precompute_flags - flags for the precomputation of the nfsft
 * @param nfsft_threshold - a parameter for the precomputation of the nfsft
 */
void texture_precompute_advanced(int N, unsigned int texture_precompute_flags,
    unsigned int nfsft_precompute_flags, double nfsft_threshold);

/** Initialisation of a plan with default values for all parameters.
 * The arguments after ths will be stored in the plan ths.
 *
 * @par ths - Points to the transformation plan.
 * @par N - the bandwidth
 * @par N1 - the number of pole figures
 * @par N2 - the number of samples per pole figure
 * @par omega - the frequencies
 * @par x - the samples
 * @par h_phi - the latitudes of the pole figures
 * @par h_theta - the longitudes of the pole figures
 * @par r - the samples of the pole figures
 *
 * @attention
 * - ::texture_init performes only a flat copy. If you change the storage
 * referenced to by any of the pointer arguments, the plan can get an
 * inconsistent
 * state.
 * - Use ::texture_finalize to free allocated memory.
 *
 * @pre
 * All pointer arguments have to point to allocated memory.
 * omega, x, h_phi, h_theta and r have to refer to arrays of appropriate
 * lengths.
 *
 * @note
 * For details about data representation see @ref texture_data_rep.
 */
void texture_init(texture_plan *ths, int N, int N1, int N2, complex* omega,
    complex* x, const double* h_phi, const double* h_theta, const double* r);

/** Initialisation of a plan.
 * The arguments after ths will be stored in the plan ths.
 *
 * @par ths - Points to the transformation plan.
 * @par N - the bandwidth
 * @par N1 - the number of pole figures
 * @par N2 - the number of samples per pole figure
 * @par omega - the frequencies
 * @par x - the samples
 * @par h_phi - the latitudes of the pole figures
 * @par h_theta - the longitudes of the pole figures
 * @par r - the nodes of the pole figures
 * @par texture_init_flags - does not have any effect
 * @par nfsft_init_flags - flags to use for the initialisation of the nfsft
 * @par nfft_cutoff - a parameter for the initialisation of the nfsft
 *
 * @attention
 * - ::texture_init_advanced performes only a flat copy. If you change the
 *   storage
 * referenced to by any of the pointer arguments, the plan can get an
 * inconsistent
 * state.
 * - Use ::texture_finalize to free allocated memory.
 *
 * @remark Use ::texture_init instead if you do not know, what you are
 * doing.
 *
 * @pre
 * All pointer arguments have to point to allocated memory.
 * omega, x, h_phi, h_theta and r have to refer to arrays of appropriate
 * lengths.
 *
 * @note
 * For details about data representation see @ref texture_data_rep.
 */
void texture_init_advanced(texture_plan *ths, int N, int N1, int N2,
    complex* omega, complex* x, const double* h_phi, const double* h_theta,
    const double *r, unsigned int texture_init_flags,
    unsigned int nfsft_init_flags, int nfft_cutoff);

/** Carries out the direct transform.
 * Maps the frequencies on the samples.
 * Therefore the samples will be changed,
 * everything else will be preserved.
 *
 * @par ths - Points to the transformation plan.
 * @pre
 * - ::texture_precompute or ::texture_precompute_advanced must have been called
 *   with appropriate arguments.
 * - The plan hast to be initialised with ::texture_init or
 *   ::texture_init_advanced.
 */
void texture_trafo(texture_plan *ths);

/** Carries out the adjoint transform.
 * Maps the samples on the frequencies.
 * Therefor the frequencies change, everything else is
 * preserved.
 *
 * @par ths - Points to the transformation plan.
 * @pre
 * - ::texture_precompute or ::texture_precompute_advanced must have been called
 *   with appropriate arguments.
 * - The plan hast to be initialised with ::texture_init or
 *   ::texture_init_advanced.
 */
void texture_adjoint(texture_plan *ths);

/** Frees all memory allocated by ::texture_init or ::texture_init_advanced.
 */
void texture_finalize(texture_plan *ths);

/** Frees all memory allocated by ::texture_precompute or
 * ::texture_precompute_advanced.
 */
void texture_forget();

/** @addtogroup texture_util
 * @{
 */

/** Convert a non-flat index of the frequencies @f$ \omega @f$ to a
 * flat index.
 * See @ref texture_data_rep for more information.
 *
 * @arg l - the first index
 * @arg m - the second index
 * @arg n - the third index
 *
 * @return - the flat index
 *
 * @pre
 * - @f$ m \in [-l \dots l] @f$
 * - @f$ n \in [-l \dots l] @f$
 */
int texture_flat_index(int l, int m, int n);

/** Determines the length of an array omega storing frequencies in a
 * given
 * bandwidth.
 *
 * @par N - the bandwidth.
 * @return the length of the corresponding omega
 */
int texture_flat_length(int N);

/** Returns the length of the frequency array stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
int texture_get_omega_length(texture_plan *ths);

/** Returns the length of the sample array stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
int texture_get_x_length(texture_plan *ths);

/** Returns the bandwidth stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
int texture_get_N(texture_plan *ths);

/** Returns the number of pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
int texture_get_N1(texture_plan *ths);

/** Returns the number of samples per pole figure stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
int texture_get_N2(texture_plan *ths);

/** Returns a pointer to the frequencies stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
const complex *texture_get_omega(texture_plan *ths);

/** Sets the frequencies in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par omega - a pointer to the new frequencies.
 * @pre omega has to point to an array of appropriate length.
 */
void texture_set_omega(texture_plan *ths, complex* omega);

/** Returns a pointer to the samples stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
const complex *texture_get_x(texture_plan *ths);

/** Sets the samples in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par x - a pointer to the new samples
 * @pre x has to point to an array of appropriate length.
 */
void texture_set_x(texture_plan *ths, complex* x);

/** Returns a pointer to the latitudes of the pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
const double *texture_get_h_phi(texture_plan *ths);

/** Sets the latitudes of the pole figures in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par h_phi - a pointer to the new latitudes
 * @pre h_phi has to point to an array of appropriate length.
 */
void texture_set_h_phi(texture_plan *ths, const double* h_phi);

/** Returns the longitudes of the pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
const double *texture_get_h_theta(texture_plan *ths);

/** Sets the longitudes of the pole figures in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par h_theta - a pointer to the longitudes
 * @pre h_theta has to point to an array of appropriate length.
 */
void texture_set_h_theta(texture_plan *ths, const double* h_theta);

/** Returns the nodes of the pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
const double *texture_get_r(texture_plan *ths);

/** Sets the nodes of the pole figures in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par r - a pointer to the nodes
 * @pre r has to point to an array of appropriate length.
 */
void texture_set_r(texture_plan *ths, const double* r);

/** Returnes the flags used for the initialisation of the nfsft stored in a
 * plan.
 *
 * @par ths - a pointer to the transformation plan
 */
unsigned int texture_get_nfsft_init_flags(texture_plan *ths);

/** Sets the flags used for the initialisation of the nfsft in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par nfsft_init_flags - the nfsft flags
 */
void texture_set_nfsft_init_flags(texture_plan *ths,
    unsigned int nfsft_init_flags);

/** Returns the nfft_cutoff parameter used for the initialisation of the nfsft
 * stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
int texture_get_nfft_cutoff(texture_plan *ths);

/** Sets the nfft_cutoff parameter used for the initialisation of the nfsft
 * in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par nfft_cutoff - the parameter
 */
void texture_set_nfft_cutoff(texture_plan *ths, int nfft_cutoff);

/** @}
 */

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup nfsft NFSFT
 * @{
 *
 * This module implements nonuniform fast spherical Fourier transforms. In the
 * following, we abbreviate the term "nonuniform fast spherical Fourier
 * transform" by NFSFT.
 *
 * \section Preliminaries
 * This section summarises basic definitions and properties related to spherical
 * Fourier transforms.
 *
 * \subsection sc Spherical Coordinates
 * Every point in \f$\mathbb{R}^3\f$ can be described in \e spherical \e
 * coordinates by a vector \f$(r,\vartheta,\varphi)^{\mathrm{T}}\f$ with the
 * radius \f$r \in \mathbb{R}^{+}\f$ and two angles \f$\vartheta \in [0,\pi]\f$,
 * \f$\varphi \in [-\pi,\pi)\f$.
 * We denote by \f$\mathbb{S}^2\f$ the two-dimensional unit sphere embedded
 * into \f$\mathbb{R}^3\f$, i.e.
 * \f[
 *   \mathbb{S}^2 := \left\{\mathbf{x} \in \mathbb{R}^{3}:\;
 *   \|\mathbf{x}\|_2=1\right\}
 * \f]
 * and identify a point from \f$\mathbb{S}^2\f$ with the corresponding vector
 * \f$(\vartheta,\varphi)^{\mathrm{T}}\f$. The
 * spherical coordinate system is illustrated in the following figure:
 * \image html sphere.png ""
 * \image latex sphere.pdf "" width=0.45\textwidth
 * For consistency with the other modules and the conventions used there, we
 * also use \e swapped \e scaled \e spherical \e coordinates \f$x_1 :=
 * \frac{\varphi}{2\pi}\f$, \f$x_2 := \frac{\vartheta}{2\pi}\f$ and identify a
 * point from \f$\mathbb{S}^2\f$ with the vector
 * \f$\mathbf{x} := \left(x_1,x_2\right) \in
 *  [-\frac{1}{2}, \frac{1}{2}) \times [0,\frac{1}{2}]\f$.
 *
 * \subsection lp Legendre Polynomials
 * The \e Legendre \e polynomials \f$P_k : [-1,1]
 * \rightarrow \mathbb{R}$, $k \in \mathbb{N}_{0}\f$ as \e classical \e
 * orthogonal \e polynomials are given by their corresponding \e Rodrigues \e
 * formula
 * \f[
 *   P_k(t) := \frac{1}{2^k k!} \frac{\text{d}^k}{\text{d} t^k}
 *   \left(t^2-1\right)^k.
 * \f]
 * The corresponding three-term recurrence relation is
 * \f[
 *   (k+1)P_{k+1}(t) = (2k+1) x P_{k}(t) - k P_{k-1}(t) \quad (k \in
 *   \mathbb{N}_0).
 * \f]
 * With
 * \f[
 *   \left< f,g \right>_{\text{L}^2\left([-1,1]\right)} :=
 *   \int_{-1}^{1} f(t) g(t) \text{d} t
 * \f]
 * being the usual \f$\text{L}^2\left([-1,1]\right)\f$ inner product,
 * the Legendre polynomials obey the orthogonality condition
 * \f[
 *   \left< P_k,P_l \right>_{\text{L}^2\left([-1,1]\right)} = \frac{2}{2k+1}
 *   \delta_{k,l}.
 * \f]
 *
 * \remark The normalisation constant \f$ c_k := \sqrt{\frac{2k+1}{2}}\f$
 * renders the scaled Legendre polynomials \f$c_k P_k\f$ orthonormal with
 * respect to the induced \f$\text{L}^2\left([-1,1]\right)\f$ norm
 * \f[
 *   \|f\|_{\text{L}^2\left([-1,1]\right)} :=
 *   \left(<f,f>_{\text{L}^2\left([-1,1]\right)}\right)^{1/2} =
 *   \left(\int_{-1}^{1} |f(t)|^2 \; \text{d} t\right)^{1/2}.
 * \f]
 *
 * \subsection alf Associated Legendre Functions
 * The \a associated \a Legendre \a functions \f$P_k^n : [-1,1] \rightarrow
 * \mathbb{R} \f$, \f$n \in \mathbb{N}_0\f$, \f$k \ge n\f$ are defined by
 * \f[
 *   P_k^n(t) := \left(\frac{(k-n)!}{(k+n)!}\right)^{1/2}
 *   \left(1-t^2\right)^{n/2} \frac{\text{d}^n}{\text{d} t^n} P_k(t).
 * \f]
 * For \f$n = 0\f$, they coincide with the Legendre polynomials, i.e.
 * \f$P_k^0 = P_k\f$.
 * The associated Legendre functions obey the three-term recurrence relation
 * \f[
 *   P_{k+1}^n(t) = v_{k}^n t P_k^n(t) + w_{k}^n P_{k-1}^n(t) \quad (k \ge n),
 * \f]
 * with \f$P_{n-1}^n(t) := 0\f$, \f$P_{n}^n(t) := \frac{\sqrt{(2n)!}}{2^n n!}
 * \left(1-t^2\right)^{n/2}\f$, and
 * \f[
 *   v_{k}^n := \frac{2k+1}{((k-n+1)(k+n+1))^{1/2}}\; ,\qquad
 *   w_{k}^n := - \frac{((k-n)(k+n))^{1/2}}{((k-n+1)(k+n+1))^{1/2}}.
 * \f]
 * For fixed \f$n\f$, the set \f$\left\{P_k^n:\: k
 * \ge n\right\}\f$ forms a complete set of orthogonal functions in
 * \f$\text{L}^2\left([-1,1]\right)\f$
 * with
 * \f[
 *   \left< P_k^n,P_l^n \right>_{\text{L}^2\left([-1,1]\right)} = \frac{2}{2k+1}
 *   \delta_{k,l} \quad (0 \le n \le k,l).
 * \f]
 *
 * \remark The normalisation constant \f$ c_k = \sqrt{\frac{2k+1}{2}}\f$
 * renders the scaled associated Legendre functions \f$c_k P_k^n\f$ orthonormal
 * with respect to the induced \f$\text{L}^2\left([-1,1]\right)\f$ norm
 * \f[
 *   \|f\|_{\text{L}^2\left([-1,1]\right)} :=
 *   \left(<f,f>_{\text{L}^2\left([-1,1]\right)}\right)^{1/2} =
 *   \left(\int_{-1}^{1} |f(t)|^2 \; \text{d} t\right)^{1/2}.
 * \f]
 *
 * \subsection sh Spherical Harmonics
 * The standard orthogonal basis of spherical harmonics for \f$\text{L}^2
 * \left(\mathbb{S}^2\right)\f$ with yet unnormalised basis functions
 * \f$\tilde{Y}_k^n : \mathbb{S}^2 \rightarrow \mathbb{C}\f$ is given by
 * \f[
 *   \tilde{Y}_k^n(\vartheta,\varphi) := P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}
 * \f]
 * with the usual \f$\text{L}^2\left(\mathbb{S}^2\right)\f$ inner product
 * \f[
 *   \left< f,g \right>_{\mathrm{L}^2\left(\mathbb{S}^2\right)} :=
 *   \int_{\mathbb{S}^2} f(\vartheta,\varphi) \overline{g(\vartheta,\varphi)}
 *   \: \mathrm{d} \mathbf{\xi} := \int_{-\pi}^{\pi} \int_{0}^{\pi}
 *   f(\vartheta,\varphi) \overline{g(\vartheta,\varphi)} \sin \vartheta
 *   \; \mathrm{d} \vartheta \; \mathrm{d} \varphi.
 * \f]
 * The normalisation constant \f$c_k^n := \sqrt{\frac{2k+1}{4\pi}}\f$ renders
 * the  scaled basis functions
 * \f[
 *   Y_k^n(\vartheta,\varphi) := c_k^n P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}
 * \f]
 * orthonormal with respect to the induced \f$\text{L}^2\left(\mathbb{S}^2
 * \right)\f$ norm
 * \f[
 *   \|f\|_{\text{L}^2\left(\mathbb{S}^2\right)} =
 *   \left(<f,f>_{\text{L}^2\left(\mathbb{S}^2\right)}\right)^{1/2} =
 *   \left(\int_{-\pi}^{\pi} \int_{0}^{\pi} |f(\vartheta,\varphi)|^2 \sin
 *   \vartheta \; \mathrm{d} \vartheta \; \mathrm{d} \varphi\right)^{1/2}.
 * \f]
 * A function \f$f \in \mathrm{L}^2\left(\mathbb{S}^2\right)\f$ has the
 * orthogonal expansion
 * \f[
 *   f = \sum_{k=0}^{\infty} \sum_{n=-k}^{k} \hat{f}(k,n) Y_k^n,
 * \f]
 * where the coefficients \f$\hat{f}(k,n) := \left< f, Y_k^{n}
 * \right>_{\mathrm{L}^2\left(\mathbb{S}^2\right)}\f$ are the \e spherical
 * \e Fourier \e coefficients and the equivalence is understood in the
 * \f$\mathrm{L}^2\f$-sense.
 *
 *
 * \section nfsfts Nonuniform Fast Spherical Fourier Transforms
 *
 * This section describes the input and output relation of the spherical
 * Fourier transform algorithms and the layout of the corresponding plan
 * structure.
 *
 * \subsection ndsft Nonuniform Discrete Spherical Fourier Transform
 * The \e nonuniform \e discrete \e spherical \e Fourier \e transform (\e NDSFT)
 * is defined as follows:
 * \f[
 *     \begin{array}{rcl}
 *       \text{\textbf{Input}} & : & \text{coefficients }
 *         \hat{f}(k,n) \in \mathbb{C} \text{ for } k=0,\ldots,N,\;n=-k,
 *         \ldots,k,\; N \in \mathbb{N}_0,\\[1ex]
 *                             &   & \text{arbitrary nodes } \mathbf{x}(m) \in
 *         [-\frac{1}{2},\frac{1}{2}] \times [0,\frac{1}{2}]
 *         \text{ for } m=0,\ldots,M-1, M \in \mathbb{N}. \\[1ex]
 *       \text{\textbf{Task}}  & : & \text{evaluate } f(m) := f\left(
 *       \mathbf{x}(m)\right) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}_k^n
 *         Y_k^n\left(\mathbf{x}(m)\right) \text{ for } m=0,\ldots,M-1.
 *         \\[1ex]
 *       \text{\textbf{Output}} & : & \text{coefficients } f(m) \in
 *         \mathbb{C} \text{ for } m=0,\ldots,M-1.\\
 *     \end{array}
 * \f]
 *
 * \subsection andsft Adjoint Nonuniform Discrete Spherical Fourier Transform
 * The \e adjoint \e nonuniform \e discrete \e spherical \e Fourier \e transform
 * (\e adjoint \e NDSFT)
 * is defined as follows:
 * \f[
 *     \begin{array}{rcl}
 *       \text{\textbf{Input}} & : & \text{coefficients } f(m) \in
 *         \mathbb{C} \text{ for } m=0,\ldots,M-1, M \in \mathbb{N},\\
 *                             &   & \text{arbitrary nodes } \mathbf{x}(m) \in
 *         [-\frac{1}{2},\frac{1}{2}] \times [0,\frac{1}{2}] \text{ for }
 *         m=0,\ldots,M-1, N \in \mathbb{N}_0.\\[1ex]
 *       \text{\textbf{Task}}  & : & \text{evaluate } \hat{f}(k,n)
 *         := \sum_{m=0}^{M-1} f(m) \overline{Y_k^n\left(\mathbf{x}(m)\right)}cd Do
 *         \text{ for } k=0,\ldots,N,\;n=-k,\ldots,k.\\[1ex]
 *       \text{\textbf{Output}} & : & \text{coefficients }
 *         \hat{f}(k,n) \in \mathbb{C} \text{ for }
 *         k=0,\ldots,N,\;n=-k,\ldots,k.\\[1ex]
 *     \end{array}
 * \f]
 *
 * \subsection dl Data Layout
 * This section describes the public  layout of the \ref nfsft_plan structure
 * which
 * contains all data for the computation of the aforementioned spherical Fourier
 * transforms. The structure contains private (no read or write allowed), public
 * read-only (only
 * read access permitted), and public read-write (read and write access allowed)
 * members. In the following, we indicate read and write access by \c read and
 * \c write. The public members are structured as follows:
 * \li \c N_total (\c read)
 *        The total number of components in \c f_hat. If the bandwidth is
 *        \f$N \in \mathbb{N}_0\f$, the total number of components in \c f_hat
 *        is \c N_total \f$ = (2N+2)^2\f$.
 * \li \c M_total (\c read)
 *        the total number of samples \f$M\f$
 * \li \c f_hat (\c read-write)
 *        The flattened array of spherical Fourier coefficents. The array
 *        has length \f$(2N+2)^2\f$ such that valid indices \f$i \in
 *        \mathbb{N}_0\f$ for array access \c f_hat \c[ \f$i\f$ \c] are
 *        \f$i=0,1,\ldots,(2N+2)^2-1\f$.
 *        However, only read and write access to indices corresponding to
 *        spherical Fourier coefficients \f$\hat{f}(k,n)\f$ is defined. The index
 *        \f$i\f$ corresponding to the spherical Fourier coefficient
 *        \f$\hat{f}(k,n)\f$ with \f$0 \le k \le M\f$, \f$-k \le n \le k\f$ is
 *        \f$i = (N+2)(N-n+1)+N+k+1\f$. For convenience, the helper macro
 *        \ref NFSFT_INDEX(k,n) provides the necessary index calculations such
 *        that
 *        one can write \c f_hat[ \c NFSFT_INDEX(\f$k,n\f$\c)] \c =
 *        \c ... to access
 *        the component corresponding to \f$\hat{f}(k,n)\f$.
 *        The data layout is due to implementation details.
 * \li \c f (\c read-write)
 *        the array of coefficients \f$f(m)\f$ for \f$m=0,\ldots,M-1\f$ such
 *        that \c f[\f$m\f$\c] = \f$f(m)\f$
 * \li \c N (\c read)
 *        the bandwidth \f$N \in \mathbb{N}_0\f$
 * \li \c x
 *        the array of nodes \f$\mathbf{x}(m) \in
 *        [-\frac{1}{2},\frac{1}{2}] \times [0,\frac{1}{2}]\f$ for \f$m = 0,
 *        \ldots,M-1\f$ such that \c f[\f$2m\f$\c] = \f$x_1\f$ and
 *        \c f[\f$2m+1\f$\c] = \f$x_2\f$
 *
 * \subsection gtn Good to know...
 * When using the routines of this module you should bear in mind the following:
 * \li The bandwidth \f$N_{\text{max}}\f$ up to which precomputation is
 *   performed is always chosen as the next power of two with respect to the
 *   specified maximum bandwidth.
 * \li By default, the NDSFT transforms (see \ref ndsft_trafo, \ref nfsft_trafo)
 *   are allowed to destroy the input \c f_hat while the input \c x is
 *   preserved. On the contrary, the adjoint NDSFT transforms
 *   (see \ref ndsft_adjoint, \ref nfsft_adjoint) do not destroy the input
 *   \c f and \c x by default. The desired behaviour can be assured by using the
 *   \ref NFSFT_PRESERVE_F_HAT, \ref NFSFT_PRESERVE_X, \ref NFSFT_PRESERVE_F and
 *   \ref NFSFT_DESTROY_F_HAT, \ref NFSFT_DESTROY_X, \ref NFSFT_DESTROY_F
 *   flags.
 */

/* Planner flags */

/**
 * By default, all computations are performed with respect to the
 * unnormalized basis functions
 * \f[
 *   \tilde{Y}_k^n(\vartheta,\varphi) = P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}.
 * \f]
 * If this flag is set, all computations are carried out using the \f$L_2\f$-
 * normalized basis functions
 * \f[
 *   Y_k^n(\vartheta,\varphi) = \sqrt{\frac{2k+1}{4\pi}} P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}.
 * \f]
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_NORMALIZED    (1U << 0)

/**
 * If this flag is set, the fast NFSFT algorithms (see \ref nfsft_trafo,
 * \ref nfsft_adjoint) will use internally the exact but usually slower direct
 * NDFT algorithm in favor of fast but approximative NFFT algorithm.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_USE_NDFT      (1U << 1)

/**
 * If this flag is set, the fast NFSFT algorithms (see \ref nfsft_trafo,
 * \ref nfsft_adjoint) will use internally the usually slower direct
 * DPT algorithm in favor of the fast FPT algorithm.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 * \warning This feature is not implemented yet!
 */
#define NFSFT_USE_DPT       (1U << 2)

/**
 * If this flag is set, the init methods (see \ref nfsft_init , \ref
 * nfsft_init_advanced , and \ref nfsft_init_guru) will allocate memory and the
 * method \ref nfsft_finalize will free the array \c x for you. Otherwise,
 * you have to assure by yourself that \c x points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_MALLOC_X      (1U << 3)

/**
 * If this flag is set, the init methods (see \ref nfsft_init , \ref
 * nfsft_init_advanced , and \ref nfsft_init_guru) will allocate memory and the
 * method \ref nfsft_finalize will free the array \c f_hat for you. Otherwise,
 * you have to assure by yourself that \c f_hat points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_MALLOC_F_HAT  (1U << 5)

/**
 * If this flag is set, the init methods (see \ref nfsft_init , \ref
 * nfsft_init_advanced , and \ref nfsft_init_guru) will allocate memory and the
 * method \ref nfsft_finalize will free the array \c f for you. Otherwise,
 * you have to assure by yourself that \c f points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_MALLOC_F      (1U << 6)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref ndsft_trafo or \ref nfsft_trafo the content of \c f_hat remains
 * unchanged.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_PRESERVE_F_HAT (1U << 7)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref ndsft_trafo, \ref nfsft_trafo or \ref ndsft_adjoint, \ref nfsft_adjoint
 * the content of \c x remains
 * unchanged.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_PRESERVE_X     (1U << 8)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref ndsft_adjoint or \ref nfsft_adjoint the content of \c f remains
 * unchanged.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_PRESERVE_F     (1U << 9)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref ndsft_trafo or \ref nfsft_trafo the content of \c f_hat may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_F_HAT    (1U << 10)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref ndsft_trafo, \ref nfsft_trafo or \ref ndsft_adjoint, \ref nfsft_adjoint
 * the content of \c x may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_X      (1U << 11)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref ndsft_adjoint or \ref nfsft_adjoint the content of \c f may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_F      (1U << 12)

/* Precomputation flags */

/**
 * If this flag is set, the transforms \ref ndsft_trafo and \ref ndsft_adjoint
 * do not work. Setting this flag saves some memory for precomputed data.
 *
 * \see nfsft_precompute
 * \see ndsft_trafo
 * \see ndsft_adjoint
 * \author Jens Keiner
 */
#define NFSFT_NO_DIRECT_ALGORITHM    (1U << 13)

/**
 * If this flag is set, the transforms \ref nfsft_trafo and \ref nfsft_adjoint
 * do not work. Setting this flag saves memory for precomputed data.
 *
 * \see nfsft_precompute
 * \see nfsft_trafo
 * \see nfsft_adjoint
 * \author Jens Keiner
 */
#define NFSFT_NO_FAST_ALGORITHM      (1U << 14)

/**
 * If this flag is set, the fast transforms \ref nfsft_trafo and
 * \ref nfsft_adjoint only
 * work in a defined bandwidth window. If \f$N_{\text{max}}\f$ is the power of
 * two up to which precomputation is performed, only fast transformations for
 * bandwidths \f$N\f$ with \f$N_{\text{max}}/2 < N \le N_{\text{max}}\f$ will
 * work. The direct but usually slow transforms \ref ndsft_trafo and
 * \ref ndsft_adjoint are unaffected. Setting this flag saves memory for
 * precomputed data.
 *
 * \see nfsft_precompute
 * \see nfsft_trafo
 * \see nfsft_adjoint
 * \author Jens Keiner
 */
#define NFSFT_BANDWIDTH_WINDOW       (1U << 15)

/**
 * If this flag is set, the transforms \ref nfsft_adjoint and
 * \ref ndsft_adjoint set all unused entries in \c f_hat not corresponding to
 * spherical Fourier coefficients to zero.
 *
 * \author Jens Keiner
 */
#define NFSFT_ZERO_F_HAT             (1U << 16)

/* */

/**
 * This helper macro expands to the index \f$i\f$
 * corresponding to the spherical Fourier coefficient
 * \f$f_hat(k,n)\f$ for \f$0 \le k \le N\f$, \f$-k \le n \le k\f$ with
 * \f[
 *   (N+2)(N-n+1)+N+k+1
 * \f]
 */
#define NFSFT_INDEX(k,n,plan)        ((2*(plan)->N+2)*((plan)->N-n+1)+(plan)->N+k+1)

/**
 * This helper macro expands to the logical size of a spherical Fourier coefficients
 * array for a bandwidth N.
 */
#define NFSFT_F_HAT_SIZE(N)          ((2*N+2)*(2*N+2))

/** Structure for a NFSFT transform plan */
typedef struct nfsft_plan_
{
  /** Inherited public members */
  MACRO_MV_PLAN(complex);

  /* Public members */
  int N;                              /**< the bandwidth \f$N\f$              */
  double *x;                          /**< the nodes \f$\mathbf{x}(m) =       *
                                           \left(x_1,x_2\right) \in           *
                                           [-\frac{1}{2},\frac{1}{2}) \times  *
                                           [0,\frac{1}{2}]\f$ for             *
                                           \f$m=0,\ldots,M-1\f$,\f$M \in      *
                                           \mathbb{N},\f$                     */

  /* Private members */
  /*int NPT;*/                        /**< the next greater power of two with *
                                           respect to \f$N\f$                 */
  int t;                              /**< the logaritm of NPT with           *
                                           respect to the basis 2             */
  unsigned int flags;                 /**< the planner flags                  */
  nfft_plan plan_nfft;                /**< the internal NFFT plan             */
  complex *f_hat_intern;              /**< Internally used pointer to         *
                                           spherical Fourier coefficients     */
} nfsft_plan;

/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \ref nfsft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 *
 * \author Jens Keiner
 */
void nfsft_init(nfsft_plan *plan, int N, int M);

/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \verbatim nfsft_plan \endverbatim structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsft_flags the NFSFT flags
 *
 * \author Jens Keiner
 */
void nfsft_init_advanced(nfsft_plan* plan, int N, int M, unsigned int
                         nfsft_flags);

/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \verbatim nfsft_plan \endverbatim structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsft_flags the NFSFT flags
 * \arg nfft_cutoff the NFFT cutoff parameter
 *
 * \author Jens Keiner
 */
void nfsft_init_guru(nfsft_plan *plan, int N, int M, unsigned int nfsft_flags,
                     int nfft_cutoff);

/**
 * Performes precomputation up to the next power of two with respect to a given
 * bandwidth \f$N \in \mathbb{N}_2\f$. The threshold parameter \f$\kappa \in
 * \mathbb{R}^{+}\f$ determines the number of stabilization steps computed in
 * the discrete polynomial transform and thereby its accuracy.
 *
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg threshold the threshold \f$\kappa \in \mathbb{R}^{+}\f$
 * \arg nfsft_precomputation_flags the precomputation flags
 *
 * \author Jens Keiner
 */
void nfsft_precompute(int N, double kappa,
                      unsigned int flags);

/**
 * Forgets all precomputed data.
 *
 * \author Jens Keiner
 */
void nfsft_forget();

/**
 * Executes a direct NDSFT, i.e. computes for \f$m = 0,\ldots,M-1\f$
 * \f[
 *   f(m) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}(k,n) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void ndsft_trafo(nfsft_plan* plan);

/**
 * Executes a direct adjoint NDSFT, i.e. computes for \f$k=0,\ldots,N;
 * n=-k,\ldots,k\f$
 * \f[
 *   \hat{f}(k,n) = \sum_{m = 0}^{M-1} f(m) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void ndsft_adjoint(nfsft_plan* plan);

/**
 * Executes a NFSFT, i.e. computes for \f$m = 0,\ldots,M-1\f$
 * \f[
 *   f(m) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}(k,n) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void nfsft_trafo(nfsft_plan* plan);

/**
 * Executes an adjoint NFSFT, i.e. computes for \f$k=0,\ldots,N;
 * n=-k,\ldots,k\f$
 * \f[
 *   \hat{f}(k,n) = \sum_{m = 0}^{M-1} f(m) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void nfsft_adjoint(nfsft_plan* plan);

/**
 * Destroys a plan.
 *
 * \arg plan the plan to be destroyed
 *
 * \author Jens Keiner
 */
void nfsft_finalize(nfsft_plan* plan);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup fpt FPT
 * @{
 *
 * This module implements fast polynomial transforms. In the following, we
 * abbreviate the term "fast polynomial transforms" by FPT.
 */

/* Flags for fpt_init() */
#define FPT_NO_STABILIZATION  (1U << 0) /**< If set, no stabilization will be
                                             used.                            */
#define FPT_BANDWIDTH_WINDOW  (1U << 1) /**< If set, TODO complete comment.   */
#define FPT_NO_FAST_TRANSFORM (1U << 2) /**< If set, TODO complete comment.   */
#define FPT_NO_SLOW_TRANSFORM (1U << 3) /**< If set, TODO complete comment.   */
#define FPT_PERSISTENT_DATA   (1U << 4) /**< If set, TODO complete comment.   */

/* Flags for fpt_trafo(), dpt_transposed(), fpt_trafo(), fpt_transposed() */
#define FPT_FUNCTION_VALUES   (1U << 5) /**< If set, the output are function
                                             values at Chebyshev nodes rather
                                             than Chebyshev coefficients.     */
#define FPT_ODD_EVEN_SYMMETRY (1U << 6) /**< TODO Don't use this flag!        */
#define FPT_NEW_STABILIZATION (1U << 7) /**< TODO Don't use this flag!        */

/* Data structures */
typedef struct fpt_set_s_ *fpt_set;     /**< A set of precomputed data for a
                                             set of DPT transforms of equal
                                             maximum length.                  */

/**
 * Initializes a set of precomputed data for DPT transforms of equal length.
 *
 * \arg M The maximum DPT transform index \f$M \in \mathbb{N}_0\f$. The
 *        individual transforms are addressed by and index number \f$m \in
 *        \mathbb{N}_0\f$ with range \f$m = 0,\ldots,M\f$. The total number
 *        of transforms is therefore \f$M+1\f$.
 * \arg t The exponent \f$t \in \mathbb{N}, t \ge 2\f$ of the transform length
 *        \f$N = 2^t \in \mathbb{N}, N \ge 4\f$
 * \arg flags A bitwise combination of the flags FPT_NO_STABILIZATION,
 *            FPT_BANDWIDTH_WINDOW
 *
 * \author Jens Keiner
 */
fpt_set fpt_init(const int M, const int t, const unsigned int flags);

/**
 * Computes the data required for a single DPT transform.
 *
 * \arg set The set of DPT transform data where the computed data will be stored.
 * \arg m The transform index \f$m \in \mathbb{N}, 0 \le m \le M\f$.
 * \arg alpha The three-term recurrence coefficients \f$\alpha_k \in
 *      \mathbb{R}\f$ for \f$k=0,\ldots,N\f$ such that \verbatim alpha[k]
 *      \endverbatim \f$=\alpha_k\f$.
 * \arg beta The three-term recurrence coefficients \f$\beta_k \in \mathbb{R}\f$
 *            for \f$k=0,\ldots,N\f$ such that \verbatim beta[k] \endverbatim
 *            \f$=\beta_k\f$.
 * \arg gamma The three-term recurrence coefficients \f$\gamma_k \in
 *            \mathbb{R}\f$ for \f$k=0,\ldots,N\f$ such that \verbatim gamma[k]
 *            \endverbatim \f$=\gamma_k\f$.
 * \arg k_start The index \f$k_{\text{start}} \in \mathbb{N}_0,
 *              0 \le k_{\text{start}} \le N\f$
 * \arg threshold The treshold \f$\kappa \in \mathbb{R}, \kappa > 0\f$.
 *
 * \author Jens Keiner
 */
void fpt_precompute(fpt_set set, const int m, const double *alpha,
                    const double *beta, const double *gamma, int k_start,
                    const double threshold);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void dpt_trafo(fpt_set set, const int m, const complex *x, complex *y,
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void fpt_trafo(fpt_set set, const int m, const complex *x, complex *y,
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void dpt_transposed(fpt_set set, const int m, complex *x, const complex *y,
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void fpt_transposed(fpt_set set, const int m, complex *x, const complex *y,
  const int k_end, const unsigned int flags);

void fpt_finalize(fpt_set set);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup solver Group
 * @{
 */

/**
 * Constant symbols for precomputation and memory usage (inverse problem)
 */
#define LANDWEBER             (1U<< 0)
#define STEEPEST_DESCENT      (1U<< 1)
#define CGNR                  (1U<< 2)
#define CGNE                  (1U<< 3)
#define NORMS_FOR_LANDWEBER   (1U<< 4)
#define PRECOMPUTE_WEIGHT     (1U<< 5)
#define PRECOMPUTE_DAMP       (1U<< 6)
/** will come in again, as method by their own, together with some new method for
    the minimal seminorm interpolation
    #define REGULARIZE_CGNR       (1U<< 7)
    #define REGULARIZE_CGNR_R_HAT   (1U<< 8)
 double *r_hat;                        regularisation weights
  double lambda;                        regularisation parameter
*/

/**
 * \brief function mangling macro
 *
 * \arg MV matrix vector multiplication type
 * \arg FLT float type, i.e., double resp. complex
 * \arg name name of the function
 * \arg ...  argument list of the function
 */
#define F(MV, FLT, name, ...) void i ## MV ## _ ## name(__VA_ARGS__)

#define MACRO_SOLVER_PLAN(MV, FLT)                                            \
typedef struct i ## MV ## _plan_                                              \
{                                                                        \
  MV ## _plan *mv;                      /**< matrix vector multiplication   */\
  unsigned flags;      /**< iteration type, ...            */\
                                                                              \
  double *w;                           /**< weighting factors              */\
  double *w_hat;                       /**< damping factors                */\
                                                                              \
  FLT *y;                               /**< right hand side, samples       */\
                                                                              \
  FLT *f_hat_iter;                      /**< iterative solution             */\
                        \
  FLT *r_iter;              /**< iterated residual vector       */\
  FLT *z_hat_iter;                       /**< residual vector of normal eq.1 */\
  FLT *p_hat_iter;                       /**< search direction               */\
  FLT *v_iter;                           /**< residual vector update         */\
                        \
  double alpha_iter;                    /**< step size for search direction */\
  double beta_iter;                     /**< step size for search correction*/\
                        \
  double dot_r_iter;                    /**< dotproductc{_w}(r_iter)        */\
  double dot_r_iter_old;                /**< old dotproductc{_w}(r_iter)    */\
  double dot_z_hat_iter;                /**< dotproductc{_w}(z_hat_iter)    */\
  double dot_z_hat_iter_old;            /**< old dotproductc{_w}(z_hat_iter)*/\
  double dot_p_hat_iter;                /**< dotproductc{_w}(p_hat_iter)    */\
  double dot_v_iter;                    /**< dotproductc{_w}(v_iter)        */\
} i ## MV ## _plan;                                                           \
                        \
F(MV, FLT, init,    i ## MV ## _plan *ths, MV ## _plan *mv);            \
F(MV, FLT, init_advanced, i ## MV ## _plan *ths, MV ## _plan *mv,             \
             unsigned i ## MV ## _flags);                  \
F(MV, FLT, before_loop,   i ## MV ## _plan *ths);                             \
F(MV, FLT, loop_one_step, i ## MV ## _plan *ths);                             \
F(MV, FLT, finalize,      i ## MV ## _plan *ths);                             \


MACRO_SOLVER_PLAN(nfft, complex)
MACRO_SOLVER_PLAN(nfct, double)
MACRO_SOLVER_PLAN(nfst, double)
MACRO_SOLVER_PLAN(nnfft, complex)
MACRO_SOLVER_PLAN(mri_inh_2d1d, complex)
MACRO_SOLVER_PLAN(mri_inh_3d, complex)
MACRO_SOLVER_PLAN(nfsft, complex)
MACRO_SOLVER_PLAN(texture, complex)
/** @}
 */

#endif
/* nfft3.h */
