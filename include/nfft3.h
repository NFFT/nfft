/*! \file nfft3.h
 *  \brief Header file for the nfft3 library.
 */
#ifndef NFFT3_H
#define NFFT3_H

#include <complex.h>
#include <fftw3.h>

/** 
 * Constant symbols for precomputation and memory usage
 */
#define PRE_PHI_HUT      (1U<< 0)
#define PRE_LIN_PSI      (1U<< 1)
#define PRE_PSI          (1U<< 2)
#define PRE_FULL_PSI     (1U<< 3)
#define MALLOC_X         (1U<< 4)
#define MALLOC_F_HAT     (1U<< 5)
#define MALLOC_F         (1U<< 6)
#define FFT_OUT_OF_PLACE (1U<< 7)
#define FFTW_INIT        (1U<< 8)

#define MALLOC_V         (1U<< 9)

#define SNDFT            (1U<< 10)

#define MACRO_MV_PLAN(float_type)			                      \
  int N_total;                          /**< total number of Fourier coeffs.*/\
  int M_total;                          /**< total number of samples        */\
  float_type *f_hat;                    /**< Fourier coefficients           */\
  float_type *f;                        /**< samples                        */\
  
/** @defgroup nfft Group
 * Direct and fast computation of the
 * discrete Fourier transform at nonequispaced knots
 * @{ 
 */
typedef struct nfft_plan_ 
{
  /** api */
  MACRO_MV_PLAN(complex);          

  int d;                                /**< dimension, rank                 */
  int *N;                               /**< cut-off-frequencies             */
  double *sigma;	                /**< oversampling-factor             */
  int *n;                               /**< fftw-length = sigma*N           */
  int n_total;                          /**< total size of fftw              */
  int m;                                /**< cut-off parameter in time-domain*/
  double *b;                            /**< shape parameters                */
  int K;                                /**< number of precomp. uniform psi  */

  unsigned nfft_flags;                  /**< flags for precomputation, malloc*/
  unsigned fftw_flags;                  /**< flags for the fftw              */

  double *x;                            /**< nodes (in time/spatial domain)  */

  /** internal*/
  fftw_plan  my_fftw_plan1;             /**< fftw_plan forward               */
  fftw_plan  my_fftw_plan2;		/**< fftw_plan backward              */

  double **c_phi_inv;                   /**< precomputed data, matrix D      */
  double *psi;				/**< precomputed data, matrix B      */
  int *psi_index_g;                     /**< only for PRE_FULL_PSI           */
  int *psi_index_f;                     /**< only for PRE_FULL_PSI           */

  complex *g;
  complex *g_hat;
  complex *g1;                          /**< input of fftw                   */
  complex *g2;                          /**< output of fftw                  */

  double *spline_coeffs;          	/**< input for de Boor algorithm, if  
					     B_SPLINE or SINC_2m is defined  */
} nfft_plan;


/** see formula (1.1), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 */
void ndft_trafo(nfft_plan *ths);

/** see formula (1.1), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 */
void ndft_conjugated(nfft_plan *ths);

/** see formula (1.2), computes
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void ndft_adjoint(nfft_plan *ths);

/** see formula (1.2), computes
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(-2(pi) k x[j])
 */
void ndft_transposed(nfft_plan *ths);

/** see formula (1.1), computes in a fast and approximative way
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 */
void nfft_trafo(nfft_plan *ths);

/** see formula (1.1), computes in a fast and approximative way
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 */
void nfft_conjugated(nfft_plan *ths);

/** see formula (1.2), computes in a fast and approximative way
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void nfft_adjoint(nfft_plan *ths);

/** see formula (1.2), computes in a fast and approximative way
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void nfft_transposed(nfft_plan *ths);

/** wrapper for nfft_init and d=1
 */
void nfft_init_1d(nfft_plan *ths, int N1, int M);

/** wrapper for nfft_init and d=2
 */
void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M);

/** wrapper for nfft_init and d=3
 */
void nfft_init_3d(nfft_plan *ths, int N1, int N2, int N3, int M);

/** initialisation of transform, simple interface
 */
void nfft_init(nfft_plan *ths, int d, int *N, int M);

/** initialisation of transform, advanced interfacegoogle.de/
 */
void nfft_init_advanced(nfft_plan *ths, int d, int *N, int M,
			unsigned nfft_flags_on, unsigned nfft_flags_off);

/** initialisation of transform, guru interface
 */
void nfft_init_guru(nfft_plan *ths, int d, int *N, int M, int *n,
	            int m, unsigned nfft_flags, unsigned fftw_flags);

/** precomputes equally spaced values of psi
 *  if PRE_LIN_PSI is set the application program has to call this routine
 */
void nfft_precompute_lin_psi(nfft_plan *ths);

/** precomputes the values psi in a tensor product form
 *  if PRE_PSI is set the application program has to call this routine
 *  after setting the nodes x
 */
void nfft_precompute_psi(nfft_plan *ths);

/** precomputes the values psi and their indices in non tensor product form
 *  if PRE_FULL_PSI is set the application program has to call this routine
 */
void nfft_precompute_full_psi(nfft_plan *ths);

/** finalisation of transform
 */
void nfft_finalize(nfft_plan *ths);

/* @} 
 */

 
/** @defgroup nnfft Group
 * Direct and fast computation of the
 * discrete Fourier transform at nonequispaced knots in time- and fourier-domain
 * @{ 
 */
typedef struct nnfft_plan_ 
{
  /** api */
  MACRO_MV_PLAN(complex);

  int d;                                /**< dimension, rank                 */
  double *sigma;                        /**< oversampling-factor 1           */
  double *a;                            /**< 1 + 2*m1/N1                     */
  int *N;                               /**< cut-off-frequencies             */
  int *N1;
  int *aN1;
  int m;                                /**< cut-off parameter in time-domain*/
  double *b;                            /**< shape parameters                */
  int K;                                /**< number of precomp. uniform psi  */

  int aN1_total;
  
  nfft_plan *direct_plan;               /**< plan for the nfft               */
  unsigned nnfft_flags;                 /**< flags for precomputation, malloc*/
  int *n;                               /**<  = a*N                          */
  
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


void nndft_trafo(nnfft_plan *ths_plan);

void nndft_conjugated(nnfft_plan *ths_plan);

void nndft_adjoint(nnfft_plan *ths_plan);

void nndft_transposed(nnfft_plan *ths_plan);

void nnfft_init(nnfft_plan *ths_plan, int d, int N_total, int M_total, int *N);

void nnfft_init_guru(nnfft_plan *ths_plan, int d, int N_total, int M_total,
                     int *N, int *N1, int m, unsigned nnfft_flags);
                        
void nnfft_finalize(nnfft_plan *ths_plan);

void nnfft_trafo(nnfft_plan *ths_plan);

void nnfft_conjugated(nnfft_plan *ths_plan);

void nnfft_adjoint(nnfft_plan *ths_plan);

void nnfft_transposed(nnfft_plan *ths_plan);

void nnfft_precompute_lin_psi(nnfft_plan *ths_plan);

void nnfft_precompute_psi(nnfft_plan *ths_plan);

void nnfft_precompute_full_psi(nnfft_plan *ths_plan);

void nnfft_precompute_phi_hut(nnfft_plan *ths_plan);
 
/* @} 
 */

/** @defgroup nsfft Group
 * @{ 
 */
typedef struct nsfft_plan_ 
{
  MACRO_MV_PLAN(complex);

  int d;                                /**< d=2,3                           */

  int J;                      
  int sigma;

  unsigned flags;

  int *index_sparse_to_full;            /**< overflow for d=3, J=9!          */

  int r_act_nfft_plan;
  nfft_plan *act_nfft_plan;
  nfft_plan *center_nfft_plan;

  fftw_plan* set_fftw_plan1;
  fftw_plan* set_fftw_plan2;

  nfft_plan *set_nfft_plan_1d;
  nfft_plan *set_nfft_plan_2d;          /**< d=3                             */

  double *x_transposed;                 /**< d=2                             */
  double *x_102,*x_201,*x_120,*x_021;   /**< d=3                             */
  
} nsfft_plan;

void nsdft_trafo(nsfft_plan *ths);
void nsdft_adjoint(nsfft_plan *ths);

void nsfft_trafo(nsfft_plan *ths);
void nsfft_adjoint(nsfft_plan *ths);

/** Copy coefficients and knots from ths to ths_nfft */
void nsfft_cp(nsfft_plan *ths, nfft_plan *ths_nfft);
void nsfft_init_random_nodes_coeffs(nsfft_plan *ths);

void nsfft_init(nsfft_plan *ths, int d, int J, int M, int m, unsigned flags);
void nsfft_finalize(nsfft_plan *ths);

/* @} 
 */

/** @defgroup mri_inh Group
 * @{ 
 */

typedef struct mri_inh_2d1d_plan_
{
  /** api */
  MACRO_MV_PLAN(complex);

  int d;                                /**< dimension, rank                 */
  int *N;                               /**< cut-off-frequencies             */
  double *sigma;                        /**< oversampling-factor             */
  int *n;                               /**< fftw-length = sigma*N           */
  int n_total;                          /**< total size of fftw              */
  int m;                                /**< cut-off parameter in time-domain*/
  double *b;                            /**< shape parameters                */
  int K;                                /**< number of precomp. uniform psi  */

  unsigned nfft_flags;                  /**< flags for precomputation, malloc*/
  unsigned fftw_flags;                  /**< flags for the fftw              */

  double *x;                            /**< nodes (in time/spatial domain)  */

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

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined  */
  int N3;
  double *t;
  double *w;
} mri_inh_2d1d_plan;

typedef struct mri_inh_3d_plan_
{
  /** api */
  MACRO_MV_PLAN(complex);

  nfft_plan plan;
  
  int N3;
  double *t;
  double *w;
} mri_inh_3d_plan;


void mri_inh_2d1d_trafo(mri_inh_2d1d_plan *ths);

void mri_inh_2d1d_adjoint(mri_inh_2d1d_plan *ths);

void mri_inh_2d1d_init_guru(mri_inh_2d1d_plan *ths, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags);

void mri_inh_2d1d_finalize(mri_inh_2d1d_plan *ths);

void mri_inh_3d_trafo(mri_inh_3d_plan *ths);

void mri_inh_3d_adjoint(mri_inh_3d_plan *ths);

void mri_inh_3d_init_guru(mri_inh_3d_plan *ths, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags);

void mri_inh_3d_finalize(mri_inh_3d_plan *ths);
/* @} 
 */

/** @defgroup texture group
 * @{
 */

/**
 * The structure for the transform plan.
 */
typedef struct texture_plan_ {

	MACRO_MV_PLAN(complex);

	int N;
	int N1;
	int N2;
	double *h_phi;
	double *h_theta;
	double **r_phi;
	double **r_theta;
} texture_plan;

/**
 * Carries out the transform.
 */
void texture_trafo(texture_plan *ths);

/**
 * The adjoint version of the transform.
 */
void texture_adjoint(texture_plan *ths);

/* @}
 */

/** @defgroup solver Group
 * @{ 
 */

/** 
 * Constant symbols for precomputation and memory usage (inverse problem)
 */
#define LANDWEBER           	(1U<< 0)
#define STEEPEST_DESCENT    	(1U<< 1)
#define CGNR                    (1U<< 2)
#define CGNE 		    	(1U<< 3)
#define NORMS_FOR_LANDWEBER 	(1U<< 4)
#define PRECOMPUTE_WEIGHT   	(1U<< 5)
#define PRECOMPUTE_DAMP     	(1U<< 6)
/** will come in again, as method by their own, together with some new method for
    the minimal seminorm interpolation
    #define REGULARIZE_CGNR     	(1U<< 7)
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

#define MACRO_SOLVER_PLAN(MV, FLT)	          		              \
typedef struct i ## MV ## _plan_				              \
{									      \
  MV ## _plan *mv;			/**< matrix vector multiplication   */\
  unsigned flags;			/**< iteration type, ...            */\
									      \
  double *w;                         	/**< weighting factors              */\
  double *w_hat;                     	/**< damping factors                */\
									      \
  FLT *y;                	        /**< right hand side, samples       */\
									      \
  FLT *f_hat_iter;		        /**< iterative solution             */\
									      \
  FLT *r_iter;			        /**< iterated residual vector       */\
  FLT *z_hat_iter;             	        /**< residual vector of normal eq.1 */\
  FLT *p_hat_iter;             	        /**< search direction               */\
  FLT *v_iter;                 	        /**< residual vector update         */\
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
F(MV, FLT, init,	  i ## MV ## _plan *ths, MV ## _plan *mv);            \
F(MV, FLT, init_advanced, i ## MV ## _plan *ths, MV ## _plan *mv,             \
		     	  unsigned i ## MV ## _flags);		              \
F(MV, FLT, before_loop,   i ## MV ## _plan *ths);                             \
F(MV, FLT, loop_one_step, i ## MV ## _plan *ths);                             \
F(MV, FLT, finalize,      i ## MV ## _plan *ths);                             \


MACRO_SOLVER_PLAN(nfft, complex)
MACRO_SOLVER_PLAN(nnfft, complex)
MACRO_SOLVER_PLAN(mri_inh_2d1d, complex)
MACRO_SOLVER_PLAN(mri_inh_3d, complex)
MACRO_SOLVER_PLAN(texture, complex)
/** @} 
 */
 
#endif
/* nfft3.h */
