#ifndef solver_h_inc
#define solver_h_inc

#include <complex.h>
#include "nfft3.h"

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
#define REGULARIZE_CGNR     	(1U<< 7)
#define REGULARIZE_CGNR_R_HAT   (1U<< 8)

typedef struct infft_plan_ 
{
  nfft_plan *mv;		/**< matrix vector multiplication     */
  unsigned flags;			/**< iteration type, damping, weights */
  
  double *w;                         	/**< weighting factors                */
  double *w_hat;                     	/**< damping factors                  */
  double *r_hat;                        /**< regularisation weights           */
  double lambda;                        /**< regularisation parameter         */

  /* right hand side */
  complex *y;                /**< right hand side, sampled values  */

  /* solution */ 
  complex *f_hat_iter;             /**< iterative solution               */

  /* vectors */
  complex *r_iter;                 /**< iterated original residual vector*/ 
  complex *z_hat_iter;             /**< residual vector of normal eq.    */
  complex *p_hat_iter;             /**< (non damped) search direction    */
  complex *v_iter;                 /**< residual vector update in CGNR   */

  /* factors */
  double alpha_iter;                    /**< step size for search direction   */
  double beta_iter;                     /**< step size for search correction  */

  /* dot products */
  double dot_r_iter;                    /**< dotproductc{_w}(r_iter)          */
  double dot_r_iter_old;                /**< old dotproductc{_w}(r_iter)      */
  double dot_z_hat_iter;                /**< dotproductc{_w}(z_hat_iter)      */
  double dot_z_hat_iter_old;            /**< old dotproductc{_w}(z_hat_iter)  */
  double dot_p_hat_iter;                /**< dotproductc{_w}(p_hat_iter)      */
  double dot_v_iter;                    /**< dotproductc{_w}(v_iter)          */  
} infft_plan;

/** initialisation for inverse transform, simple interface
 */
void infft_init(infft_plan *ths, nfft_plan *mv);

/** initialisation for inverse transform, specific interface
 */
void infft_init_specific(infft_plan *ths, nfft_plan *mv, int infft_flags);

/** computes start residuals
 */
void infft_before_loop(infft_plan *ths);

/** iterates one step
 */
void infft_loop_one_step(infft_plan *ths);

/** finalisation for inverse transform
 */
void infft_finalize(infft_plan *ths);
/** @} 
 */ 

#endif
