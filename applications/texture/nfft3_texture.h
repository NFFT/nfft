#include<nfft3.h>

#ifndef NFFT3_TEXTURE_H
#define NFFT3_TEXTURE_H

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
	 
/** Default value for fpt_precompute_flags.
 * @see texture_precompute_advanced
 */
#define TEXTURE_DEF_FPT_PRECOMPUTE_FLAGS 0U

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

/** Default value for nfft_init_flags.
 * @see texture_init_advanced
 */
#define TEXTURE_DEF_NFFT_INIT_FLAGS (PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE)

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
  MACRO_MV_PLAN(double complex)

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

  /** The flags for the initialisation of the nfft.
   * @see texture_init_advanced
   */
  unsigned int nfft_init_flags;

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
  double complex *nfsft_f_hat;

  /** Stores the samples for the nfsft transformation.
   */
  double complex *nfsft_f;

  /** Stores the nodes for the nfsft transformation.
   */
  double *nfsft_angles;
} texture_plan;

MACRO_SOLVER_PLAN(texture, complex, double complex)

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
 * @param fpt_precompute_flags - flags for the precomputation of the fpt
 * @param nfsft_threshold - a parameter for the precomputation of the nfsft
 */
void texture_precompute_advanced(int N, unsigned int texture_precompute_flags,
    unsigned int nfsft_precompute_flags, unsigned int fpt_precompute_flags,
		double nfsft_threshold);

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
void texture_init(texture_plan *ths, int N, int N1, int N2, double complex* omega,
    double complex* x, const double* h_phi, const double* h_theta, const double* r);

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
 * referenced to by any of the pointer arguments, the plan can get into an
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
    double complex* omega, double complex* x, const double* h_phi, const double* h_theta,
    const double *r, unsigned int texture_init_flags,
    unsigned int nfsft_init_flags, unsigned int nfft_init_flags, int nfft_cutoff);

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
inline int texture_flat_index(int l, int m, int n);

/** Determines the length of an array omega storing frequencies in a
 * given
 * bandwidth.
 *
 * @par N - the bandwidth.
 * @return the length of the corresponding omega
 */
inline int texture_flat_length(int N);

/** Returns the length of the frequency array stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline int texture_get_omega_length(texture_plan *ths);

/** Returns the length of the sample array stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline int texture_get_x_length(texture_plan *ths);

/** Returns the bandwidth stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline int texture_get_N(texture_plan *ths);

/** Returns the number of pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline int texture_get_N1(texture_plan *ths);

/** Returns the number of samples per pole figure stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline int texture_get_N2(texture_plan *ths);

/** Returns a pointer to the frequencies stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline const double complex *texture_get_omega(texture_plan *ths);

/** Sets the frequencies in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par omega - a pointer to the new frequencies.
 * @pre omega has to point to an array of appropriate length.
 */
inline void texture_set_omega(texture_plan *ths, double complex* omega);

/** Returns a pointer to the samples stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline const double complex *texture_get_x(texture_plan *ths);

/** Sets the samples in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par x - a pointer to the new samples
 * @pre x has to point to an array of appropriate length.
 */
inline void texture_set_x(texture_plan *ths, double complex* x);

/** Returns a pointer to the latitudes of the pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline const double *texture_get_h_phi(texture_plan *ths);

/** Sets the latitudes of the pole figures in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par h_phi - a pointer to the new latitudes
 * @pre h_phi has to point to an array of appropriate length.
 */
inline void texture_set_h_phi(texture_plan *ths, const double* h_phi);

/** Returns the longitudes of the pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline const double *texture_get_h_theta(texture_plan *ths);

/** Sets the longitudes of the pole figures in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par h_theta - a pointer to the longitudes
 * @pre h_theta has to point to an array of appropriate length.
 */
inline void texture_set_h_theta(texture_plan *ths, const double* h_theta);

/** Returns the nodes of the pole figures stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline const double *texture_get_r(texture_plan *ths);

/** Sets the nodes of the pole figures in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par r - a pointer to the nodes
 * @pre r has to point to an array of appropriate length.
 */
inline void texture_set_r(texture_plan *ths, const double* r);

/** Returnes the flags used for the initialisation of the nfsft stored in a
 * plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline unsigned int texture_get_nfsft_init_flags(texture_plan *ths);

/** Sets the flags used for the initialisation of the nfsft in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par nfsft_init_flags - the nfsft flags
 */
inline void texture_set_nfsft_init_flags(texture_plan *ths,
    unsigned int nfsft_init_flags);

/** Sets the flags used for the initialisation of the nfft in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par nfft_init_flags - the nfft flags
 */
inline void texture_set_nfft_init_flags(texture_plan *ths,
    unsigned int nfft_init_flags);

/** Returns the nfft_cutoff parameter used for the initialisation of the nfsft
 * stored in a plan.
 *
 * @par ths - a pointer to the transformation plan
 */
inline int texture_get_nfft_cutoff(texture_plan *ths);

/** Sets the nfft_cutoff parameter used for the initialisation of the nfsft
 * in a plan.
 *
 * @par ths - a pointer to the transformation plan
 * @par nfft_cutoff - the parameter
 */
inline void texture_set_nfft_cutoff(texture_plan *ths, int nfft_cutoff);

/** @}
 */

/** @}
 */
#endif
