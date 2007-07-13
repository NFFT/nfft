#include<stdio.h>
#include<assert.h>

#include<nfft3_texture.h>

/** @defgroup texture_utility Texture: Utility Functions (advanced)
 * This library provides useful functions for reading, writing, creating and
 * analysing data.
 *
 * @section file_formats File Formats
 * This library provides functions for reading and writing pole figure, node,
 * sample and frequency vectors.
 * We describe the file format.
 * 
 * Lines are seperated by newline characters '\\n'.
 * 
 * The first line contains the key word to identify the type of the data.
 * If the file contains a pole figure vector, the key word is "Polefigures",
 * if, a node vector, it is "Nodes", if, a sample vector, it is "Samples" and,
 * if, a frequency vector, it is "Omega".
 *
 * It follow zero or more lines starting with '#' that contain comments.
 *
 * The next line is empty.
 *
 * The next line contains one or two numbers (separated by one or more blanks),
 * which determine the size of the vector.
 * In the case of a pole figure vector the number is N1, the number of 
 * different pole figures, in the case of a node vector, 
 * it is N2, the number of points on a single pole figure, in the case of a
 * sample vector, it is N1 followed by N2, and in the case of a frequency 
 * vector it is N, the bandwidth.
 *
 * Finally, the elements of the vector are listed.
 * Complex numbers have to be compatible with the format string "%lg + %lgi"
 * of fscanf.
 * If the imaginary part is zero, it can be ommited (including the '+').
 * 
 * @author Matthias Schmalz
 * @ingroup texture_examples
 * @{
 */

// macros

/**
 * Prevents the user from calling malloc instead of ::smart_malloc.
 */
#define malloc(a) call_smart_malloc_instead(a)

/** 
 * Prevents the user from calling calloc instad of ::smart_calloc.
 */
#define calloc(a, b) call_smart_calloc_instead(a)

/**
 * Maximal line length in input files.
 */
#define MAX_LINE 1000

// data types

/**
 * A data structure for the size of samples.
 */
typedef struct sample_dim_ {
	/**
	 * Number of pole figures.
	 */
	int N1;

	/** 
	 * Number of nodes per pole figure.
	 */
	int N2;
} sample_dim;

/**
 * A data structure for the size of a grid with equidistant angles.
 */
typedef struct angle_dim_ {
	/**
	 * The number of different azimuth angles of the pole figure vector.
	 */
	int h_phi_count;

	/**
	 * The number of different elevation angles of the pole figure vector.
	 */
	int h_theta_count;

	/**
	 * The number of different azimuth angles of the node vector.
	 */
	int r_phi_count;

	/**
	 * The number of different elevation angles of the node vector.
	 */
	int r_theta_count;
} angle_dim;

/**
 * A data structur for the size of a grid which can be equidistant but does
 * not have to be.
 */
typedef union grid_dim_ {
	/**
	 * The size, if it is not equidistant.
	 */
	sample_dim samples;

	/**
	 * The size, if it is equidistant.
	 */
	angle_dim angles;
} grid_dim;

/**
 * The parameter set for ::texture_itrafo.
 */
typedef struct itexture_params_ {
	// data input
	/**
	 * The expected result.
	 */
	double complex *omega_ref;

	// input parameters
	/**
	 * The maximum number of epochs.
	 */
	int max_epochs;

	/** 
	 * If after some epoch a file with the specified name appears, the solver 
	 * stops.
	 */
	char stop_file_name[100];

	/**
	 * If the relative @f$l_2@f$ residuum reaches this value, the solver stops.
	 */
	double residuum_goal;

	/**
	 * Let @f$res@f$ be the current relative @f$l_2@f$ residuum and 
	 * @f$res_{min}@f$ the minimal observed @f$l_2@f$ residuum after some epoch.
	 * If @f$res \leq (1 - min\_improve) \cdot res_{min}@f$, we say that
	 * the residuum has improved.
	 */
	double min_improve;

	/**
	 * If the residuum does not improve for the specified number of epochs, the 
	 * solver stops.
	 * Note that the epochs without improvement have to be contiguous.
	 */
	int max_epochs_without_improve;

	/**
	 * If the relative @f$l_2@f$ residuum degrades for the specified 
	 * number of (contiguous) epochs, the solver stops.
	 */
	int max_fail;

	/**
	 * The number of iterations in one epoch.
	 */
	int steps_per_epoch;

	/**
	 * If this is set to true, the solver does not calculate the residuum and 
	 * uses dot_r_iter instead.
	 * However, this could result in undesired behaviour due to numerical errors.
	 */
	int use_updated_residuum;

	/**
	 * If this is set to true, the relative @f$l_2@f$ error between omega_ref 
	 * and f_hat_iter will be monitored.
	 */
	int monitor_error;

	// parameters concerning status messages
	/**
	 * If this is set to true, the solver prints status messages on stderr.
	 */
	int messages_on;

	/**
	 * The interval between status messages (in epochs).
	 */
	int message_interval;

	// data output
	/**
	 * The frequencies at the minimal residuum.
	 */
	double complex *omega_min_res;

	/**
	 * The number of epochs after which the minimum residuum occured.
	 */
	int epochs_until_min_res;

	/**
	 * The relative @f$l_2@f$ value of the minimal residuum.
	 */
	double min_residuum;

	/**
	 * The relative @f$l_2@f$ error at the minimal residuum.
	 */
	double error_during_min_residuum;

	/**
	 * The frequencies at the minimal error.
	 */
	double complex *omega_min_err;

	/**
	 * The relative @f$l_2@f$ value of the minimal error.
	 */
	double min_error;

	/**
	 * The number of epochs after which the minimal error occured.
	 */
	int epochs_until_min_err;

	/**
	 * The status of the itexture_params.
	 * ("initialized", "destroyed",
	 * "residuum goal reached", "degradation", "stagnation", 
	 * "anormal numbers (nan, inf, subnormal)",
	 * "max epochs reached", "stopped")
	 * The solver sets status to the first applying value in the list.
	 */
	const char *status;

} itexture_params;

// input

/**
 * Reads a pole figure vector from in.
 * Allocates storage for the vector.
 * Comments in the file are passed to out.
 * Aborts execution if a syntax error is encountered.
 * @see file_formats
 */
void read_h(int *N1_ptr, double **h_phi_ptr, double **h_theta_ptr, FILE * in,
						FILE * out);

/**
 * Reads a node vector from in.
 * The function reads the nodes for one pole figure, and copies them N1 - 1
 * times such that finally there are nodes for N1 pole figures.
 * Allocates storage for the vector.
 * Comments in the file are passed to out.
 * Aborts execution if a syntax error is encountered.
 * @see file_formats
 */
void read_r(int N1, int *N2_ptr, double **r_ptr, FILE * in, FILE * out);

/**
 * Reads both a pole figure and a node vector from in.
 * Allocates storage for the vectors.
 * Comments in the files are passed to out.
 * Aborts execution if a syntax error is encountered.
 * @see file_formats
 */
void read_grid(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
							 double **h_theta_ptr, double **r_ptr, FILE * h_in, FILE * r_in,
							 FILE * out);

/**
 * Reads a samples vector from in.
 * Allocates storage for the vector.
 * Comments in the file are passed to out.
 * Aborts execution if a syntax error is encountered.
 * @see file_formats
 */
void read_x(int *N1_ptr, int *N2_ptr, complex ** x_ptr, FILE * in,
						FILE * out);

/**
 * Reads a pole figure, a node and a sample vector from in.
 * Allocates storage for the vectors.
 * Comments in the files are passed to out.
 * Aborts execution if a syntax error is encountered, or if the sizes of grid
 * and samples do not fit.
 * @see file_formats
 */
void read_samples(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
									double **h_theta_ptr, double **r_ptr, complex ** x_ptr,
									FILE * h_in, FILE * r_in, FILE * x_in, FILE * out);

/**
 * Reads a frequency vector from in.
 * Allocates storage for the vector.
 * Comments in the file are passed to out.
 * Aborts execution if a syntax error is encountered.
 * @see file_formats
 */
void read_omega(int *N_ptr, complex ** omega_ptr, FILE * in, FILE * out);

/**
 * Reads the bandwidth of a frequency vector from in.
 * Aborts execution if a syntax error is encountered.
 * @see file_formats
 */
void read_N(int *N_ptr, FILE * in);

// output

/**
 * Writes a node vector to a file.
 */
void write_r(int N2, double *r, FILE * out);

/**
 * Writes a pole figure vector to a file.
 */
void write_h(int N1, double *h_phi, double *h_theta, FILE * out);

/**
 * Writes a frequency vector to a file.
 */
void write_omega(int N, complex * omega, FILE * out);

/**
 * Writes a sample vector to a file.
 */
void write_x(int N1, int N2, complex * x, FILE * out);

// basic helper

/**
 * Creates a random number in @f$[0, RAND\_MAX] + i [0, RAND\_MAX]@f$.
 */
inline complex crand();

/**
 * Creates a uniformly distributed random number in the circle with center 
 * @f$0@f$ and radius @f$2^30@f$.
 */
inline complex crand_circle();

/**
 * Creates a random number in @f$[0, RAND\_MAX]@f$.
 */
inline double drand();

/**
 * Creates a random number in @f$[0, 1]@f$.
 */
inline double drand1();

/** Returns if @f$ | x - y | \leq delta @f$.
 */
inline int equal(complex x, complex y, double delta);

/** Returns some value on a equidistant one dimensional grid.
 *
 * @par start - the lower margin
 * @par end - the upper margin
 * @par i - the index of the point that will be returned
 * Ranges from 0 to n-1.
 * @par n - the number of points on the grid
 * @par incl determines if max_value is a point on the grid.
 */
inline double equidist(double start, double end, int i, int n, int incl);

/**
 * Prints msg on stderr and aborts execution.
 */
inline void error(const char *msg);

/** Returns @f$e^{i \cdot phi}@f$
 */
inline complex expi(double phi);

/** Returns @f$ \| x - y \|_2 @f$.
 *
 * @par length - the length of x and y.
 */
inline double l_2_dist(const complex * x, const complex * y,
											 unsigned int length);

/**
 * Returns @f$ \| vec - ref \|_2 / \|ref\|_2 @f$.
 */
inline double l_2_rel_dist(const complex * vec, const complex * ref,
													 unsigned int length);

/** Returns @f$ \| vec \|_2@f$.
 *
 * @par length - the length of vec.
 */
inline double l_2_norm(const complex * vec, unsigned int length);

/**
 * Returns @f$ \| vec \|_2 / \|ref\|_2 @f$.
 */
inline double l_2_rel_norm(const complex * vec, const complex * ref,
													 unsigned int length);

/**
 * Returns @f$ \| vec - ref\|_M / \|ref\|_M @f$, where 
 * @f$ M = diag(m) @f$ with @f$ m(\ell, m, n) = \dfrac{1}{2 \ell + 1}@f$.
 */
inline double compute_ren(const complex * vec, const complex * ref,
													unsigned int N);

/**
 * Returns @f$ \| vec - ref\|_M / \|ref\|_M @f$, where 
 * @f$ M = diag(m) @f$ with @f$ m(\ell, m, n) = \dfrac{1}{(2 \ell + 1)^2}@f$.
 */
inline double compute_rwen(const complex * vec, const complex * ref,
													 unsigned int N);

/**
 * Does the same as calloc, but raises an error if calloc returns 0.
 */
inline void *smart_calloc(size_t nmemb, size_t size);

/**
 * Does the same as malloc, but raises an error if malloc return 0.
 */
inline void *smart_malloc(size_t size);

/** Returns @f$ Y_k^n(phi, theta)@f$.
 * See @ref sh for the definition of spherical harmonics.
 *
 */
complex spherical_harmonic(int k, int n, double phi, double theta);

/**
 * Splits a number @f$n \geq 2 @f$ in two factors t1 and t2 such that t1 and t2
 * have more or less the same size and
 * @f$(t1 \cdot (t2-2)) = n - 2@f$.
 */
inline void split(int n, int *t1, int *t2);

/**
 * Prints a warning message to stdout.
 */
inline void warning(const char *message);

// helper for special purposes

/**
 * Multiplies each pole figure in x which a random number.
 * The random number is fixed for each pole figure but different for
 * different pole figures.
 * It belongs to [min_err, max_err].
 */
void block_mult_error(int N1, int N2, complex * x, double min_err,
											double max_err);

/**
 * Stores the descriptions for different grid types.
 */
const char *grid_descr[3];

/**
 * Creates a grid.
 * If grid = 0 both the pole figure vector and node vector have equidistant
 * angles.
 * That means, the number of points having the same elevation angle is 
 * independent of that elevation angle.
 * There are no duplicated points.
 * grid = 1 is not supported.
 * If grid = 2 the elevation angles are equidistant, and the number of points
 * having the same elevation angle depends on that elevation angle.
 * Near to the poles, their are less such points, near to the equator, more.
 */
void calculate_grid(grid_dim dims, double *h_phi, double *h_theta,
										double *r, int grid);

/**
 * Stores the descriptions for different omega policies.
 */
const char *omega_policy_descr[7];

/**
 * Creates a random frequency vector:
 * omega_policy = 0: Each component is created by ::crand.
 * omega_policy = 1: Each component is created by ::crand/(i+1), where i is its 
 * flat index.
 * omega_policy = 2: Each component is 1. 
 * omega_policy = 3: Each component is created by ::crand/(i+1)^2, where i is 
 * its flat index.
 */
void init_omega(complex * omega, int N, int omega_policy);

/**
 * Multiplies each component of x by a random number in [min_err, max_err].
 */
void mult_error(int N1, int N2, complex * x, double min_err, double max_err);

/**
 * Descriptions for the different types of solver algorithms.
 */
const char *solver_algo_descr[2];

/**
 * Creates the flags for ::itexture_init_advanced corresponding to a
 * weight_policy and a solver algorithm.
 */
unsigned int solver_flags(int solver_algo, int weight_policy);

/**
 * Descriptions for the different weight policies.
 */
const char *weight_policy_descr[6];

/**
 * Initializes the damping factors of iplan.
 * iplan has to be initialized before by itexture_init_advanced with the flags
 * generated by ::solver_flags.
 * weight_policy = 0: Each component is one.
 * weight_policy = 1: Each component is 1/(i+1).
 * weight_policy = 2: Each component is 1/(i+1)^2.
 * weight_policy = 3: Each component is 1/(i+1)^3.
 * In each of the above cases, i is the flat index of the current component.
 */
void set_weights(itexture_plan * iplan, int weight_policy);

// helper for the solver
/**
 * Initializes a parameter set for ::texture_itrafo with default values.
 * Allocates memory for itexture_params::omega_min_res.
 */
void initialize_itexture_params(itexture_params * pars, int N);

/**
 * Requests the user to input the fields of pars that determin the abort 
 * criterion and verbosity of ::texture_itrafo.
 * See the code for details.
 */
void read_itexture_params(itexture_params * pars);

/**
 * Executes ::itexture_before_loop and some iterations of 
 * ::itexture_loop_one_step according to the parameters in pars.
 * @par is the plan passed to the solver functions. Each parameter in par has 
 * to be initialized before invoking texture_itrafo.
 */
void texture_itrafo(itexture_plan * iplan, itexture_params * pars);

/**
 * Frees the memory allocated by ::initialize_itexture_params.
 */
void destroy_itexture_params(itexture_params * pars);
/**
 * @}
 */
