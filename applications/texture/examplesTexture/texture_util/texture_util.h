#include<stdio.h>
#include<nfft3_texture.h>

// macros

#define malloc(a) call_smart_malloc_instead(a)
#define calloc(a, b) call_smart_calloc_instead(a)

#define MAX_LINE 1000

// data types

typedef struct sample_dim_ {
	int N1;
	int N2;
} sample_dim;

typedef struct angle_dim_ {
	int h_phi_count;
	int h_theta_count;
	int r_phi_count;
	int r_theta_count;
} angle_dim;

typedef union grid_dim_ {
	sample_dim samples;
	angle_dim angles;
} grid_dim;

// input

void read_h(int *N1_ptr, double **h_phi_ptr, double **h_theta_ptr, FILE * in,
						FILE * out);

void read_r(int N1, int *N2_ptr, double **r_ptr, FILE * in, FILE * out);

void read_grid(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
							 double **h_theta_ptr, double **r_ptr, FILE * h_in, FILE * r_in,
							 FILE * out);

void read_x(int *N1_ptr, int *N2_ptr, complex ** x_ptr, FILE * in,
						FILE * out);

void read_samples(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
									double **h_theta_ptr, double **r_ptr, complex ** x_ptr,
									FILE * h_in, FILE * r_in, FILE * x_in, FILE * out);

void read_omega(int *N_ptr, complex ** omega_ptr, FILE * in, FILE * out);

// output

void write_r(int N2, double *r, FILE * out);

void write_h(int N1, double *h_phi, double *h_theta, FILE * out);

void write_omega(int N, complex * omega, FILE * out);

void write_x(int N1, int N2, complex * x, FILE * out);

// basic helper

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

inline double l_2_rel_dist(const complex * vec, const complex * ref,
													 unsigned int length);

/** Returns @f$ \| vec \|_2.
 *
 * @par length - the length of vec.
 */
inline double l_2_norm(const complex * vec, unsigned int length);

inline double l_2_rel_norm(const complex * vec, const complex * ref,
													 unsigned int length);

inline void *smart_calloc(size_t nmemb, size_t size);

inline void *smart_malloc(size_t size);

/** Returns @f$ Y_k^n(phi, theta)@f$.
 * See @ref sh for the definition of spherical harmonics.
 *
 */
complex spherical_harmonic(int k, int n, double phi, double theta);

inline void split(int n, int *t1, int *t2);

inline void warning(const char *message);

// helper for special purposes

void block_mult_error(int N1, int N2, complex * x, double min_err,
											double max_err);

const char *grid_descr[3];

void calculate_grid(grid_dim dims, double *h_phi, double *h_theta,
										double *r, int grid);

const char *omega_policy_descr[4];

void init_omega(complex * omega, int N, int omega_policy);

void mult_error(int N1, int N2, complex * x, double min_err, double max_err);

const char *solver_algo_descr[2];

unsigned int solver_flags(int solver_algo, int weight_policy);

const char *weight_policy_descr[4];

void set_weights(itexture_plan * iplan, int weight_policy);
