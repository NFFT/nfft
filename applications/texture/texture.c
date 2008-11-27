#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#include "util.h"

#include "texture.h"
#include "solver.c"

/** @file texture.c
 * Provides the implementation of all methods of the texture transform.
 * @ingroup texture_private
 */

void texture_precompute(int N)
{
	texture_precompute_advanced(N, TEXTURE_DEF_PRECOMPUTE_FLAGS,
															TEXTURE_DEF_NFSFT_PRECOMPUTE_FLAGS,
															TEXTURE_DEF_FPT_PRECOMPUTE_FLAGS,
															TEXTURE_DEF_NFSFT_THRESHOLD);
}

void texture_precompute_advanced(int N, unsigned int texture_precompute_flags,
																 unsigned int nfsft_precompute_flags,
																 unsigned int fpt_precompute_flags,
																 double nfsft_threshold)
{

	// Forget nfsft wisdom if necessary.
	if (is_precomputed &&
			(N > precompute_params.N ||
			 texture_precompute_flags != precompute_params.texture_precompute_flags
			 || nfsft_precompute_flags != precompute_params.nfsft_precompute_flags
			 || abs(nfsft_threshold - precompute_params.nfsft_threshold) > 1E-15)) {

		nfsft_forget();
	}
	// Precompute nfsft wisdom.
	nfsft_precompute(N, nfsft_threshold, nfsft_precompute_flags,
									 fpt_precompute_flags);

	// Store precompute state.
	precompute_params.N = N;
	precompute_params.texture_precompute_flags = texture_precompute_flags;
	precompute_params.nfsft_precompute_flags = nfsft_precompute_flags;
	precompute_params.nfsft_threshold = nfsft_threshold;

	is_precomputed = 1;
}

void texture_init(texture_plan * ths, int N, int N1, int N2,
									double _Complex * omega, double _Complex * x,
									const double *h_phi, const double *h_theta, const double *r)
{

	texture_init_advanced(ths, N, N1, N2, omega, x, h_phi, h_theta, r,
												TEXTURE_DEF_INIT_FLAGS, TEXTURE_DEF_NFSFT_INIT_FLAGS,
												TEXTURE_DEF_NFFT_INIT_FLAGS, TEXTURE_DEF_NFFT_CUTOFF);
}

void texture_init_advanced(texture_plan * ths, int N, int N1, int N2,
													 double _Complex * omega, double _Complex * x,
													 const double *h_phi, const double *h_theta,
													 const double *r, unsigned int texture_init_flags,
													 unsigned int nfsft_init_flags,
													 unsigned int nfft_init_flags, int nfft_cutoff)
{
	// Set sizes.
	ths->N = N;
	ths->N1 = N1;
	ths->N2 = N2;
	ths->N_total = texture_flat_length(N);
	ths->M_total = N1 * N2;

	// Allocate storage for internal computations.
	ths->nfsft_f_hat =nfft_malloc(sizeof(double _Complex) * NFSFT_F_HAT_SIZE(N));
	ths->nfsft_f =
		(double _Complex *)nfft_malloc(sizeof(double _Complex) * ths->M_total);

	ths->nfsft_angles = (double *)nfft_malloc(sizeof(double) * N1 * N2 * 2);
	ths->cos_h_theta = (double *)nfft_malloc(sizeof(double) * N1);
	ths->sin_h_theta = (double *)nfft_malloc(sizeof(double) * N1);

	// Set Parameters.
	texture_set_nfsft_init_flags(ths, nfsft_init_flags);
	texture_set_nfft_init_flags(ths, nfft_init_flags);
	texture_set_nfft_cutoff(ths, nfft_cutoff);

	// Set input / output.
	texture_set_omega(ths, omega);
	texture_set_x(ths, x);

	// Set nodes.
	texture_set_h_phi(ths, h_phi);
	texture_set_h_theta(ths, h_theta);
	texture_set_r(ths, r);
}

/**
 * Carries out the transform.
 */
void texture_trafo(texture_plan * ths)
{
	int i;
	// new
	// texture_plan p2;
	// prepare_nfsft_plan(&p2, 0, i);



	for (i = 0; i < ths->N1; i++) {
		nfsft_plan nfsft;
		// printf("%d\n", i);
		prepare_nfsft_plan(ths, &nfsft, i);
		// prepare_nfsft_plan(&p2, 0, i);

		nfsft_precompute_x(&nfsft);
		nfsft_trafo(&nfsft);

		nfsft_finalize(&nfsft);
	}
}

/**
 * The adjoint version of the transform.
 */
void texture_adjoint(texture_plan * ths)
{
	int i;

	memcpy(ths->nfsft_f, ths->f, sizeof(double _Complex) * ths->M_total);
	memset(ths->f_hat, 0, sizeof(double _Complex) * ths->N_total);

	for (i = 0; i < ths->N1; i++) {
		nfsft_plan nfsft;

		nfsft_init_guru(&nfsft, ths->N, ths->N2, ths->nfsft_init_flags,
										ths->nfft_init_flags, ths->nfft_cutoff);
		nfsft.f_hat = ths->nfsft_f_hat;
		nfsft.x = &(ths->nfsft_angles[2 * ths->N2 * i]);
		nfsft.f = &(ths->nfsft_f[ths->N2 * i]);

		nfsft_precompute_x(&nfsft);
		nfsft_adjoint(&nfsft);

		process_results(ths, &nfsft, i);

		nfsft_finalize(&nfsft);
	}
}

/**
 * Frees all memory allocated by texture_init.
 */
void texture_finalize(texture_plan * ths)
{
	free(ths->cos_h_theta);
	free(ths->sin_h_theta);

	free(ths->nfsft_f_hat);
	free(ths->nfsft_f);
	free(ths->nfsft_angles);
}

void texture_forget()
{
	nfsft_forget();
	is_precomputed = 0;
}

inline int texture_flat_index(int l, int m, int n)
{
	return m + n + ((l * (5 + 6 * l + 4 * l * l + 6 * m)) / 3);
}

inline int texture_flat_length(int N)
{
	return texture_flat_index(N, N, N) + 1;
}

inline int texture_get_omega_length(texture_plan * ths)
{
	return ths->N_total;
}

inline int texture_get_x_length(texture_plan * ths)
{
	return ths->M_total;
}

inline int texture_get_N(texture_plan * ths)
{
	return ths->N;
}

inline int texture_get_N1(texture_plan * ths)
{
	return ths->N1;
}

inline int texture_get_N2(texture_plan * ths)
{
	return ths->N2;
}

inline double _Complex *texture_get_omega(texture_plan * ths)
{
	return ths->f_hat;
}

inline void texture_set_omega(texture_plan * ths, double _Complex * omega)
{
	ths->f_hat = omega;
}

inline double _Complex *texture_get_x(texture_plan * ths)
{
	return ths->f;
}

inline void texture_set_x(texture_plan * ths, double _Complex * x)
{
	ths->f = x;
}

inline const double *texture_get_h_phi(texture_plan * ths)
{
	return ths->h_phi;
}

inline void texture_set_h_phi(texture_plan * ths, const double *h_phi)
{
	ths->h_phi = h_phi;
}

inline const double *texture_get_h_theta(texture_plan * ths)
{
	return ths->h_theta;
}

inline void texture_set_h_theta(texture_plan * ths, const double *h_theta)
{
	int i;

	ths->h_theta = h_theta;

	for (i = 0; i < ths->N1; i++) {
		ths->cos_h_theta[i] = cos(ths->h_theta[i]);
		ths->sin_h_theta[i] = sin(ths->h_theta[i]);
	}
}

inline const double *texture_get_r(texture_plan * ths)
{
	return ths->r;
}

inline void texture_set_r(texture_plan * ths, const double *r)
{
	int i;

	ths->r = r;

	for (i = 0; i < ths->N1 * ths->N2; i++) {
		ths->nfsft_angles[2 * i] = ths->r[2 * i] / TEXTURE_MAX_ANGLE;
		ths->nfsft_angles[2 * i + 1] = ths->r[2 * i + 1] / TEXTURE_MAX_ANGLE;
	}
}

inline unsigned int texture_get_nfsft_init_flags(texture_plan * ths)
{
	return ths->nfsft_init_flags;
}

inline void texture_set_nfsft_init_flags(texture_plan * ths,
																				 unsigned int nfsft_init_flags)
{
	ths->nfsft_init_flags = nfsft_init_flags;
}

inline void texture_set_nfft_init_flags(texture_plan * ths,
																				unsigned int nfft_init_flags)
{
	ths->nfft_init_flags = nfft_init_flags;
}

inline int texture_get_nfft_cutoff(texture_plan * ths)
{
	return ths->nfft_cutoff;
}

inline void texture_set_nfft_cutoff(texture_plan * ths, int nfft_cutoff)
{
	ths->nfft_cutoff = nfft_cutoff;
}

void prepare_nfsft_plan(texture_plan * ths, nfsft_plan * nfsft, int i)
{
	// Stores the associated legendre polynomial of order and degree |n|.
	// p_diag = p_n^n(cos_h_theta[i])
	double p_diag = 1;
	int n;

	// Initialise the frequencies for the nfsft_plan with 0.
	memset(ths->nfsft_f_hat, 0,
				 NFSFT_F_HAT_SIZE(ths->N) * sizeof(double _Complex));

	// Build the nfsft_plan.
	nfsft_init_guru(nfsft, ths->N, ths->N2, ths->nfsft_init_flags,
									ths->nfft_init_flags, ths->nfft_cutoff);

	// Sum up the frequencies for n = 0.
	{
		// Stores the values with lower index in the three term recurrency
		// relation for the associated legendre polynomials.
		double p_old2 = 0;
		double p_old1 = p_diag;
		int l;

		// Add the frequency with m = l = 0.
		ths->nfsft_f_hat[NFSFT_INDEX(0, 0, nfsft)] +=
			ths->f_hat[texture_flat_index(0, 0, 0)];

		// Add the frequencies with 0 < l <= N.
		for (l = 1; l <= ths->N; l++) {
			int m;
			// Stores the associated legendre polynomial of order 0 and degree l.
			// p = p_l^0(cos_h_theta[i])
			double p =
				(double) (2 * l - 1) / (double) (l) * ths->cos_h_theta[i] * p_old1
				- (double) (l - 1) / (double) (l) * p_old2;

			p_old2 = p_old1;
			p_old1 = p;

			// Add the frequencies with -l <= m <= l.
			for (m = -l; m <= l; m++) {
				ths->nfsft_f_hat[NFSFT_INDEX(l, m, nfsft)]
					+= p * ths->f_hat[texture_flat_index(l, m, 0)];
			}
		}
	}

	// Sum up the frequencies for 0 < |n| <= N.
	for (n = 1; n <= ths->N; n++) {
		// Stores the values with lower index in the three term recurrency
		// relation for the associated legendre polynomials.
		double p_old2 = 0;
		double p_old1;
		int l, m;

		p_diag *= sqrt((double) (2 * n - 1) / (double) (2 * n))
			* ths->sin_h_theta[i];

		p_old1 = p_diag;

		// Add the frequencies for l = |n|, -l <= m <= l.
		for (m = -n; m <= n; m++) {
			ths->nfsft_f_hat[NFSFT_INDEX(n, m, nfsft)]
				+= p_diag * (ths->f_hat[texture_flat_index(n, m, n)]
										 * (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]))
										 + ths->f_hat[texture_flat_index(n, m, -n)]
										 * (cos(n * ths->h_phi[i]) + _Complex_I * sin(n * ths->h_phi[i])));
		}

		// Add the frequencies for |n| < l <= N.
		for (l = n + 1; l <= ths->N; l++) {
			// Stores the associated legendre polynomial of order |n| and degree l.
			// p = p_l^{|n|}(cos_h_theta[i])
			double p = (double) (2 * l - 1) / sqrt((double) (l * l - n * n))
				* ths->cos_h_theta[i] * p_old1
				-
				sqrt((double) ((l - n - 1) * (l + n - 1)) / (double) (l * l - n * n))
				* p_old2;

			p_old2 = p_old1;
			p_old1 = p;

			// Add the frequencies for -l <= m <= l.
			for (m = -l; m <= l; m++) {
				ths->nfsft_f_hat[NFSFT_INDEX(l, m, nfsft)] += p
					* (ths->f_hat[texture_flat_index(l, m, n)]
						 * (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]))
						 + ths->f_hat[texture_flat_index(l, m, -n)]
						 * (cos(n * ths->h_phi[i]) + _Complex_I * sin(n * ths->h_phi[i])));
			}
		}
	}

	// Put the values in the nfsft_plan.
	nfsft->f_hat = ths->nfsft_f_hat;
	nfsft->x = &(ths->nfsft_angles[2 * ths->N2 * i]);
	nfsft->f = &(ths->f[ths->N2 * i]);
}

void process_results(texture_plan * ths, nfsft_plan * nfsft, int i)
{
	// Stores the associated legendre polynomial of order and degree |n|.
	// p_diag = p_n^n(cos_h_theta[i])
	double p_diag = 1;
	// Stores the values with lower index in the three term recurrency
	// relation for the associated legendre polynomials.
	double p_old2 = 0;
	double p_old1 = p_diag;
	int l, n;

	// Add the frequency with n = m = l = 0.
	ths->f_hat[texture_flat_index(0, 0, 0)] +=
		ths->nfsft_f_hat[NFSFT_INDEX(0, 0, nfsft)];

	// Add the frequencies for n = 0, 0 < l <= N.
	for (l = 1; l <= ths->N; l++) {
		int m;
		// Stores the associated legendre polynomial of order 0 and degree l.
		// p = p_l^0(cos_h_theta[i])
		double p =
			(double) (2 * l - 1) / (double) (l) * ths->cos_h_theta[i] * p_old1
			- (double) (l - 1) / (double) (l) * p_old2;

		p_old2 = p_old1;
		p_old1 = p;

		// Add the frequencies for -l <= m <= l.
		for (m = -l; m <= l; m++) {
			ths->f_hat[texture_flat_index(l, m, 0)]
				+= p * ths->nfsft_f_hat[NFSFT_INDEX(l, m, nfsft)];
		}
	}

	// Add the frequencies for 0 < |n| <= N.
	for (n = 1; n <= ths->N; n++) {
		int m;
		p_diag *=
			sqrt((double) (2 * n - 1) / (double) (2 * n)) * ths->sin_h_theta[i];
		p_old2 = 0;
		p_old1 = p_diag;

		// Add the frequencies for l = |n|, -l <= m <= l.
		for (m = -n; m <= n; m++) {
			ths->f_hat[texture_flat_index(n, m, n)]
				+= p_diag * ths->nfsft_f_hat[NFSFT_INDEX(n, m, nfsft)]
				* (cos(n * ths->h_phi[i]) + _Complex_I * sin(n * ths->h_phi[i]));

			ths->f_hat[texture_flat_index(n, m, -n)]
				+= p_diag * ths->nfsft_f_hat[NFSFT_INDEX(n, m, nfsft)]
				* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]));
		}

		// Add the freqeuncies for |n| < l <= N.
		for (l = n + 1; l <= ths->N; l++) {
			// Stores the associated legendre polynomial of order |n| and degree l.
			// p = p_l^{|n|}(cos_h_theta[i])
			double p = (double) (2 * l - 1) / sqrt((double) (l * l - n * n))
				* ths->cos_h_theta[i] * p_old1
				-
				sqrt((double) ((l - n - 1) * (l + n - 1)) / (double) (l * l - n * n))
				* p_old2;

			p_old2 = p_old1;
			p_old1 = p;

			// Add the frequencies for -l <= m <= l.
			for (m = -l; m <= l; m++) {
				ths->f_hat[texture_flat_index(l, m, n)]
					+= p * ths->nfsft_f_hat[NFSFT_INDEX(l, m, nfsft)]
					* (cos(n * ths->h_phi[i]) + _Complex_I * sin(n * ths->h_phi[i]));
				ths->f_hat[texture_flat_index(l, m, -n)]
					+= p * ths->nfsft_f_hat[NFSFT_INDEX(l, m, nfsft)]
					* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]));
			}
		}
	}
}

MACRO_SOLVER_IMPL(texture, complex, double _Complex)
