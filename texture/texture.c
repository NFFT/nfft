#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include <util.h>

#include "texture.h"

void texture_precompute(int N) {
	texture_precompute_advanced(N, TEXTURE_DEF_PRECOMPUTE_FLAGS, 
			TEXTURE_DEF_NFSFT_PRECOMPUTE_FLAGS, TEXTURE_DEF_NFSFT_THRESHOLD);
}

void texture_precompute_advanced(int N, unsigned int texture_precompute_flags, 
		unsigned int nfsft_precompute_flags, double nfsft_threshold) {

	if (is_precomputed && 
		 	(N > precompute_params.N ||
		 	 texture_precompute_flags != precompute_params.texture_precompute_flags ||
			 nfsft_precompute_flags != precompute_params.nfsft_precompute_flags ||
			 abs(nfsft_threshold - precompute_params.nfsft_threshold) > 1E-15)) {

 		nfsft_forget_old();
	}		
	
	nfsft_precompute_old(N, nfsft_threshold, nfsft_precompute_flags);
	
	precompute_params.N = N;
	precompute_params.texture_precompute_flags = texture_precompute_flags;
	precompute_params.nfsft_precompute_flags = nfsft_precompute_flags;
	precompute_params.nfsft_threshold = nfsft_threshold;

	is_precomputed = 1;
}

void texture_init(texture_plan *ths, int N,	int N1, int N2, complex* omega, 
		complex* x, const double* h_phi, const double* h_theta, const double* r) {
	
	texture_init_advanced(ths, N, N1, N2, omega, x, h_phi, h_theta, r, 
			TEXTURE_DEF_INIT_FLAGS, TEXTURE_DEF_NFSFT_INIT_FLAGS, 
			TEXTURE_DEF_NFFT_CUTOFF);
}

void texture_init_advanced(texture_plan *ths, int N, int N1, int N2, 
		complex* omega, complex* x, const double* h_phi, const double* h_theta, 
		const double* r, unsigned int texture_init_flags, 
		unsigned int nfsft_init_flags, int nfft_cutoff) {

	int n;
	int N_tilde = next_power_of_2(N);

	ths->N = N;	
	ths->N1 = N1;
	ths->N2 = N2;
	ths->N_total = texture_get_omega_length(ths);
	ths->M_total = texture_get_x_length(ths);

	ths->nfsft_f_hat = (complex**) malloc(sizeof(complex*) * (2*N + 1));
	for(n = -N; n <= N; n++)
	{
		ths->nfsft_f_hat[n+N] = (complex*) malloc(sizeof(complex) * (N_tilde + 1));
	}
	ths->nfsft_f = (complex*) malloc(sizeof(complex) * ths->M_total);

	ths->nfsft_angles = (double*) malloc(sizeof(double) * N1 * N2 * 2);
	ths->cos_h_theta = (double*) malloc(sizeof(double) * N1);
	ths->sin_h_theta = (double*) malloc(sizeof(double) * N1);

	texture_set_nfsft_init_flags(ths, nfsft_init_flags);
	texture_set_nfft_cutoff(ths, nfft_cutoff);

	texture_set_omega(ths, omega);
	texture_set_x(ths, x);
	
	texture_set_h_phi(ths, h_phi);
	texture_set_h_theta(ths, h_theta);
 	texture_set_r(ths, r);	
}

/**
 * Carries out the transform.
 */
void texture_trafo(texture_plan *ths) {
	int i;

	for(i = 0; i < ths->N1; i++) {
		nfsft_plan_old nfsft;
		nfsft = prepare_nfsft_plan(ths, i);
		
		nfsft_trafo_old(nfsft);
		
		nfsft_finalize_old(nfsft);
	}
}

/**
 * The adjoint version of the transform.
 */
void texture_adjoint(texture_plan *ths) {
	int i;

	memcpy(ths->nfsft_f, ths->f, sizeof(complex) * ths->M_total);
	memset(ths->f_hat, 0, sizeof(complex) * ths->N_total);
	
	for(i = 0; i < ths->N1; i++) {
		nfsft_plan_old nfsft;

		nfsft = nfsft_init_guru_old(ths->N, ths->N2, ths->nfsft_f_hat, 
				&(ths->nfsft_angles[2 * ths->N2 * i]), &(ths->nfsft_f[ths->N2 * i]), 
				ths->nfsft_init_flags, ths->nfft_cutoff);
		nfsft_adjoint_old(nfsft);
		nfsft_finalize_old(nfsft);

		process_results(ths, i);
	}
}

/**
 * Frees all memory allocated by texture_init.
 */
void texture_finalize(texture_plan *ths) {
	int n;
	
	free(ths->cos_h_theta);
	free(ths->sin_h_theta);

	for(n = -ths->N; n <= ths->N; n++)
	{
		free(ths->nfsft_f_hat[n+ths->N]);
	}
	free(ths->nfsft_f_hat);
	free(ths->nfsft_f);
	free(ths->nfsft_angles);
}

void texture_forget() {
	nfsft_forget_old();
	is_precomputed = 0;
}

/**
 * Convert a non-flat index to a flat index.
 * \arg l the frequence
 * \arg m ranges from -l to l
 * \arg n ranges from -l to l
 */
inline int texture_flat_index(int l, int m, int n) {
	return m + n + ((l*(5 + 6*l + 4*l*l + 6*m)) / 3);
}

inline int texture_flat_length(int N) {
	return texture_flat_index(N, N, N) + 1;
}

inline int texture_get_omega_length(texture_plan *ths) {
	return texture_flat_length(ths->N);
}

inline int texture_get_x_length(texture_plan *ths) {
	return ths->N1 * ths->N2;
}

inline int texture_get_N(texture_plan *ths) {
	return ths->N;
}

inline int texture_get_N1(texture_plan *ths) {
	return ths->N1;
}

inline int texture_get_N2(texture_plan *ths) {
	return ths->N2;
}

inline const complex *texture_get_omega(texture_plan *ths) {
	return ths->f_hat;
}

inline void texture_set_omega(texture_plan *ths, complex *omega) {
	ths->f_hat = omega;
}

inline const complex *texture_get_x(texture_plan *ths) {
	return ths->f;
}

inline void texture_set_x(texture_plan *ths, complex* x) {
	ths->f = x;
}

inline const double *texture_get_h_phi(texture_plan *ths) {
	return ths->h_phi;
}

inline void texture_set_h_phi(texture_plan *ths, const double* h_phi) {
	ths->h_phi = h_phi;
}

inline const double *texture_get_h_theta(texture_plan *ths) {
	return ths->h_theta;
}

inline void texture_set_h_theta(texture_plan *ths, const double* h_theta) {
	int i;

	ths->h_theta = h_theta;

	for(i = 0; i < ths->N1; i++) {
		ths->cos_h_theta[i] = cos(ths->h_theta[i]);
		ths->sin_h_theta[i] = sin(ths->h_theta[i]);
	}
}

inline const double *texture_get_r(texture_plan *ths) {
	return ths->r;
}

inline void texture_set_r(texture_plan *ths, const double* r) {
	int i;
	
	ths->r = r;

	for(i = 0; i < ths->N1 * ths->N2; i++) {
		ths->nfsft_angles[2*i] = -ths->r[2*i] / TEXTURE_MAX_ANGLE;
		ths->nfsft_angles[2*i + 1] = ths->r[2*i + 1] / TEXTURE_MAX_ANGLE;
	}
}

inline unsigned int texture_get_nfsft_init_flags(texture_plan *ths) {
	return ths->nfsft_init_flags;
}

inline void texture_set_nfsft_init_flags(texture_plan *ths, 
		unsigned int nfsft_init_flags) {
	ths->nfsft_init_flags = nfsft_init_flags;
}

inline int texture_get_nfft_cutoff(texture_plan *ths) {
	return ths->nfft_cutoff;
}

inline void texture_set_nfft_cutoff(texture_plan *ths, int nfft_cutoff) {
	ths->nfft_cutoff = nfft_cutoff;
}

static nfsft_plan_old prepare_nfsft_plan(texture_plan *ths, int i) {
	double p_diag = 1;
	int n;

	for(n = -ths->N; n <= ths->N; n++)
	{
		memset(&(ths->nfsft_f_hat[n + ths->N][abs(n)]), 0, 
				(ths->N - abs(n) + 1) * sizeof(complex));
	}
	
	{
		double p_old2 = 0;
		double p_old1 = p_diag;
		int l;

		ths->nfsft_f_hat[0 + ths->N][0] += ths->f_hat[texture_flat_index(0, 0, 0)];

		for(l = 1; l <= ths->N; l++) {
			int m;
			double p = 
				(double) (2*l - 1) / (double) (l) *	ths->cos_h_theta[i] * p_old1 
				-	(double) (l - 1) / (double) (l)	*	p_old2;

			p_old2 = p_old1;
			p_old1 = p;

			for(m = -l; m <= l; m++) {
				ths->nfsft_f_hat[m + ths->N][l] 
					+= p * ths->f_hat[texture_flat_index(l, m, 0)]; 
			}
		}
	}
	
	for(n = 1; n <= ths->N; n++) {
		double p_old2 = 0;
		double p_old1;
		int l, m;
		
		p_diag *= sqrt((double) (2*n - 1) / (double) (2*n)) 
			* ths->sin_h_theta[i];

		p_old1 = p_diag;
		
		for(m = -n; m <= n; m++) {
			ths->nfsft_f_hat[m + ths->N][n] 
				+= p_diag 
					* (ths->f_hat[texture_flat_index(n, m, n)] 
							* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i])) 
						+ ths->f_hat[texture_flat_index(n, m, -n)] 
							*	(cos(n * ths->h_phi[i]) + I * sin(n * ths->h_phi[i]))); 
		}
		
		for(l = n+1; l <= ths->N; l++) {
			double p = 
				(double) (2*l - 1) / sqrt((double) (l*l - n*n)) 
					* ths->cos_h_theta[i] * p_old1 
				-	sqrt((double) ((l - n - 1) * (l + n - 1)) / (double) (l*l - n*n))
					*	p_old2;
		
			p_old2 = p_old1;
			p_old1 = p;

			for(m = -l; m <= l; m++) {
				ths->nfsft_f_hat[m + ths->N][l] += p 
					* (ths->f_hat[texture_flat_index(l, m, n)] 
								* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]))
							+ ths->f_hat[texture_flat_index(l, m, -n)]
								* (cos(n * ths->h_phi[i]) + I * sin(n * ths->h_phi[i])));
			}
		}
	}

	return nfsft_init_guru_old(ths->N, ths->N2, ths->nfsft_f_hat, 
			&(ths->nfsft_angles[2 * ths->N2 * i]), &(ths->f[ths->N2 * i]), 
			ths->nfsft_init_flags, ths->nfft_cutoff);
}

static void process_results(texture_plan *ths, int i) {
	double p_diag = 1;
	double p_old2 = 0;
	double p_old1 = p_diag;
	int l, n;

	ths->f_hat[texture_flat_index(0, 0, 0)] += ths->nfsft_f_hat[0 + ths->N][0];

	for(l = 1; l <= ths->N; l++) {
		int m;
		double p = 
			(double) (2*l - 1) / (double) (l) * ths->cos_h_theta[i] * p_old1 
			- (double) (l - 1) / (double) (l) * p_old2;
		
		p_old2 = p_old1;
		p_old1 = p;

		for(m = -l; m <= l; m++) {
			ths->f_hat[texture_flat_index(l, m, 0)] 
				+= p * ths->nfsft_f_hat[m + ths->N][l];
		}
	}

	for(n = 1; n <= ths->N; n++) {
		int m;
		p_diag *= sqrt((double) (2*n - 1) / (double) (2*n)) * ths->sin_h_theta[i];
		p_old2 = 0;
		p_old1 = p_diag;

		for(m = -n; m <= n; m++) {
			ths->f_hat[texture_flat_index(n, m, n)] 
				+= p_diag * ths->nfsft_f_hat[m + ths->N][n] 
					* (cos(n * ths->h_phi[i]) + I * sin(n * ths->h_phi[i]));
			
			ths->f_hat[texture_flat_index(n, m, -n)] 
				+= p_diag * ths->nfsft_f_hat[m + ths->N][n]
					* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]));
		}

		for(l = n + 1; l <= ths->N; l++) {
			double p =
				(double) (2*l - 1) / sqrt((double) (l*l - n*n)) 
					* ths->cos_h_theta[i] * p_old1
				- sqrt((double) ((l - n - 1) * (l + n - 1)) / (double) (l*l - n*n))
					* p_old2;

			p_old2 = p_old1;
			p_old1 = p;

			for(m = -l; m <= l; m++) {
				ths->f_hat[texture_flat_index(l, m, n)] 
					+= p * ths->nfsft_f_hat[m + ths->N][l] 
						* (cos(n * ths->h_phi[i]) + I * sin(n * ths->h_phi[i]));
				ths->f_hat[texture_flat_index(l, m, -n)] 
					+= p * ths->nfsft_f_hat[m + ths->N][l]
						* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]));
			}
		}
	}
}
