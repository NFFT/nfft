#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include <util.h>

#include "texture.h"

void texture_init(texture_plan *ths, int N,	int N1, int N2) {
	texture_init_advanced(ths, N, N1, N2);
}

void texture_init_advanced(texture_plan *ths, int N, int N1, int N2) {
	int n;
	int N_tilde = next_power_of_2(N);

	ths->N = N;
	ths->N1 = N1;
	ths->N2 = N2;
	
	ths->N_total = texture_flat_index(N, N, N) + 1;
	ths->f_hat = (complex*) malloc(sizeof(complex) * ths->N_total);
	ths->M_total = N1 * N2;
	ths->f = (complex*) malloc(sizeof(complex) * ths->M_total);
	
	ths->h_phi = (double*) malloc(sizeof(double) * N1);
	ths->h_theta = (double*) malloc(sizeof(double) * N1);
	ths->r = (double*) malloc(sizeof(double) * N1 * N2 * 2);

	ths->cos_h_theta = (double*) malloc(sizeof(double) * N1);
	ths->sin_h_theta = (double*) malloc(sizeof(double) * N1);

	ths->nfsft_f_hat = (complex**) malloc(sizeof(complex*) * (2*N + 1));
	for(n = -N; n <= N; n++)
	{
		ths->nfsft_f_hat[n+N] = (complex*) malloc(sizeof(complex) * (N_tilde + 1));
	}
	ths->nfsft_f = (complex*) malloc(sizeof(complex) * ths->M_total);
	ths->nfsft_x = (double*) malloc(sizeof(double) * N1 * N2 * 2);
}

/**
 * Carries out the transform.
 */
#ifdef TEXTURE_DUMMY_TRAFO
void texture_trafo(texture_plan *ths) {
	int i, j;
	
	for(i = 0, j = 0; i <= ths->N_total - 2 && j < ths->M_total; i += 2, j++) {
		ths->f[j] = ths->f_hat[i] + ths->f_hat[i+1];
	}
	
	for(; j < ths->M_total; j++) {
		ths->f[j] = 0;
	}
}
#else
void texture_trafo(texture_plan *ths) {
	int i;
	int m;

	precompute(ths);

	//TODO
	nfsft_precompute_old(ths->N, 1000, 0U);

	memset(ths->f, 1, ths->M_total * sizeof(complex));
	vpr_complex(ths->f, ths->M_total, "f before");
	
	for(i = 0; i < ths->N1; i++) {
		nfsft_plan_old nfsft;
		nfsft = prepare_nfsft_plan(ths, i);
		
		for(m = -ths->N; m <= ths->N; m++) {
			vpr_complex(ths->nfsft_f_hat[m+ths->N], next_power_of_2(ths->N) + 1, "nfsft_f_hat");
		}
		vpr_double(ths->nfsft_x, ths->N1 * ths->N2 * 2, "nfsft_x");
		
		ndsft_trafo_old(nfsft);
		
		vpr_complex(ths->f, ths->M_total, "f");
		
		nfsft_finalize_old(nfsft);
	}
}
#endif

/**
 * The adjoint version of the transform.
 */
#ifdef TEXTURE_DUMMY_TRAFO
void texture_adjoint(texture_plan *ths) {
	int i, j;

	for(i = 0, j = 0; i <= ths->N_total - 2 && j < ths->M_total; i += 2, j++) {
		ths->f_hat[i] = ths->f[j];
		ths->f_hat[i+1] = ths->f[j];
	}

	for(; i < ths->N_total; i++) {
		ths->f_hat[i] = 0;
	}
}
#else
void texture_adjoint(texture_plan *ths) {
	int i;

	precompute(ths);
	texture_vec_init(ths->f_hat, ths->N_total, 0);
	memcpy(ths->nfsft_f, ths->f, sizeof(complex) * ths->M_total);
	
	for(i = 0; i < ths->N1; i++) {
		nfsft_plan_old nfsft;

		nfsft = nfsft_init_old(ths->N, ths->N2, ths->nfsft_f_hat, 
				&(ths->nfsft_x[2 * ths->N2 * i]), &(ths->nfsft_f[ths->N2 * i]), 0U);
		ndsft_adjoint_old(nfsft);
		nfsft_finalize_old(nfsft);

		process_results(ths, i);
	}
}
#endif

/**
 * Frees all memory allocated by texture_init.
 */
void texture_finalize(texture_plan *ths) {
	int n;
	
	free(ths->f_hat);
	free(ths->f);
	free(ths->h_phi);
	free(ths->h_theta);
	free(ths->r);

	free(ths->cos_h_theta);
	free(ths->sin_h_theta);

	for(n = -ths->N; n <= ths->N; n++)
	{
		free(ths->nfsft_f_hat[n+ths->N]);
	}
	free(ths->nfsft_f_hat);
	free(ths->nfsft_f);
	free(ths->nfsft_x);
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

static inline void precompute(texture_plan *ths) {
	int i;

	for(i = 0; i < ths->N1; i++) {
		ths->cos_h_theta[i] = cos(ths->h_theta[i]);
		ths->sin_h_theta[i] = sin(ths->h_theta[i]);
	}
	for(i = 0; i < ths->N1 * ths->N2; i++) {
		ths->nfsft_x[2*i] = -ths->r[2*i] / TEXTURE_MAX_ANGLE;
		ths->nfsft_x[2*i + 1] = ths->r[2*i + 1] / TEXTURE_MAX_ANGLE;
	}
}

static nfsft_plan_old prepare_nfsft_plan(texture_plan *ths, int i) {
	double p_diag = 1;
	int n;

	//TODO optimization
	for(n = 0; n < 2*ths->N + 1; n++)
	{
		texture_vec_init(ths->nfsft_f_hat[n], next_power_of_2(ths->N) + 1, 0);
	}
	
	{
		double p_old2 = 0;
		double p_old1 = p_diag;
		int l;

		ths->nfsft_f_hat[0][0] += ths->f_hat[texture_flat_index(0, 0, 0)];

		for(l = 1; l <= ths->N; l++) {
			int m;
			double p = 
				(double) (2*l - 1) / (double) (l) *	ths->cos_h_theta[i] * p_old1 
				-	(double) (l - 1) / (double) (l)	*	p_old2;

			p_old2 = p_old1;
			p_old1 = p;

			for(m = -l; m <= l; m++) {
				ths->nfsft_f_hat[m+l][l] += p * ths->f_hat[texture_flat_index(l, m, 0)]; 
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
			ths->nfsft_f_hat[n+m][n] += 
				p_diag 
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
				ths->nfsft_f_hat[m+l][l] += p 
					* (ths->f_hat[texture_flat_index(l, m, n)] 
								* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]))
							+ ths->f_hat[texture_flat_index(l, m, -n)]
								* (cos(n * ths->h_phi[i]) + I * sin(n * ths->h_phi[i])));
			}
		}
	}

/*	{
		int test;
		complex** test2;
		test = ths->N;
		test = ths->N2;
		test2 = ths->nfsft_f_hat;
		
	}

	for (n = 0; n < ths->N * 2 + 1; n++) {
		vpr_complex(ths->nfsft_f_hat[n], next_power_of_2(ths->N) + 1, "nfsft_f_hat");
	}
	printf("\n\n\n");*/

	return nfsft_init_old(ths->N, ths->N2, ths->nfsft_f_hat, 
			&(ths->nfsft_x[2 * ths->N2 * i]), &(ths->f[ths->N2 * i]), 0U);
}

static void process_results(texture_plan *ths, int i) {
	double p_diag = 1;
	double p_old2 = 0;
	double p_old1 = p_diag;
	int l, n;

	ths->f_hat[texture_flat_index(0, 0, 0)] = ths->nfsft_f_hat[0][0];

	for(l = 1; l <= ths->N; l++) {
		int m;
		double p = 
			(double) (2*l - 1) / (double) (l) * ths->cos_h_theta[i] * p_old1 
			- (double) (l - 1) / (double) (l) * p_old2;
		
		p_old2 = p_old1;
		p_old1 = p;

		for(m = -l; m <= l; m++) {
			ths->f_hat[texture_flat_index(l, m, 0)] += p * ths->nfsft_f_hat[m+l][l];
		}
	}

	for(n = 1; n <= ths->N; n++) {
		int m;
		p_diag *= sqrt((double) (2*n - 1) / (double) (2*n) * ths->sin_h_theta[i]);
		p_old2 = 0;
		p_old1 = p_diag;

		for(m = -n; m <= n; m++) {
			ths->f_hat[texture_flat_index(n, m, n)] += p_diag * ths->nfsft_f_hat[m+l][l] 
				* (cos(n * ths->h_phi[i]) + I * sin(n * ths->h_phi[i]));
			ths->f_hat[texture_flat_index(n, m, -n)] += p_diag * ths->nfsft_f_hat[m+l][l]
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
				ths->f_hat[texture_flat_index(l, m, n)] += p * ths->nfsft_f_hat[m+l][l] 
					* (cos(n * ths->h_phi[i]) + I * sin(n * ths->h_phi[i]));
				ths->f_hat[texture_flat_index(l, m, -n)] += p * ths->nfsft_f_hat[m+l][l]
					* (cos(n * ths->h_phi[i]) - I * sin(n * ths->h_phi[i]));
			}
		}
	}
}

inline void texture_vec_init(complex* vec, int n, complex value) {
	int i;
	
	for(i = 0; i < n; i++)
	{
		vec[i] = value;
	}
}
