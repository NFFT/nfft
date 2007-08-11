#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include<math.h>

#include<nfft3_texture.h>
#include<texture_util.h>

/**
 * @defgroup texture_matrix Texture: Matrix
 * This program calculates the matrix of the texture transform.
 *
 * @section Inp Input
 * The requested input is explained by messages on stderr.
 * h_phi, h_theta, r_phi and r_theta are the number of different
 * @f$\phi@f$s and @f$\theta@f$s of the pole figures (h_phi and h_theta) and
 * nodes (r_phi and r_theta).
 * 
 * @section ProcOutp Processing and Output
 * The program prints some message to stdout cointaining the raw dimensions 
 * of the matrix.
 *
 * If the user choosed to calculate the adjoint matrix, the adjoint matrix is
 * calculated using @ref texture_trafo and printed to stdout.
 * If weight factors where set row @f$r@f$ is divided through 
 * @f$(r+1)^{weight/2}@f$.
 *
 * If the user choosed to calculate the matrix, it is calculated using 
 * @ref texture_adjoint and printed to stdout.
 * If weight factors where set column @f$r@f$ is divided through 
 * @f$(r+1)^{weight/2}@f$.
 * 
 * @author Matthias Schmalz
 * @ingroup texture_examples
 */

int N;
int h_phi_count, h_theta_count, r_phi_count, r_theta_count;
int adjoint;
int N1, N2;
double weight;

void read_input()
{
	fprintf(stderr, "N: ");
	scanf("%d", &N);
	fprintf(stderr, "h_phi, h_theta, r_phi, r_theta: ");
	scanf("%d%d%d%d", &h_phi_count, &h_theta_count, &r_phi_count,
				&r_theta_count);
	fprintf(stderr, "Output adjoint matrix (0 = no, 1 = yes)? ");
	scanf("%d", &adjoint);
	fprintf(stderr, "weight factors (input i -> 1/n^i): ");
	scanf("%lg", &weight);

	N1 = h_phi_count * (h_theta_count - 2) + 2;
	N2 = r_phi_count * (r_theta_count - 2) + 2;
}

double *h_phi, *h_theta, *r;

void make_grid()
{
	grid_dim dim;
	dim.angles.h_phi_count = h_phi_count;
	dim.angles.h_theta_count = h_theta_count;
	dim.angles.r_phi_count = r_phi_count;
	dim.angles.r_theta_count = r_theta_count;

	h_phi = smart_malloc(N1 * sizeof(double));
	h_theta = smart_malloc(N1 * sizeof(double));
	r = smart_malloc(N1 * N2 * 2 * sizeof(double));

	calculate_grid(dim, h_phi, h_theta, r, 0);
}

void output_info()
{
	fprintf(stderr, "N1=%d, N2=%d\n", N1, N2);
	fprintf(stderr, "N1*N2=%d\n", N1 * N2);
	fprintf(stderr, "#J_N=%d\n", texture_flat_length(N));
	fflush(0);
}

int main()
{
	double _Complex *omega, *x;
	texture_plan plan;

	read_input();

	make_grid();

	output_info();

	omega = (double _Complex *) smart_calloc(texture_flat_length(N), sizeof(double _Complex));
	x = (double _Complex *) smart_calloc(N1 * N2, sizeof(double _Complex));

	texture_precompute(N);

	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	{
		int r;

		if (adjoint) {
			printf("%d %d\n", texture_flat_length(N), N1 * N2);

			for (r = 0; r < texture_flat_length(N); r++) {
				int s;

				omega[r] = ((double) 1) / (pow(r + 1, weight / 2));
				texture_set_omega(&plan, omega);

				texture_trafo(&plan);

				for (s = 0; s < N1 * N2; s++) {
					printf("%lg + %lgi ", creal(x[s]), -cimag(x[s]));
				}
				printf("\n");

				omega[r] = 0;
			}
		} else {
			printf("%d %d\n", N1 * N2, texture_flat_length(N));

			for (r = 0; r < N1 * N2; r++) {
				int s;

				x[r] = 1;
				texture_set_x(&plan, x);

				texture_adjoint(&plan);

				for (s = 0; s < texture_flat_length(N); s++) {
					omega[s] /= (pow(s + 1, weight / 2));
					printf("%lg + %lgi ", creal(omega[s]), -cimag(omega[s]));
				}
				printf("\n");

				x[r] = 0;
			}
		}
	}

	texture_finalize(&plan);

	texture_forget();

	free(omega);
	free(x);

	free(h_phi);
	free(h_theta);
	free(r);

	return 0;
}
