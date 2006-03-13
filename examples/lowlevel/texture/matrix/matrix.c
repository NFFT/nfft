#include<stdio.h>
#include<complex.h>
#include<stdlib.h>

#include<nfft3.h>
#include<texture_util.h>

int N;
int h_phi_count, h_theta_count, r_phi_count, r_theta_count;
int adjoint;
int N1, N2;

void read_input()
{
	fprintf(stderr, "N: ");
	scanf("%d", &N);
	fprintf(stderr, "h_phi, h_theta, r_phi, r_theta: ");
	scanf("%d%d%d%d", &h_phi_count, &h_theta_count, &r_phi_count,
				&r_theta_count);
	fprintf(stderr, "Output adjoint matrix (0 = no, 1 = yes)? ");
	scanf("%d", &adjoint);

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
	fprintf(stderr, "N1*N2=%d\n", N1*N2);
	fprintf(stderr, "#J_N=%d\n", texture_flat_length(N));
	fflush(0);
}

int main()
{
	complex *omega, *x;
	texture_plan plan;

	read_input();

	make_grid();

	output_info();

	omega = (complex *) smart_calloc(texture_flat_length(N), sizeof(complex));
	x = (complex *) smart_calloc(N1 * N2, sizeof(complex));

	texture_precompute(N);

	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	{
		int r;

		if (adjoint) {
			printf("%d %d\n", texture_flat_length(N), N1 * N2);

			for (r = 0; r < texture_flat_length(N); r++) {
				int s;

				omega[r] = 1;
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
