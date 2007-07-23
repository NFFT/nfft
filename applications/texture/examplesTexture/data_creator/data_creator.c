#include<stdlib.h>
#include<stdio.h>
#include<complex.h>

#include<nfft3.h>
#include<texture_util.h>

/** @defgroup data_creator Texture: Data Creator
 * This program was designed to create some samples based on random frequencies.
 * 
 * @section CLA Command Line Arguments
 * -# The name of the pole figure file. (default: grid_h)
 * -# The name of the node file. (default: grid_r)
 * -# The name of the omega file. (default: omega)
 * -# The name of the property file. (default: propfile)
 *
 * @section PF Property File
 * -# Grid type.
 * -# Omega policy
 * -# minimal error
 * -# maximal error
 * -# minimal alpha
 * -# maximal alpha
 *
 * @section Inp Input
 * If the grid type is 0 (equidistant angles), the program reads
 * -# The bandwidth N.
 * -# The number of different phis and thetas of the pole figures.
 * -# The number of different phis and thetas of the nodes of one pole figure.
 *
 * If the grid type is 1 (read from file), the program reads only the bandwidth
 * N.
 *
 * @section ProcOutp Processing and Output
 * First the program prints a brief description of the input to stderr.
 * 
 * If the grid type is 0 (equidistant angles), the program creates a grid and
 * frequencies omega using the according functions of @ref texture_utility.
 * The grid and omega will be written to the files defined by the command line.
 * 
 * The intensities x are calculated and afterwards randomly modified:
 * Each component is multiplied by a random number between the minimal and 
 * maximal error.
 * Each pole figure is multiplied by a random number between the minimal and
 * maximal alpha.
 * Afterwards x is printed to stdout.
 * 
 * If the grid type is 1 (read from file), the program reades the grid from
 * the files defined by the command line and does not write them.
 * 
 * @author Matthias Schmalz
 * @ingroup texture_examples
 */

complex *omega, *x;

double *h_phi, *h_theta, *r;

int N, N1, N2;

int h_phi_count, h_theta_count, r_phi_count, r_theta_count;

int grid;
int omega_policy;
double min_error, max_error;
double min_alpha, max_alpha;

void init()
{
	// unsigned short int seed[] = { 1, 2, 3 };

	// seed48(seed);
	srand(0);
}

void read_propfile(const char *file)
{
	FILE *f = fopen(file, "r");

	fscanf(f, "%d%d", &grid, &omega_policy);
	fprintf(stderr, "grid: %d (%s)\n", grid, grid_descr[grid]);
	fprintf(stderr, "omega_policy: %d (%s)\n", omega_policy,
					omega_policy_descr[omega_policy]);

	fscanf(f, "%lg%lg%lg%lg", &min_error, &max_error, &min_alpha, &max_alpha);
	fprintf(stderr, "rounding error in [%lg, %lg) * 100%%\n", min_error,
					max_error);
	fprintf(stderr, "alpha in [%lg, %lg)\n", min_alpha, max_alpha);
	fclose(f);
}

void read_params()
{
	switch (grid) {
		case 0:
		{
			scanf("%d%d%d%d%d", &N, &h_phi_count, &h_theta_count, &r_phi_count,
						&r_theta_count);
			break;
		}
		case 1:
		{
			scanf("%d", &N);
			break;
		}
		default:
			error("illegal grid\n");
	}
}

void create_grid(const char *h_file, const char *r_file)
{
	switch (grid) {
		case 0:
		{
			FILE *hf = fopen(h_file, "w");
			FILE *rf = fopen(r_file, "w");
			grid_dim dims;

			dims.angles.h_phi_count = h_phi_count;
			dims.angles.h_theta_count = h_theta_count;
			dims.angles.r_phi_count = r_phi_count;
			dims.angles.r_theta_count = r_theta_count;

			N1 = h_phi_count * (h_theta_count - 2) + 2;
			N2 = r_phi_count * (r_theta_count - 2) + 2;

			h_phi = (double *) smart_malloc(N1 * sizeof(double));
			h_theta = (double *) smart_malloc(N1 * sizeof(double));
			r = (double *) smart_malloc(N1 * N2 * 2 * sizeof(double));

			calculate_grid(dims, h_phi, h_theta, r, grid);

			fprintf(hf, "Polefigures\n");
			fprintf(hf, "# This grid has been created by the data creator.\n");
			fprintf(hf, "# h_phi_count = %d, h_theta_count = %d\n", h_phi_count,
							h_theta_count);
			write_h(N1, h_phi, h_theta, hf);
			fprintf(rf, "Nodes\n");
			fprintf(rf, "# This grid has been created by the data creator.\n");
			fprintf(rf, "# r_phi_count = %d, r_theta_count = %d\n", r_phi_count,
							r_theta_count);
			write_r(N2, r, rf);

			fclose(hf);
			fclose(rf);
			break;
		}
		case 1:
		{
			FILE *hf = fopen(h_file, "r");
			FILE *rf = fopen(r_file, "r");
			read_grid(&N1, &N2, &h_phi, &h_theta, &r, hf, rf, stdout);
			fclose(hf);
			fclose(rf);
			break;
		}
		default:
			error("illegal grid\n");
	}
}

void create_omega(const char *file)
{
	FILE *f = fopen(file, "w");

	omega = (complex *) smart_malloc(texture_flat_length(N) * sizeof(complex));
	init_omega(omega, N, omega_policy);
	fprintf(f, "Omega\n");
	fprintf(f, "# This omega has been created by data creator.\n");
	write_omega(N, omega, f);

	fclose(f);
}

void create_x()
{
	texture_plan plan;

	x = (complex *) smart_malloc(N1 * N2 * sizeof(complex));

	texture_precompute(N);
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);
	texture_trafo(&plan);
	texture_finalize(&plan);
	texture_forget();

	mult_error(N1, N2, x, min_error, max_error);
	block_mult_error(N1, N2, x, min_alpha, max_alpha);

	write_x(N1, N2, x, stdout);
}

void cleanup()
{
	free(h_phi);
	free(h_theta);
	free(r);
	free(omega);
	free(x);
}

void usage()
{
	fprintf(stderr, "Wrong parameters!\n");
}

int main(int argc, char *argv[])
{
	const char *propfile_name = "propfile";
	const char *grid_h_def = "grid_h";
	const char *grid_r_def = "grid_r";
	const char *omega_def = "omega";

	if (argc > 1) {
		grid_h_def = argv[1];
	}
	if (argc > 2) {
		grid_r_def = argv[2];
	}
	if (argc > 3) {
		omega_def = argv[3];
	}
	if (argc > 4) {
		propfile_name = argv[4];
	}

	if (argc <= 5) {
		init();
		read_propfile(propfile_name);
		read_params();

		printf("Samples\n");
		printf("# These samples are created by the data creator.\n");

		create_grid(grid_h_def, grid_r_def);
		create_omega(omega_def);
		create_x();
		cleanup();
	} else {
		usage();
	}

	return 0;
}
