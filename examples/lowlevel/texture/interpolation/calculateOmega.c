#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <nfft3.h>

#define MAX_LINE 100

typedef struct properties_ {
	double res_delta;
	int iterations_without_check;
	double min_improve;
} properties;

properties prop;

double *h_phi, *h_theta, *r;
complex *x, *omega;
int N, N1, N2;

void init()
{
}

void read_N()
{
	scanf("%d", &N);
	printf("Omega\n");
}

void read_properties(const char *propfile)
{
	FILE *f = fopen(propfile, "r");

	fscanf(f, "%lg%d%lg", &prop.res_delta, &prop.iterations_without_check,
				 &prop.min_improve);

	fclose(f);

	printf("# res_delta: %lg\n", prop.res_delta);
	printf("# iterations_without_check: %d\n", prop.iterations_without_check);
	printf("# min_improve: %lg\n", prop.min_improve);
	printf("#\n");
}

void read_samples(const char *file)
{
	int i, j;
	FILE *f = fopen(file, "r");
	char line[MAX_LINE];
	char blank;

	fgets(line, MAX_LINE, f);
	if (strcmp(line, "Samples\n")) {
		fprintf(stderr, "Invalid sample file!\n");
		fflush(0);
		exit(-1);
	}

	printf("# From the sample file:\n");

	fgets(line, MAX_LINE, f);
	while (line[0] == '#') {
		printf("%s", line);
		fgets(line, MAX_LINE, f);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in sample file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(f, "%d%d", &N1, &N2);

	printf("# N1=%d N2=%d\n", N1, N2);
	printf("#\n");

	x = (complex *) malloc(N1 * N2 * sizeof(complex));
	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			double r, im;
			if (!fscanf(f, "%lg", &r)) {
				fprintf(stderr, "Parse error in sample file!\n");
				fflush(0);
				exit(-1);
			} else {
				fscanf(f, " + %lgi", &im);
				x[i * N2 + j] = r + I * im;
			}
		}
	}

	do {
		blank = fgetc(f);
	} while (isspace(blank));

	if (!feof(f)) {
		fprintf(stderr, "Parse error at the end of sample file!\n");
		fflush(0);
		exit(-1);
	}

	fclose(f);
}

void read_grid_h(const char *file)
{
	int i;
	int N1_new;
	FILE *f = fopen(file, "r");
	char line[MAX_LINE];
	char blank;

	fgets(line, MAX_LINE, f);
	if (strcmp(line, "Polefigures\n")) {
		fprintf(stderr, "Invalid pole figure file!\n");
		fflush(0);
		exit(-1);
	}

	printf("# From the pole figure file:\n");

	fgets(line, MAX_LINE, f);
	while (line[0] == '#') {
		printf("%s", line);
		fgets(line, MAX_LINE, f);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in pole figure file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(f, "%d", &N1_new);

	printf("# N1=%d\n", N1_new);
	printf("#\n");

	if (N1 != N1_new) {
		fprintf(stderr, "Inconsistent N1 in pole figure file!\n");
		fflush(0);
		exit(-1);
	}

	h_phi = (double *) malloc(N1 * sizeof(double));
	h_theta = (double *) malloc(N1 * sizeof(double));

	for (i = 0; i < N1; i++) {
		if (!fscanf(f, "%lg%lg", &(h_phi[i]), &(h_theta[i]))) {
			fprintf(stderr, "Parse error in pole figure file!\n");
			fflush(0);
			exit(-1);
		}
	}

	do {
		blank = fgetc(f);
	} while (isspace(blank));

	if (!feof(f)) {
		fprintf(stderr, "Parse error in pole figure file!\n");
		fflush(0);
		exit(-1);
	}

	fclose(f);
}

void read_grid_r(const char *file)
{
	int i, j;
	int N2_new;
	FILE *f = fopen(file, "r");
	char line[MAX_LINE];
	char blank;

	fgets(line, MAX_LINE, f);
	if (strcmp(line, "Nodes\n")) {
		fprintf(stderr, "Invalid node file!\n");
		fflush(0);
		exit(-1);
	}

	printf("# From the node file:\n");

	fgets(line, MAX_LINE, f);
	while (line[0] == '#') {
		printf("%s", line);
		fgets(line, MAX_LINE, f);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in node file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(f, "%d", &N2_new);

	printf("# N2=%d\n", N2_new);
	printf("#\n");

	if (N2 != N2_new) {
		fprintf(stderr, "Inconsistent N2 in node file!\n");
		fflush(0);
		exit(-1);
	}

	r = (double *) malloc(N1 * N2 * 2 * sizeof(double));

	for (j = 0; j < 2 * N2; j++) {
		if (!fscanf(f, "%lg", &(r[j]))) {
			fprintf(stderr, "Parse error in node file!\n");
			fflush(0);
			exit(-1);
		}
	}

	for (i = 1; i < N1; i++) {
		memcpy(&(r[2 * N2 * i]), r, 2 * N2 * sizeof(double));
	}

	do {
		blank = fgetc(f);
	} while (isspace(blank));

	if (!feof(f)) {
		fprintf(stderr, "Parse error in node file!\n");
		fflush(0);
		exit(-1);
	}

	fclose(f);
}

inline double l_2_norm(const complex * vec, unsigned int length)
{
	double sum = 0.0;
	int i;

	for (i = 0; i < length; i++) {
		double x = cabs(vec[i]);
		sum += x * x;
	}

	return sqrt(sum);
}

inline double l_2_dist(const complex * x, const complex * y,
											 unsigned int length)
{
	double sum = 0.0;
	int i;

	for (i = 0; i < length; i++) {
		double diff = cabs(x[i] - y[i]);
		sum += diff * diff;
	}

	return sqrt(sum);
}

double res_dist(const complex * vec, const complex * ref, int length)
{
	double abs_dist = l_2_dist(vec, ref, length);
	double base = l_2_norm(ref, length);
	if (base != 0.0 || abs_dist != 0.0) {
		return abs_dist / base;
	} else {
		return 0.0;
	}
}

void calculate_omega()
{
	itexture_plan iplan;
	texture_plan plan, test_plan;

	omega = (complex *) malloc(texture_flat_length(N) * sizeof(complex));

	texture_precompute(N);
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);
	texture_init(&test_plan, N, N1, N2, omega, x, h_phi, h_theta, r);
	itexture_init(&iplan, &plan);
	memcpy(iplan.y, x, N1 * N2 * sizeof(complex));
	memset(iplan.f_hat_iter, 0, texture_flat_length(N) * sizeof(complex));
	itexture_before_loop(&iplan);
	{
		double new_res, old_res;

		texture_set_omega(&test_plan, iplan.f_hat_iter);
		texture_trafo(&test_plan);
		new_res =
			res_dist(texture_get_x(&test_plan), iplan.y,
							 texture_get_x_length(&test_plan));

		do {
			int count;

			for (count = 0; count < prop.iterations_without_check; count++) {
				itexture_loop_one_step(&iplan);
			}

			old_res = new_res;
			texture_set_omega(&test_plan, iplan.f_hat_iter);
			texture_trafo(&test_plan);
			new_res =
				res_dist(texture_get_x(&test_plan), iplan.y,
								 texture_get_x_length(&test_plan));
#ifdef DEBUG_RESIDUUM
			fprintf(stderr, "residuum: %lg\n", new_res);
			fflush(0);
#endif
		} while (old_res > 0 && new_res > prop.res_delta
						 && (old_res - new_res) / old_res >= prop.min_improve);

#ifdef DEBUG_CONVERGENCE_QUALITY
		if (new_res > prop.res_delta) {
			fprintf(stderr, "No convergence, residuum was %lg!\n", new_res);
			fflush(0);
		}
#endif
		printf("# residuum: %lg\n", new_res);
		printf("#\n");

		memcpy(omega, iplan.f_hat_iter, texture_flat_length(N) * sizeof(complex));
		itexture_finalize(&iplan);
		texture_finalize(&plan);
		texture_finalize(&test_plan);
	}

	texture_forget();
}

void output_omega()
{
	int i;

	printf("\n");
	printf("%d\n", N);
	for (i = 0; i < texture_flat_length(N); i++) {
		printf("%lg + %lgi\n", creal(omega[i]), cimag(omega[i]));
	}
}

void cleanup()
{
	free(h_phi);
	free(h_theta);
	free(r);
	free(x);
	free(omega);
}

void usage()
{
}

int main(int argc, char *argv[])
{
	const char *sample_file = "samples.in";
	const char *grid_h_file = "grid_h.in";
	const char *grid_r_file = "grid_r.in";

	if (argc > 1) {
		sample_file = argv[1];
	}
	if (argc > 2) {
		grid_h_file = argv[2];
	}
	if (argc > 3) {
		grid_r_file = argv[3];
	}

	if (argc <= 4) {

		init();
		read_N();
		read_properties("propfile_omega");
		read_samples(sample_file);
		read_grid_h(grid_h_file);
		read_grid_r(grid_r_file);
		fflush(0);
		calculate_omega();
		output_omega();
		cleanup();
	} else {
		usage();
	}

	return 0;
}
