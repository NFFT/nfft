#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>

#include<nfft3.h>
#include"texture_util.h"

const char *omega_policy_descr[] = { "flat", "1/n" };

const char *grid_descr[] = { "equidistant angles", "file" };

static void check_eof(FILE * f, const char *location)
{
	char blank;

	do {
		blank = fgetc(f);
	} while (isspace(blank));

	if (!feof(f)) {
		fprintf(stderr, "Parse error at the end of %s!\n", location);
		fflush(0);
		exit(-1);
	}
}

void read_h(int *N1_ptr, double **h_phi_ptr, double **h_theta_ptr, FILE * in,
						FILE * out)
{
	int i;
	int N1;
	double *h_phi, *h_theta;
	char line[MAX_LINE];

	fgets(line, MAX_LINE, in);
	if (strcmp(line, "Polefigures\n")) {
		fprintf(stderr, "Invalid pole figure file!\n");
		fflush(0);
		exit(-1);
	}

	fprintf(out, "# From the pole figure file:\n");

	fgets(line, MAX_LINE, in);
	while (line[0] == '#') {
		fprintf(out, "%s", line);
		fgets(line, MAX_LINE, in);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in pole figure file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(in, "%d", N1_ptr);
	N1 = *N1_ptr;

	fprintf(out, "# N1=%d\n", N1);
	fprintf(out, "#\n");

	*h_phi_ptr = (double *) malloc(N1 * sizeof(double));
	*h_theta_ptr = (double *) malloc(N1 * sizeof(double));
	h_phi = *h_phi_ptr;
	h_theta = *h_theta_ptr;

	for (i = 0; i < N1; i++) {
		if (!fscanf(in, "%lg%lg", &(h_phi[i]), &(h_theta[i]))) {
			fprintf(stderr, "Parse error in pole figure file!\n");
			fflush(0);
			exit(-1);
		}
	}

	check_eof(in, "pole figure file");
}

void read_r(int N1, int *N2_ptr, double **r_ptr, FILE * in, FILE * out)
{
	int i, j;
	char line[MAX_LINE];
	int N2;
	double *r;

	fgets(line, MAX_LINE, in);
	if (strcmp(line, "Nodes\n")) {
		fprintf(stderr, "Invalid node file!\n");
		fflush(0);
		exit(-1);
	}

	fprintf(out, "# From the node file:\n");

	fgets(line, MAX_LINE, in);
	while (line[0] == '#') {
		fprintf(out, "%s", line);
		fgets(line, MAX_LINE, in);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in node file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(in, "%d", N2_ptr);
	N2 = *N2_ptr;

	fprintf(out, "# N2=%d\n", N2);
	fprintf(out, "#\n");

	*r_ptr = (double *) malloc(N1 * N2 * 2 * sizeof(double));
	r = *r_ptr;

	for (j = 0; j < 2 * N2; j++) {
		if (!fscanf(in, "%lg", &(r[j]))) {
			fprintf(stderr, "Parse error in node file!\n");
			fflush(0);
			exit(-1);
		}
	}

	for (i = 1; i < N1; i++) {
		memcpy(&(r[2 * N2 * i]), r, 2 * N2 * sizeof(double));
	}

	check_eof(in, "node file");
}

void read_grid(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
							 double **h_theta_ptr, double **r_ptr, FILE * h_in, FILE * r_in,
							 FILE * out)
{
	read_h(N1_ptr, h_phi_ptr, h_theta_ptr, h_in, out);
	read_r(*N1_ptr, N2_ptr, r_ptr, r_in, out);
}

void read_x(int *N1_ptr, int *N2_ptr, complex ** x_ptr, FILE * in, FILE * out)
{
	int i, j;
	int N1, N2;
	complex *x;
	char line[MAX_LINE];

	fgets(line, MAX_LINE, in);
	if (strcmp(line, "Samples\n")) {
		fprintf(stderr, "Invalid sample file!\n");
		fflush(0);
		exit(-1);
	}

	fprintf(out, "# From the sample file:\n");

	fgets(line, MAX_LINE, in);
	while (line[0] == '#') {
		fprintf(out, "%s", line);
		fgets(line, MAX_LINE, in);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in sample file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(in, "%d%d", N1_ptr, N2_ptr);
	N1 = *N1_ptr;
	N2 = *N2_ptr;

	fprintf(out, "# N1=%d N2=%d\n", N1, N2);
	fprintf(out, "#\n");

	*x_ptr = (complex *) malloc(N1 * N2 * sizeof(complex));
	x = *x_ptr;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			double r, im = 0;
			if (!fscanf(in, "%lg + %lgi", &r, &im)) {
				fprintf(stderr, "Parse error in sample file!\n");
				fflush(0);
				exit(-1);
			} else {
				x[i * N2 + j] = r + I * im;
			}
		}
	}

	check_eof(in, "sample file");
}

void read_samples(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
									double **h_theta_ptr, double **r_ptr, complex ** x_ptr,
									FILE * h_in, FILE * r_in, FILE * x_in, FILE * out)
{
	int N1, N2;

	read_grid(N1_ptr, N2_ptr, h_phi_ptr, h_theta_ptr, r_ptr, h_in, r_in, out);
	N1 = *N1_ptr;
	N2 = *N2_ptr;

	read_x(N1_ptr, N2_ptr, x_ptr, x_in, out);
	if (N1 != *N1_ptr || N2 != *N2_ptr) {
		fprintf(stderr, "Conflicting N1 = %d | %d, N2 = %d | %d\n", N1, *N1_ptr,
						N2, *N2_ptr);
		fflush(0);
		exit(-1);
	}
}

void read_omega(int *N_ptr, complex ** omega_ptr, FILE * in, FILE * out)
{
	int i;
	int N;
	complex *omega;
	char line[MAX_LINE];

	fgets(line, MAX_LINE, in);
	if (strcmp(line, "Omega\n")) {
		fprintf(stderr, "Invalid omega file!\n");
		fflush(0);
		exit(-1);
	}

	fprintf(out, "# From the omega file:\n");

	fgets(line, MAX_LINE, in);
	while (line[0] == '#') {
		fprintf(out, "%s", line);
		fgets(line, MAX_LINE, in);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in omega file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(in, "%d", N_ptr);
	N = *N_ptr;

	fprintf(out, "# N=%d\n", N);
	fprintf(out, "#\n");

	*omega_ptr = (complex *) malloc(texture_flat_length(N) * sizeof(complex));
	omega = *omega_ptr;

	for (i = 0; i < texture_flat_length(N); i++) {
		double r, im = 0;

		if (!fscanf(in, "%lg + %lgi", &r, &im)) {
			fprintf(stderr, "Error parsing input!\n");
			fflush(0);
			exit(-1);
		} else {
			omega[i] = r + I * im;
		}
	}

	check_eof(in, "omega file");
}

void write_r(int N2, double *r, FILE * out)
{
	int j;

	fprintf(out, "\n");
	fprintf(out, "%d\n", N2);

	for (j = 0; j < N2; j++) {
		fprintf(out, "%lg %lg\n", r[2 * j], r[2 * j + 1]);
	}
}

void write_h(int N1, double *h_phi, double *h_theta, FILE * out)
{
	int i;

	fprintf(out, "\n");
	fprintf(out, "%d\n", N1);

	for (i = 0; i < N1; i++) {
		fprintf(out, "%lg %lg\n", h_phi[i], h_theta[i]);
	}
}

void write_omega(int N, complex * omega, FILE * out)
{
	int i;

	fprintf(out, "\n");
	fprintf(out, "%d\n", N);
	for (i = 0; i < texture_flat_length(N); i++) {
		fprintf(out, "%lg + %lgi\n", creal(omega[i]), cimag(omega[i]));
	}
}

void write_x(int N1, int N2, complex * x, FILE * out)
{
	int i, j;

	fprintf(out, "\n");
	fprintf(out, "%d %d\n", N1, N2);
	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			fprintf(out, "%lg + %lgi\n", creal(x[i * N2 + j]),
							cimag(x[i * N2 + j]));
		}
	}
}

double equidist(double start, double end, int i, int n, int incl)
{
	if (!incl) {
		n++;
	}
	return start + (end - start) * (double) i / (double) (n - 1);
}

void calculate_grid(int h_phi_count, int h_theta_count, int r_phi_count,
										int r_theta_count, double *h_phi, double *h_theta,
										double *r, int type)
{
	int N1, N2;

	N1 = h_phi_count * h_theta_count;
	N2 = r_phi_count * r_theta_count;

	switch (type) {
		case 0:
		{
			int i, j;
			int s, t;

			for (s = 0; s < h_theta_count; s++) {
				double theta =
					equidist(0, TEXTURE_MAX_ANGLE / 2, s, h_theta_count, 1);
				for (t = 0; t < h_phi_count; t++) {
					double phi =
						equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, t,
										 h_phi_count, 0);
					h_phi[s * h_phi_count + t] = phi;
					h_theta[s * h_phi_count + t] = theta;
				}
			}

			j = 0;
			for (s = 0; s < r_theta_count; s++) {
				double theta =
					equidist(0, TEXTURE_MAX_ANGLE / 2, s, r_theta_count, 1);
				for (t = 0; t < r_phi_count; t++) {
					double phi =
						equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, t,
										 r_phi_count, 0);
					r[j++] = phi;
					r[j++] = theta;
				}
			}

			for (i = 1; i < N1; i++) {
				memcpy(&(r[2 * i * N2]), r, 2 * N2 * sizeof(double));
			}
			break;
		}
		case 1:
			error("Cannot calculate a grid from type file!\n");
		default:
		{
			fprintf(stderr, "Illegal grid type (%d)!", type);
			fflush(0);
			exit(-1);
		}
	}
}

void error(const char *msg)
{
	fprintf(stderr, "%s", msg);
	fflush(0);
	exit(-1);
}

void init_omega(complex * omega, int N, int policy)
{
	int i;
	switch (policy) {
		case 0:
		{
			for (i = 0; i < texture_flat_length(N); i++) {
				omega[i] = drand48() * rand() + I * drand48() * rand();
			}
			break;
		}
		case 1:
		{
			for (i = 0; i < texture_flat_length(N); i++) {
				omega[i] = (drand48() * rand() + I * drand48() * rand()) / (i + 1);
			}
			break;
		}
		default:
			error("illegal omega policy\n");
	}
}

void mult_error(int N1, int N2, complex *x, double min_err, double max_err)
{
	int i, j;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			x[i * N2 + j] *= min_err + drand48() * (max_err - min_err);
		}
	}
}

void block_mult_error(int N1, int N2, complex *x, double min_err,
											double max_err)
{
	int i, j;

	for (i = 0; i < N1; i++) {
		double err = min_err + drand48() * (max_err - min_err);
		for (j = 0; j < N2; j++) {
			x[i * N2 + j] *= err;
		}
	}
}
