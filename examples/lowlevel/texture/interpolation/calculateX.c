#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <nfft3.h>
#include <util.h>

#define MAX_LINE 100

int grid_type;
const char *grid_type_descr[] = { "equidistant", "file" };
char grid_h_file[100];
char grid_r_file[100];

int N, h_theta_count, h_phi_count, r_theta_count, r_phi_count;
int N1, N2;
complex *omega, *x;
double *h_phi, *h_theta, *r;

void init()
{
}

void read_properties(const char *propfile_name)
{
	FILE *f = fopen(propfile_name, "r");

	fscanf(f, "%d", &grid_type);
	printf("Samples\n");
	printf("# grid_type: %d (%s)\n", grid_type, grid_type_descr[grid_type]);

	if (grid_type == 0) {
		scanf("%d%d%d%d", &h_phi_count, &h_theta_count, &r_phi_count,
					&r_theta_count);
		printf("# h_phi_count: %d h_theta_count: %d\n", h_phi_count,
					 h_theta_count);
		printf("# r_phi_count: %d r_theta_count: %d\n", r_phi_count,
					 r_theta_count);
	} else {
		fscanf(f, "%s%s", grid_h_file, grid_r_file);
	}
	printf("#\n");
	fclose(f);
}

void read_omega(const char *file)
{
	int i;
	FILE *f = fopen(file, "r");
	char line[MAX_LINE];
	char blank;

	fgets(line, MAX_LINE, f);
	if (strcmp(line, "Omega\n")) {
		fprintf(stderr, "Invalid omega file!\n");
		fflush(0);
		exit(-1);
	}

	printf("# From the omega file:\n");

	fgets(line, MAX_LINE, f);
	while (line[0] == '#') {
		printf("%s", line);
		fgets(line, MAX_LINE, f);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in omega file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(f, "%d", &N);

	printf("# N=%d\n", N);
	printf("#\n");

	omega = (complex *) malloc(texture_flat_length(N) * sizeof(complex));
	for (i = 0; i < texture_flat_length(N); i++) {
		double r, im;

		if (!fscanf(f, "%lg + %lgi", &r, &im)) {
			fprintf(stderr, "Error parsing input!\n");
			fflush(0);
			exit(-1);
		} else {
			omega[i] = r + I * im;
		}
	}

	do {
		blank = fgetc(f);
	} while (isspace(blank));

	if (!feof(f)) {
		fprintf(stderr, "Parse error in omega file!\n");
		fflush(0);
		exit(-1);
	}

	fclose(f);
}

void read_grid_h(const char *file)
{
	int i;
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

	fscanf(f, "%d", &N1);

	printf("# N1=%d\n", N1);
	printf("#\n");

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

	fscanf(f, "%d", &N2);

	printf("# N2=%d\n", N2);
	printf("#\n");

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

double equidist(double start, double end, int i, int n, int incl)
{
	if (!incl) {
		n++;
	}
	return start + (end - start) * (double) i / (double) (n - 1);
}

void create_grid()
{
	if (grid_type == 0) {
		int i, j;
		int s, t;

		N1 = h_phi_count * h_theta_count;
		N2 = r_phi_count * r_theta_count;

		h_phi = (double *) malloc(N1 * sizeof(double));
		h_theta = (double *) malloc(N1 * sizeof(double));
		r = (double *) malloc(N1 * N2 * 2 * sizeof(double));

		for (s = 0; s < h_theta_count; s++) {
			double theta = equidist(0, PI, s, h_theta_count, 1);
			for (t = 0; t < h_phi_count; t++) {
				double phi = equidist(-PI, PI, t, h_phi_count, 0);
				h_phi[s * h_phi_count + t] = phi;
				h_theta[s * h_phi_count + t] = theta;
			}
		}

		j = 0;
		for (s = 0; s < r_theta_count; s++) {
			double theta = equidist(0, PI, s, r_theta_count, 1);
			for (t = 0; t < r_phi_count; t++) {
				double phi = equidist(-PI, PI, t, r_phi_count, 0);
				r[j++] = phi;
				r[j++] = theta;
			}
		}

		for (i = 1; i < N1; i++) {
			memcpy(&(r[2 * i * N2]), r, 2 * N2 * sizeof(double));
		}
	} else {
		read_grid_h(grid_h_file);
		read_grid_r(grid_r_file);
	}
}

void calculate_x()
{
	texture_plan plan;

	x = (complex *) malloc(N1 * N2 * sizeof(complex));

	texture_precompute(N);
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);
	texture_trafo(&plan);
	texture_finalize(&plan);

	texture_forget();
}

void output_x()
{
	int i, j;

	printf("\n");
	printf("%d %d\n", N1, N2);
	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			printf("%lg + %lgi\n", creal(x[i * N2 + j]), cimag(x[i * N2 + j]));
		}
	}
}

void cleanup()
{
	free(omega);
	free(x);
	free(h_phi);
	free(h_theta);
	free(r);
}

void usage()
{
	// TODO
}

int main(int argc, char *argv[])
{
	const char *omega_file = "omega";
	const char *propfile_name = "propfile_x";

	if (argc > 1) {
		omega_file = argv[1];
	}

	if (argc <= 2) {
		init();
		read_properties(propfile_name);
		read_omega(omega_file);
		create_grid();
		fflush(0);
		calculate_x();
		output_x();
		cleanup();
	} else {
		usage();
	}

	return 0;
}
