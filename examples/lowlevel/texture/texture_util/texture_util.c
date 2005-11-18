#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<math.h>

#include<nfft3.h>
#include<util.h>
#include"texture_util.h"

// macros

// descriptors
const char *grid_descr[] =
	{ "equidistant angles", "file", "uniformly distributed" };

const char *omega_policy_descr[] = { "flat", "1/n" };

const char *solver_algo_descr[] = { "CGNR", "CGNE" };

const char *weight_policy_descr[] = { "flat", "1/n", "1/n^2" };

// internal functions

inline double determine_node_count(int theta_count)
{
	if (theta_count <= 2) {
		return theta_count;
	} else {
		return 2 +
			2 * (theta_count - 1) *
			cos(PI / (theta_count - 1)) *
			sin(theta_count * PI / (2 * (theta_count - 1))) /
			sin(PI / (2 * (theta_count - 1)));
	}
}

inline int determine_theta_count(int node_count)
{
	int t_min = 0, t_max = node_count;

	do {
		int theta_count = (t_min + t_max) / 2;
		double x = determine_node_count(theta_count);

		if (x >= node_count) {
			t_max = theta_count;
		}
		if (x <= node_count) {
			t_min = theta_count;
		}
	} while (t_max - t_min > 1);

	return t_max;
}

inline double determine_total_circle_length(int theta_count)
{
	if (theta_count <= 2) {
		return 0;
	} else {
		return TEXTURE_MAX_ANGLE *
			sin(theta_count * PI / (2 * (theta_count - 1))) /
			sin(PI / (2 * (theta_count - 1)));
	}
}

inline double determine_circle_length(double theta)
{
	return TEXTURE_MAX_ANGLE * sin(theta * 2 * PI / TEXTURE_MAX_ANGLE);
}

// ordinary functions in alphabetical order

void block_mult_error(int N1, int N2, complex * x, double min_err,
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

void calculate_grid(int h_phi_count, int h_theta_count, int r_phi_count,
										int r_theta_count, double *h_phi, double *h_theta,
										double *r, int grid)
{
	int N1, N2;

	N1 = h_phi_count * h_theta_count;
	N2 = r_phi_count * r_theta_count;

	switch (grid) {
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
		case 2:
		{
			int theta_count = determine_theta_count(N1);
			int nodes_defined = 0;
			double remaining_circle_length =
				determine_total_circle_length(theta_count);
			int i;

			if (nodes_defined < N1) {
				h_phi[nodes_defined] = 0;
				h_theta[nodes_defined] = 0;
				nodes_defined++;
			}

			if (nodes_defined < N1) {
				h_phi[nodes_defined] = 0;
				h_theta[nodes_defined] = TEXTURE_MAX_ANGLE / 2;
				nodes_defined++;
			}

			for (i = 1; i < theta_count / 2; i++) {
				double theta_1 =
					equidist(0, TEXTURE_MAX_ANGLE / 2, i, theta_count, 1);
				double theta_2 = TEXTURE_MAX_ANGLE / 2 - theta_1;
				double circle_length = determine_circle_length(theta_1);
				double phi_count_double =
					(2 * circle_length / remaining_circle_length) * (N1 -
																													 nodes_defined);
				int phi_count =
					(phi_count_double - floor(phi_count_double) >=
					 0.5) ? ceil(phi_count_double) : floor(phi_count_double);
				int phi_count_1 = phi_count / 2;
				int phi_count_2 = phi_count - phi_count_1;
				int j;

				remaining_circle_length -= 2 * circle_length;

				for (j = 0; j < phi_count_1; j++) {
					h_phi[nodes_defined] =
						equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, j,
										 phi_count_1, 0);
					h_theta[nodes_defined] = theta_1;
					nodes_defined++;
				}

				for (j = 0; j < phi_count_2; j++) {
					h_phi[nodes_defined] =
						equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, j,
										 phi_count_2, 0);
					h_theta[nodes_defined] = theta_2;
					nodes_defined++;
				}
			}

			if (theta_count % 2 == 1 && theta_count >= 3) {
				double theta = TEXTURE_MAX_ANGLE / 4;
				int phi_count = N1 - nodes_defined;
				int j;

				for (j = 0; j < phi_count; j++) {
					h_phi[nodes_defined] =
						equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, j,
										 phi_count, 0);
					h_theta[nodes_defined] = theta;
					nodes_defined++;
				}
			}

			if (N1 > 0) {
				int theta_count = determine_theta_count(N2);
				int nodes_defined = 0;
				double remaining_circle_length =
					determine_total_circle_length(theta_count);
				int i;

				if (nodes_defined < N2) {
					r[2 * nodes_defined] = 0;
					r[2 * nodes_defined + 1] = 0;
					nodes_defined++;
				}

				if (nodes_defined < N2) {
					r[2 * nodes_defined] = 0;
					r[2 * nodes_defined + 1] = TEXTURE_MAX_ANGLE / 2;
					nodes_defined++;
				}

				for (i = 1; i < theta_count / 2; i++) {
					double theta_1 =
						equidist(0, TEXTURE_MAX_ANGLE / 2, i, theta_count, 1);
					double theta_2 = TEXTURE_MAX_ANGLE / 2 - theta_1;
					double circle_length = determine_circle_length(theta_1);
					double phi_count_double =
						(2 * circle_length / remaining_circle_length) * (N2 -
																														 nodes_defined);
					int phi_count =
						(phi_count_double - floor(phi_count_double) >=
						 0.5) ? ceil(phi_count_double) : floor(phi_count_double);
					int phi_count_1 = phi_count / 2;
					int phi_count_2 = phi_count - phi_count_1;
					int j;

					remaining_circle_length -= 2 * circle_length;

					for (j = 0; j < phi_count_1; j++) {
						r[2 * nodes_defined] =
							equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, j,
											 phi_count_1, 0);
						r[2 * nodes_defined + 1] = theta_1;
						nodes_defined++;
					}

					for (j = 0; j < phi_count_2; j++) {
						r[2 * nodes_defined] =
							equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, j,
											 phi_count_2, 0);
						r[2 * nodes_defined + 1] = theta_2;
						nodes_defined++;
					}
				}

				if (theta_count % 2 == 1 && theta_count >= 3) {
					double theta = TEXTURE_MAX_ANGLE / 4;
					int phi_count = N2 - nodes_defined;
					int j;

					for (j = 0; j < phi_count; j++) {
						r[2 * nodes_defined] =
							equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, j,
											 phi_count, 0);
						r[2 * nodes_defined + 1] = theta;
						nodes_defined++;
					}
				}

				{
					int i;
					for (i = 1; i < N1; i++) {
						memcpy(&(r[2 * i * N2]), r, 2 * N2 * sizeof(double));
					}
				}
			}
			break;
		}
		default:
		{
			fprintf(stderr, "Illegal grid type (%d)!", grid);
			fflush(0);
			exit(-1);
		}
	}
}

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

inline int equal(complex x, complex y, double delta)
{
	return cabs(x - y) <= delta;
}

inline double equidist(double start, double end, int i, int n, int incl)
{
	if (!incl) {
		n++;
	}
	return start + (end - start) * (double) i / (double) (n - 1);
}

inline void error(const char *msg)
{
	fprintf(stderr, "%s", msg);
	fflush(0);
	exit(-1);
}

inline complex expi(double phi)
{
	return cos(phi) + I * sin(phi);
}

void init_omega(complex * omega, int N, int omega_policy)
{
	int i;
	switch (omega_policy) {
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

double l_2_rel_dist(const complex * vec, const complex * ref,
										unsigned int length)
{
	double abs_dist = l_2_dist(vec, ref, length);
	double base = l_2_norm(ref, length);
	if (base != 0.0 || abs_dist != 0.0) {
		return abs_dist / base;
	} else {
		return 0.0;
	}
}

void mult_error(int N1, int N2, complex * x, double min_err, double max_err)
{
	int i, j;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			x[i * N2 + j] *= min_err + drand48() * (max_err - min_err);
		}
	}
}

void read_grid(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
							 double **h_theta_ptr, double **r_ptr, FILE * h_in, FILE * r_in,
							 FILE * out)
{
	read_h(N1_ptr, h_phi_ptr, h_theta_ptr, h_in, out);
	read_r(*N1_ptr, N2_ptr, r_ptr, r_in, out);
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

	*h_phi_ptr = (double *) smart_malloc(N1 * sizeof(double));
	*h_theta_ptr = (double *) smart_malloc(N1 * sizeof(double));
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

	*omega_ptr =
		(complex *) smart_malloc(texture_flat_length(N) * sizeof(complex));
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

	*r_ptr = (double *) smart_malloc(N1 * N2 * 2 * sizeof(double));
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

	*x_ptr = (complex *) smart_malloc(N1 * N2 * sizeof(complex));
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

#undef malloc
#undef calloc

inline void *smart_calloc(size_t nmemb, size_t size)
{
	void *buf = calloc(nmemb, size);
	if (!buf) {
		error("smart_calloc: out of memory!\n");
	}
	return buf;
}

inline void *smart_malloc(size_t size)
{
	void *buf = malloc(size);
	if (!buf) {
		error("smart_calloc: out of memory!\n");
	}
	return buf;
}

#define malloc(a) arglasdf(a)
#define calloc(a) arglasdf(a)

complex spherical_harmonic(int k, int n, double phi, double theta)
{
	static double p;
	static double p_old;
	static double old_theta;
	static int k0;
	static int n0 = -1;

	if (abs(n) != n0 || k < k0 || !equal(theta, old_theta, 1E-15)) {
		k0 = abs(n);
		old_theta = theta;
		p_old = 0;
		p = 1;

		for (n0 = 0; n0 < abs(n); n0++) {
			p *= sqrt((double) (2 * n0 + 1) / (double) (2 * n0 + 2))
				* sin(theta * 2 * PI / TEXTURE_MAX_ANGLE);
		}
	}

	for (; k0 < k; k0++) {
		double p_new =
			(double) (2 * k0 + 1) / sqrt((double) ((k0 - n0 + 1) * (k0 + n0 + 1)))
			* cos(theta * 2 * PI / TEXTURE_MAX_ANGLE) * p
			- sqrt((double) ((k0 - n0) * (k0 + n0)) /
						 (double) ((k0 - n0 + 1) * (k0 + n0 + 1)))
			* p_old;
		p_old = p;
		p = p_new;
	}

	return p * expi(n * phi * 2 * PI / TEXTURE_MAX_ANGLE);
}

inline void split(int n, int *t1, int *t2)
{
	if (n == 0) {
		*t1 = 0;
		*t2 = 0;
	} else if (n > 0) {
		*t1 = ceil(sqrt(n));
		while (n % *t1 != 0) {
			(*t1)--;
		}
		*t2 = n / *t1;
	} else {
		error("Try to split a negative number!\n");
	}
}

void set_weights(itexture_plan * iplan, int weight_policy)
{
	int i;
	switch (weight_policy) {
		case 0:
			break;
		case 1:
		{
			for (i = 1; i <= iplan->mv->N_total; i++) {
				iplan->w_hat[i - 1] = (double) (1) / (double) (i);
			}
			break;
		}
		case 2:
		{
			for (i = 1; i <= iplan->mv->N_total; i++) {
				iplan->w_hat[i - 1] = (double) (1) / ((double) (i) * (double) (i));
			}
			break;
		}
		default:
			error("Illegal weight policy!\n");
	}
}

unsigned int solver_flags(int solver_algo, int weight_policy)
{
	if (solver_algo == 0) {
		return CGNR;
	} else {
		if (weight_policy == 0) {
			return CGNE;
		} else {
			return CGNE | PRECOMPUTE_DAMP;
		}
	}
}

inline void warning(const char *message)
{
	printf("Warning: %s\n", message);
	fflush(0);
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

void write_r(int N2, double *r, FILE * out)
{
	int j;

	fprintf(out, "\n");
	fprintf(out, "%d\n", N2);

	for (j = 0; j < N2; j++) {
		fprintf(out, "%lg %lg\n", r[2 * j], r[2 * j + 1]);
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
