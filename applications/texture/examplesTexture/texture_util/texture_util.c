#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<math.h>

#include"texture_util.h"
#include<util.h>
// macros

// descriptors
const char *grid_descr[] =
	{ "equidistant angles", "file", "uniformly distributed" };

const char *omega_policy_descr[] = { "flat", "1/n", "one", "1/n^2" };

const char *solver_algo_descr[] = { "CGNR", "CGNE" };

const char *weight_policy_descr[] = { "flat", "1/n", "1/n^2", "1/n^3" };

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

void calculate_grid(grid_dim dims, double *h_phi, double *h_theta,
										double *r, int grid)
{
	switch (grid) {
		case 0:
		{
			int i, j;
			int s, t;
			int N1, N2;
			int h_phi_count = dims.angles.h_phi_count;
			int h_theta_count = dims.angles.h_theta_count;
			int r_phi_count = dims.angles.r_phi_count;
			int r_theta_count = dims.angles.r_theta_count;

			N1 = h_phi_count * (h_theta_count - 2) + 2;
			N2 = r_phi_count * (r_theta_count - 2) + 2;

			for (s = 0; s < h_theta_count - 2; s++) {
				double theta =
					equidist(0, TEXTURE_MAX_ANGLE / 2, s + 1, h_theta_count, 1);
				for (t = 0; t < h_phi_count; t++) {
					double phi =
						equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, t,
										 h_phi_count, 0);
					h_phi[s * h_phi_count + t] = phi;
					h_theta[s * h_phi_count + t] = theta;
				}
			}

			h_phi[N1 - 2] = 0;
			h_theta[N1 - 2] = 0;
			h_phi[N1 - 1] = 0;
			h_theta[N1 - 1] = TEXTURE_MAX_ANGLE / 2;

			j = 0;
			for (s = 0; s < r_theta_count - 2; s++) {
				double theta =
					equidist(0, TEXTURE_MAX_ANGLE / 2, s + 1, r_theta_count, 1);
				for (t = 0; t < r_phi_count; t++) {
					double phi =
						equidist(-TEXTURE_MAX_ANGLE / 2, TEXTURE_MAX_ANGLE / 2, t,
										 r_phi_count, 0);
					r[j++] = phi;
					r[j++] = theta;
				}
			}

			r[j++] = 0;
			r[j++] = 0;
			r[j++] = 0;
			r[j++] = TEXTURE_MAX_ANGLE / 2;

			for (i = 1; i < N1; i++) {
				memcpy(&(r[2 * i * N2]), r, 2 * N2 * sizeof(double));
			}
			break;
		}
		case 1:
			error("Cannot calculate a grid from type file!\n");
		case 2:
		{
			int N1 = dims.samples.N1;
			int N2 = dims.samples.N2;
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

void destroy_itexture_params(itexture_params * pars)
{
	free(pars->omega_min_res);
	pars->status = "destroyed";
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
				omega[i] =
					(drand48() - 0.5) * rand() + I * (drand48() - 0.5) * rand();
			}
			break;
		}
		case 1:
		{
			for (i = 0; i < texture_flat_length(N); i++) {
				omega[i] =
					((drand48() - 0.5) * rand() + I * (drand48() - 0.5) * rand()) / (i +
																																					 1);
			}
			break;
		}
		case 2:
		{
			for (i = 0; i < texture_flat_length(N); i++) {
				omega[i] = 1;
			}
			break;
		}
		case 3:
		{
			for (i = 0; i < texture_flat_length(N); i++) {
				omega[i] =
					((drand48() - 0.5) * rand() + I * (drand48() - 0.5) * rand());
				omega[i] /= (i + 1);
				omega[i] /= (i + 1);
			}
			break;
		}
		default:
			error("illegal omega policy\n");
	}
}

void initialize_itexture_params(itexture_params * pars, int N)
{
	pars->omega_ref = 0;

	pars->max_epochs = 10000;
	pars->stop_file_name = "stop_solver";
	pars->residuum_goal = 1e-16;
	pars->updated_residuum_limit = 1e-16;
	pars->min_improve = 0.01;
	pars->max_epochs_without_improve = 10;
	pars->max_fail = 0;
	pars->steps_per_epoch = 1;
	pars->use_updated_residuum = 0;
	pars->monitor_error = 0;

	pars->messages_on = 1;
	pars->message_interval = 1;

	pars->omega_min_res =
		smart_malloc(texture_flat_length(N) * sizeof(double complex));
	assert(pars->omega_min_res != 0);
	pars->epochs_until_min_res = -1;

	pars->omega_min_err = 0;
	pars->epochs_until_min_err = -1;

	pars->status = "initialized";
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

inline double l_2_rel_norm(const complex * vec, const complex * ref,
													 unsigned int length)
{
	double dist = l_2_norm(vec, length);
	double base = l_2_norm(ref, length);
	if (base != 0.0 || dist != 0.0) {
		return dist / base;
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
		case 3:
		{
			for (i = 1; i <= iplan->mv->N_total; i++) {
				iplan->w_hat[i - 1] = pow(i, -3);
			}
			break;
		}
		default:
			error("Illegal weight policy!\n");
	}
}

unsigned int solver_flags(int solver_algo, int weight_policy)
{
	unsigned int flags;
	if (solver_algo == 0) {
		flags = CGNR;
	} else {
		flags = CGNE;
	}
	if (weight_policy != 0) {
		flags |= PRECOMPUTE_DAMP;
	}
	return flags;
}

inline void split(int n, int *t1, int *t2)
{
	n -= 2;
	if (n == 0) {
		*t1 = 0;
		*t2 = 2;
	} else if (n > 0) {
		*t1 = ceil(sqrt(n));
		while (n % *t1 != 0) {
			(*t1)++;
		}
		*t2 = n / *t1;

		(*t2) += 2;
	} else {
		error("Try to split a number less than 2!\n");
	}
}

void texture_itrafo(itexture_plan * iplan, itexture_params * pars)
{
	int stop = 0;
	int epoch = 0, epochs_without_improve = 0, epochs_with_fail = 0;
	texture_plan *test_plan = iplan->mv;
	double y_norm = l_2_norm(iplan->y, texture_get_x_length(test_plan));

	memset(iplan->f_hat_iter, 0,
				 texture_get_omega_length(test_plan) * sizeof(double complex));

	itexture_before_loop(iplan);

	memcpy(pars->omega_min_res, iplan->f_hat_iter,
				 texture_get_omega_length(test_plan) * sizeof(double complex));
	pars->epochs_until_min_res = 0;
	pars->min_residuum = 1;

	if (pars->monitor_error) {
		pars->error_during_min_residuum = 1;
		memcpy(pars->omega_min_err, iplan->f_hat_iter,
					 texture_get_omega_length(test_plan) * sizeof(double complex));
		pars->min_error = 1;
		pars->epochs_until_min_err = 0;
	}

	if (pars->messages_on) {
		fprintf(stderr, "Create \"%s\" to exit.\n", pars->stop_file_name);
	}
	// Loop performing iterations of the solver.
	do {
		int count;
		double err = 10, res_upd;
		static double res = 1, res_old;
		FILE *stop_file;
		int better_err = 0, better_res = 0;
		const char *marker[] = { "-", "*", ":" };

		for (count = 0; count < pars->steps_per_epoch; count++) {
			itexture_loop_one_step(iplan);
		}
		epoch++;

		// Calculate residuum and error.

		res_old = res;
		res_upd = sqrt(iplan->dot_r_iter) / y_norm;

		if (pars->use_updated_residuum) {
			res = res_upd;
		} else {
			double complex * tmp = test_plan->f_hat;
			
			texture_set_omega(test_plan, iplan->f_hat_iter);
			texture_trafo(test_plan);
			texture_set_omega(test_plan, tmp);
			
			res =
				l_2_rel_dist(texture_get_x(test_plan), iplan->y,
										 texture_get_x_length(test_plan));
		}

		if ((pars->min_residuum - res) / pars->min_residuum < pars->min_improve) {
			epochs_without_improve++;
		} else {
			better_res = 2;
			epochs_without_improve = 0;
		}

		if (res_old < res) {
			epochs_with_fail++;
		} else {
			epochs_with_fail = 0;
		}

		if (pars->monitor_error) {
			err =
				l_2_rel_dist(iplan->f_hat_iter, pars->omega_ref,
										 texture_get_omega_length(test_plan));
		}

		if (res < pars->min_residuum) {
			better_res = MAX(1, better_res);
			memcpy(pars->omega_min_res, iplan->f_hat_iter,
						 texture_get_omega_length(test_plan) * sizeof(double complex));
			pars->epochs_until_min_res = epoch;
			pars->min_residuum = res;
			if (pars->monitor_error) {
				pars->error_during_min_residuum = err;
			}
		}

		if (pars->monitor_error && err < pars->min_error) {
			better_err = 2;
			memcpy(pars->omega_min_err, iplan->f_hat_iter,
						 texture_get_omega_length(test_plan) * sizeof(double complex));
			pars->min_error = err;
			pars->epochs_until_min_err = epoch;
		}
		// Messages

		if (pars->messages_on && epoch % pars->message_interval == 0) {
			fprintf(stderr, "step: %4d, ", pars->steps_per_epoch * epoch);
			if (!pars->use_updated_residuum) {
				fprintf(stderr, "res%s %.2e, ", marker[better_res], res);
			}
			fprintf(stderr, "upd res%s %.2e", marker[better_res], res_upd);
			if (pars->monitor_error) {
				fprintf(stderr, ", err%s %.2e", marker[better_err], err);
			}
			fprintf(stderr, "\n");
		}
		// Check abort conditions.

		stop_file = fopen(pars->stop_file_name, "r");
		if (stop_file) {
			fclose(stop_file);
			stop = 1;
			pars->status = "stopped";
		}

		if (epoch >= pars->max_epochs) {
			pars->status = "max epochs reached";
			stop = 1;
		}

		if (res_upd <= pars->updated_residuum_limit) {
			pars->status = "stopped to prevent an underflow";
			stop = 1;
		}

		if (epochs_without_improve > pars->max_epochs_without_improve) {
			pars->status = "stagnation";
			stop = 1;
		}

		if (epochs_with_fail > pars->max_fail) {
			pars->status = "degradation";
			stop = 1;
		}

		if (res <= pars->residuum_goal) {
			pars->status = "residuum goal reached";
			stop = 1;
		}

	} while (!stop);

	if (pars->messages_on) {
		fprintf(stderr, "Abort: %s\n", pars->status);
		if (pars->monitor_error) {
			fprintf(stderr, "minimal residuum: %.2e, error: %.2e (steps: %d)\n",
							pars->min_residuum, pars->error_during_min_residuum,
							pars->epochs_until_min_res * pars->steps_per_epoch);
			fprintf(stderr, "minimal error: %.2e, (steps: %d)\n", pars->min_error,
							pars->epochs_until_min_err * pars->steps_per_epoch);
		} else {
			fprintf(stderr, "minimal residuum: %.2e, (steps: %d)\n",
							pars->min_residuum,
							pars->epochs_until_min_res * pars->steps_per_epoch);
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
		fprintf(out, "%24.16e %24.16e\n", h_phi[i], h_theta[i]);
	}
}

void write_omega(int N, complex * omega, FILE * out)
{
	int i;

	fprintf(out, "\n");
	fprintf(out, "%d\n", N);
	for (i = 0; i < texture_flat_length(N); i++) {
		fprintf(out, "%24.16e + %24.16ei\n", creal(omega[i]), cimag(omega[i]));
	}
}

void write_r(int N2, double *r, FILE * out)
{
	int j;

	fprintf(out, "\n");
	fprintf(out, "%d\n", N2);

	for (j = 0; j < N2; j++) {
		fprintf(out, "%24.16e %24.16e\n", r[2 * j], r[2 * j + 1]);
	}
}

void write_x(int N1, int N2, complex * x, FILE * out)
{
	int i, j;

	fprintf(out, "\n");
	fprintf(out, "%d %d\n", N1, N2);
	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			fprintf(out, "%24.16e + %24.16ei\n", creal(x[i * N2 + j]),
							cimag(x[i * N2 + j]));
		}
	}
}
