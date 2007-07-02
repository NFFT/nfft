#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<ctype.h>
#include<assert.h>
#include<string.h>
#include<math.h>

#include<texture_util.h>

// input parameters
int testcases = 100, firstcase = 0;
int N;
char h_file_name[100], r_file_name[100];
char omega_file_name[100];
const char *def_omega_file_name = "omega";
char output_file_name[100];
const char *def_output_file_name = "test";
int solver_algo;
itexture_params pars;

// fixed parameters
const char *grid_prefix = "grid/";
const char *omega_prefix = "omega/";
const char *output_prefix = "output/";
// 4: (2l+1)^2
int weight_policy = 4;

// program data
char cur_omega_file_path[200];
char h_file_path[200];
char r_file_path[200];
char output_file_path[200];

FILE *output_file;

int N1, N2;
double *h_phi, *h_theta, *r;
double complex *omega;

texture_plan plan;
itexture_plan iplan;

double *rrn, *ren, *rwen;
double rrn_min, rrn_max, rrn_mean;
double ren_min, ren_max, ren_mean;
double rwen_min, rwen_max, rwen_mean;

void read_data()
{
	char next;
	FILE *data_file;
	int cur_N, i;

	// testcases
	fprintf(stderr, "testcases (%d): ", testcases);
	next = getchar();
	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%d", &testcases);
		getchar();
	}
	// firstcase
	fprintf(stderr, "firstcase (%d): ", firstcase);
	next = getchar();
	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%d", &firstcase);
		getchar();
	}
	// N
	fprintf(stderr, "N: ");
	scanf("%d", &N);
	getchar();

	// omega_name
	fprintf(stderr, "omega_file_name (%s): ", def_omega_file_name);
	next = getchar();
	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%s%*[ ]", omega_file_name);
		getchar();
	} else {
		sprintf(omega_file_name, "%s", def_omega_file_name);
	}

	for (i = firstcase; i < firstcase + testcases; i++) {
		sprintf(cur_omega_file_path, "%s%s%03d", omega_prefix, omega_file_name,
						i);
		data_file = fopen(cur_omega_file_path, "r");

		read_N(&cur_N, data_file);
		assert(cur_N >= N);

		fclose(data_file);
	}

	// h_file_name
	fprintf(stderr, "h_file_name: ");
	scanf("%s%*[ ]", h_file_name);
	getchar();

	sprintf(h_file_path, "%s%s", grid_prefix, h_file_name);
	data_file = fopen(h_file_path, "r");
	read_h(&N1, &h_phi, &h_theta, data_file, stderr);
	fclose(data_file);

	// r_file_name
	fprintf(stderr, "r_file_name: ");
	scanf("%s%*[ ]", r_file_name);
	getchar();

	sprintf(r_file_path, "%s%s", grid_prefix, r_file_name);
	data_file = fopen(r_file_path, "r");
	read_r(N1, &N2, &r, data_file, stderr);
	fclose(data_file);
}

void read_params()
{
	char next;

	// output_file
	fprintf(stderr, "output_file_name (%s): ", def_output_file_name);
	next = getchar();
	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%s%*[ ]", output_file_name);
		getchar();
	} else {
		sprintf(output_file_name, "%s", def_output_file_name);
	}

	// solver_algo
	do {
		fprintf(stderr, "solver_algo (0 = %s, 1 = %s): ", solver_algo_descr[0],
						solver_algo_descr[1]);
		scanf("%d", &solver_algo);
		getchar();
	} while (solver_algo != 0 && solver_algo != 1);

	// solver_params
	initialize_itexture_params(&pars, N);
	read_itexture_params(&pars);

	pars.monitor_error = 1;
	pars.omega_min_err =
		smart_malloc(texture_flat_length(N) * sizeof(double complex));
}

void output_params()
{
	sprintf(output_file_path, "%s%s.N_%d.N1_%d.N2_%d.%s.%s.%02d-%02d",
					output_prefix, output_file_name, N, N1, N2,
					solver_algo_descr[solver_algo], weight_policy_descr[weight_policy],
					firstcase, firstcase + testcases - 1);
	output_file = fopen(output_file_path, "w");

	fprintf(output_file, "# Name: %s\n", output_file_name);
	fprintf(output_file, "# N: %d, dim(J_N): %d\n", N, texture_flat_length(N));
	fprintf(output_file, "# omega_file_name: %s\n", omega_file_name);
	fprintf(output_file, "# N1: %d, N2: %d, N1*N2: %d\n", N1, N2, N1 * N2);
	fprintf(output_file, "# h_file_name: %s, r_file_name: %s\n", h_file_name,
					r_file_name);
	fprintf(output_file, "# solver_algo: %s\n", solver_algo_descr[solver_algo]);
	fprintf(output_file, "# testcases: %d\n", testcases);
	fprintf(output_file, "# firstcase: %d\n", firstcase);
	fprintf(output_file, "# weight_policy: %s\n",
					weight_policy_descr[weight_policy]);
	fprintf(output_file, "#\n");

	fprintf(output_file, "# Solver parameters:\n");
	fprintf(output_file, "# max_epochs: %d\n", pars.max_epochs);
	fprintf(output_file,
					"# residuum_goal: %.2e\n",
					pars.residuum_goal);
	fprintf(output_file,
					"# min_improve: %lg, max_epochs_without_improve: %d, max_fail: %d\n",
					pars.min_improve, pars.max_epochs_without_improve, pars.max_fail);
	fprintf(output_file, "# steps_per_epoch: %d\n", pars.steps_per_epoch);
	fprintf(output_file, "# use_updated_residuum: %d\n",
					pars.use_updated_residuum);
	fprintf(output_file, "#\n");

	fflush(output_file);
}

void init()
{
	texture_init(&plan, N, N1, N2,
							 smart_malloc(texture_flat_length(N) * sizeof(double complex)),
							 smart_malloc(N1 * N2 * sizeof(double complex)), h_phi, h_theta,
							 r);
	texture_precompute(N);
	rrn = smart_malloc(testcases * sizeof(double));
	ren = smart_malloc(testcases * sizeof(double));
	rwen = smart_malloc(testcases * sizeof(double));
}

void processing()
{
	int count;

	for (count = 0; count < testcases; count++) {
		FILE *omega_file;
		int cur_N;
		double complex *tmp;

		// load omeag_ref
		sprintf(cur_omega_file_path, "%s%s%03d", omega_prefix, omega_file_name,
						firstcase + count);
		omega_file = fopen(cur_omega_file_path, "r");
		read_omega(&cur_N, &omega, omega_file, stderr);
		pars.omega_ref = omega;

		// compute right hand side
		tmp = texture_get_omega(&plan);
		texture_set_omega(&plan, omega);
		texture_trafo(&plan);
		texture_set_omega(&plan, tmp);

		// initialize solver
		itexture_init_advanced(&iplan, &plan,
													 solver_flags(solver_algo, weight_policy));
		set_weights(&iplan, weight_policy);
		memcpy(iplan.y, texture_get_x(&plan), N1 * N2 * sizeof(double complex));

		// run solver
		texture_itrafo(&iplan, &pars);

		// store information
		rrn[count] = pars.min_residuum;
		ren[count] = compute_ren(pars.omega_min_res, omega, N);
		rwen[count] = compute_rwen(pars.omega_min_res, omega, N);

		// cleanup
		itexture_finalize(&iplan);

		free(omega);
		fclose(omega_file);
	}
}

void calculate_output_helper(double *data, double *min_val, double *max_val,
														 double *mean_val)
{
	int i;

	if (testcases == 0) {
		return;
	} else {
		*min_val = data[0];
		*max_val = data[0];
		*mean_val = data[0];
		for (i = 1; i < testcases; i++) {
			if (data[i] < *min_val) {
				*min_val = data[i];
			} else if (data[i] > *max_val) {
				*max_val = data[i];
			}
			*mean_val *= data[i];
		}
	}
	*mean_val = pow(*mean_val, 1.0/(double) testcases);
}

void calculate_output()
{
	calculate_output_helper(rrn, &rrn_min, &rrn_max, &rrn_mean);
	calculate_output_helper(ren, &ren_min, &ren_max, &ren_mean);
	calculate_output_helper(rwen, &rwen_min, &rwen_max, &rwen_mean);
}

void output()
{
	int i;

	fprintf(output_file, "# Relative residuum norm:\n");
	fprintf(output_file, "#%10s%10s%10s\n", "Minimum", "Geomean", "Maximum");
	fprintf(output_file, " %10.2e%10.2e%10.2e\n", rrn_min, rrn_mean, rrn_max);
	fprintf(output_file, "# Relative L^2(SO(3)) error norm:\n");
	fprintf(output_file, "#%10s%10s%10s\n", "Minimum", "Geomean", "Maximum");
	fprintf(output_file, " %10.2e%10.2e%10.2e\n", ren_min, ren_mean, ren_max);
	fprintf(output_file, "# Relative (un)weighted error norm:\n");
	fprintf(output_file, "#%10s%10s%10s\n", "Minimum", "Geomean", "Maximum");
	fprintf(output_file, " %10.2e%10.2e%10.2e\n", rwen_min, rwen_mean,
					rwen_max);
	fprintf(output_file, "#\n");

	fprintf(output_file, "# Rawdata:\n");
	fprintf(output_file, "#%10s%10s%10s\n", "rrn", "ren", "rwen");
	for (i = 0; i < testcases; i++) {
		fprintf(output_file, " %10.2e%10.2e%10.2e\n", rrn[i], ren[i], rwen[i]);
	}
}

void clean_up()
{
	free(rwen);
	free(ren);
	free(rrn);
	texture_forget();
	free(texture_get_omega(&plan));
	free(texture_get_x(&plan));
	texture_finalize(&plan);

	fclose(output_file);

	free(pars.omega_min_err);
	destroy_itexture_params(&pars);

	free(h_phi);
	free(h_theta);
	free(r);
}

int main()
{
	read_data();

	read_params();

	output_params();

	init();

	processing();

	calculate_output();

	output();

	clean_up();

	return 0;
}
