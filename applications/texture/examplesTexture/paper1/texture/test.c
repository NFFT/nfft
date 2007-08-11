#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<math.h>

#include<texture_util.h>

/** 
 * @defgroup texture_paper1_texture Texture
 * 
 * The program reads a grid h, r and frequencies omega_ref.
 * Parameters are requested interactively.
 * First, the program applies ::texture_trafo to calculate the samples x_ref at
 * the points of h, r.
 * Afterwards, it runs the solver on x_ref, using ::texture_itrafo of 
 * @ref texture_utility.
 * We call the resulting frequencies omega.
 * For each l, the program outputs 
 * @f$\|omega\_ref(l,\cdot,\cdot) - omega(l,\cdot,\cdot)\|_2 / \|omega\_ref(l,\cdot,\cdot)\|_2@f$
 * to some output file.
 * The output file has the form 
 * \<prefix>\<identifier>.N_<N>.N1_<N1>.N2_<N2>.newN_<newN>.\<weight_policy>.
 *
 * @author Matthias Schmalz
 * @ingroup texture_paper1
 */

// fixed parameters

// 4: (2l+1)^2
int weight_policy = 4;
int enforce_even_interpolation = 1;
const char *output_prefix = "output/";
const char *default_name = "test";
// CGNR
int solver_algo = 0;
const char *def_h_file = "h28";
const char *def_r_file = "r9791";
const char *def_omega_file = "omega";

// input parameters

char h_file[100], r_file[100], omega_file[100];
int new_N;
char name[1000];

// program variables 

double _Complex *omega_ref, *omega_ref_even, *omega, *x_ref, *x;
double *h_phi, *h_theta, *r;
int N, N1, N2;

texture_plan plan, plan2;
itexture_plan iplan;
itexture_params pars;

FILE *output;

double *error_harm;

void read_data()
{
	FILE *f1, *f2, *f3;
	int l, m, n;
	char next;

	fprintf(stderr, "1. grid_h_file (%s): ", def_h_file);
	next = getchar();

	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%s%*[ ]", h_file);
		getchar();
	} else {
		sprintf(h_file, "%s", def_h_file);
	}

	fprintf(stderr, "2. grid_r_file (%s): ", def_r_file);
	next = getchar();

	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%s%*[ ]", r_file);
		getchar();
	} else {
		sprintf(r_file, "%s", def_r_file);
	}

	fprintf(stderr, "3. omega_file (%s): ", def_omega_file);
	next = getchar();

	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%s%*[ ]", omega_file);
		getchar();
	} else {
		sprintf(omega_file, "%s", def_omega_file);
	}

	f1 = fopen(h_file, "r");
	f2 = fopen(r_file, "r");
	f3 = fopen(omega_file, "r");

	read_grid(&N1, &N2, &h_phi, &h_theta, &r, f1, f2, stderr);
	read_omega(&N, &omega_ref, f3, stderr);

	omega_ref_even =
		smart_malloc(texture_flat_length(N) * sizeof(double _Complex));
	memcpy(omega_ref_even, omega_ref,
				 texture_flat_length(N) * sizeof(double _Complex));
	for (l = 1; l <= N; l += 2) {
		for (m = -l; m <= l; m++) {
			for (n = -l; n <= l; n++) {
				omega_ref_even[texture_flat_index(l, m, n)] = 0;
			}
		}
	}

	fclose(f1);
	fclose(f2);
	fclose(f3);
}

void read_params()
{
	char next;

	do {
		fprintf(stderr, "Bandwidth for interpolation new_N = ");
		scanf("%d%*[ ]", &new_N);
		getchar();
	} while (new_N > N);

	initialize_itexture_params(&pars, new_N);
	read_itexture_params(&pars);

	pars.monitor_error = 1;
	pars.omega_ref = omega_ref_even;
	pars.omega_min_err =
		smart_malloc(texture_flat_length(new_N) * sizeof(double _Complex));

	fprintf(stderr, "Input a name for the data set (default=%s): ",
					default_name);
	next = getchar();

	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%s%*[ ]", name);
		getchar();
	} else {
		strcpy(name, default_name);
	}
}

void output_params()
{
	char output_file[1000];

	sprintf(output_file, "%s%s.N_%d.N1_%d.N2_%d.newN_%d.%s", output_prefix,
					name, N, N1, N2, new_N, weight_policy_descr[weight_policy]);
	output = fopen(output_file, "w");

	fprintf(output, "# Name: %s\n", name);
	fprintf(output, "# N: %d, dim(J_N): %d\n", N, texture_flat_length(N));
	fprintf(output, "# N1: %d, N2: %d, N1*N2: %d\n", N1, N2, N1 * N2);
	fprintf(output, "# new_N: %d, dim(J_new_N): %d\n", new_N,
					texture_flat_length(new_N));
	fprintf(output, "# Weight policy: %d (%s)\n", weight_policy,
					weight_policy_descr[weight_policy]);
	fprintf(output, "# enforce even: %d\n", enforce_even_interpolation);
	fprintf(output, "#\n");

	fprintf(output, "# Solver parameters:\n");
	fprintf(output, "# max_epochs: %d\n", pars.max_epochs);
	fprintf(output, "# residuum_goal: %.2e\n", pars.residuum_goal);
	fprintf(output,
					"# min_improve: %lg, max_epochs_without_improve: %d, max_fail: %d\n",
					pars.min_improve, pars.max_epochs_without_improve, pars.max_fail);
	fprintf(output, "# steps_per_epoch: %d\n", pars.steps_per_epoch);
	fprintf(output, "# use_updated_residuum: %d\n", pars.use_updated_residuum);
	fflush(output);
}

void calculate_x()
{
	x_ref = smart_malloc(N1 * N2 * sizeof(double _Complex));
	texture_precompute(N);
	texture_init(&plan, N, N1, N2, omega_ref_even, x_ref, h_phi, h_theta, r);
	texture_trafo(&plan);
}

void calculate_omega()
{
	int m, n, l;

	omega = smart_malloc(texture_flat_length(new_N) * sizeof(double _Complex));
	x = smart_malloc(N1 * N2 * sizeof(double _Complex));
	texture_init(&plan2, new_N, N1, N2, omega, x, h_phi, h_theta, r);
	itexture_init_advanced(&iplan, &plan2,
												 solver_flags(solver_algo, weight_policy));

	set_weights(&iplan, weight_policy);
	if (enforce_even_interpolation) {
		for (l = 1; l <= new_N; l += 2) {
			for (n = -l; n <= l; n++) {
				for (m = -l; m <= l; m++) {
					iplan.w_hat[texture_flat_index(l, m, n)] = 0;
				}
			}
		}
	}

	memcpy(iplan.y, x_ref,
				 texture_get_x_length(&plan2) * sizeof(double _Complex));

	texture_itrafo(&iplan, &pars);
}

void calculate_results()
{
	int l;
	// double sum = 0;

	error_harm = smart_malloc((N + 1) * sizeof(double));
	for (l = 0; l <= new_N; l++) {
		int start = texture_flat_index(l, -l, -l);
		int end = texture_flat_index(l, l, l);
		int len = end - start + 1;
		error_harm[l] =
			l_2_rel_dist(&(pars.omega_min_res[start]), &(omega_ref[start]), len);
	}
	for (l = new_N + 1; l <= N; l++) {
		error_harm[l] = 1;
	}

	/* for (l = 0; l <= new_N; l++) { sum += error_harm[l] * error_harm[l]; }
	   sum = sqrt(sum); sum *= ref; sum /= l_2_norm(omega_ref,
	   texture_flat_length(new_N));

	   fprintf(stderr, "Error sum: %.4e, Minimal error: %.4e, Difference:
	   %.4e\n", sum, pars.min_error, sum - pars.min_error); fprintf(stderr,
	   "new error sum: %.4e\n", l_2_rel_dist(pars.omega_min_err, omega_ref,
	   texture_flat_length(new_N))); */
}

void output_results()
{
	int width1 = 8, width2 = 9;
	int l;

	fprintf(output, "#\n");
	fprintf(output, "# Reason for abort: %s\n", pars.status);
	fprintf(output, "# After %d iterations:\n",
					pars.epochs_until_min_res * pars.steps_per_epoch);
	fprintf(output, "# Best residuum: %.2e, Error at the same time: %.2e.\n",
					pars.min_residuum, pars.error_during_min_residuum);
	fprintf(output, "# After %d iterations:\n",
					pars.epochs_until_min_err * pars.steps_per_epoch);
	fprintf(output, "# Best error: %.2e\n", pars.min_error);

	fprintf(output, "#\n");

	fprintf(output,
					"# best frequency vector: omega_min_res (not omega_min_err)\n");
	fprintf(output,
					"# For each degree l we give the relative l_2 error between the best frequency\n");
	fprintf(output,
					"# vector omega_min and the reference vector omega_ref projected to Harm_l, i.e.\n");
	fprintf(output,
					"# sqrt(sum_{m,n} |omega_min[l,m,n]-omega_ref[l,m,n]|^2) /\n");
	fprintf(output, "# sqrt(sum_{m,n} |omega_ref[l,m,n]|^2).\n");

	fprintf(output, "# %*s %*s\n", width1 - 2, "degree", width2, "error");
	for (l = 0; l <= N; l++) {
		fprintf(output, "%*d %*.2e\n", width1, l, width2, error_harm[l]);
	}
}

void cleanup()
{
	// calculate_results
	free(error_harm);

	// calculate_omega
	itexture_finalize(&iplan);
	texture_finalize(&plan2);
	free(x);
	free(omega);

	// calculate_x
	texture_finalize(&plan);
	texture_forget();
	free(x_ref);

	// output_params
	fclose(output);

	// read_params
	free(pars.omega_min_err);
	destroy_itexture_params(&pars);

	// read_data
	free(h_phi);
	free(h_theta);
	free(r);
	free(omega_ref);
	free(omega_ref_even);
}

int main()
{
	read_data();

	read_params();

	output_params();

	calculate_x();

	fprintf(stderr, "Calculating omega ...\n");
	fflush(stderr);

	calculate_omega();

	calculate_results();

	output_results();

	cleanup();

	return 0;
}
