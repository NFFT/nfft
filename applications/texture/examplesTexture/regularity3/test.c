#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<math.h>

#include<texture_util.h>

#define NFFT_MAX(a, b) (((a) > (b)) ? (a) : (b))

/** @defgroup texture_regularity3 Texture: Regularity 3
 * TODO
 * @author Matthias Schmalz
 * @ingroup texture_examples
 * @{
 */

double complex *omega_ref, *omega, *x_ref, *x;
double *h_phi, *h_theta, *r;
int N, N1, N2, new_N;
const char *h_file = "grid_h.in", *r_file = "grid_r.in", *omega_file =
	"omega.in";

texture_plan plan, plan2;
itexture_plan iplan;
itexture_params pars;

int weight_policy;

const char *output_prefix = "output/";
char name[1000];
const char *default_name = "testset";
FILE *output;

double *error_harm;

void read_data()
{
	FILE *f1, *f2, *f3;
	char next;
	do {
		int choice;
		char file_name[100];

		fprintf(stderr, "1. grid_h_file = %s\n", h_file);
		fprintf(stderr, "2. grid_r_file = %s\n", r_file);
		fprintf(stderr, "3. omega_file = %s\n", omega_file);
		next = getchar();

		if (!isspace(next)) {
			ungetc(next, stdin);
			scanf("%d %s%*[ ]", &choice, file_name);
			getchar();

			switch (choice) {
				case 1:
					h_file = file_name;
					break;
				case 2:
					r_file = file_name;
					break;
				case 3:
					omega_file = file_name;
					break;
			}
		}

	} while (!isspace(next));

	f1 = fopen(h_file, "r");
	f2 = fopen(r_file, "r");
	f3 = fopen(omega_file, "r");

	read_grid(&N1, &N2, &h_phi, &h_theta, &r, f1, f2, stderr);
	read_omega(&N, &omega_ref, f3, stderr);

	fclose(f1);
	fclose(f2);
	fclose(f3);
}

void read_params()
{
	char next;

	fprintf(stderr, "Bandwidth for interpolation new_N = ");
	scanf("%d%*[ ]", &new_N);
	getchar();
	initialize_itexture_params(&pars, new_N);
	do {
		int choice;

		fprintf(stderr, "1. max_epochs = %d\n", pars.max_epochs);
		fprintf(stderr, "2. stop_file_name = %s\n", pars.stop_file_name);
		fprintf(stderr, "3. residuum_goal = %.2e\n", pars.residuum_goal);
		fprintf(stderr, "4. updated_residuum_limit = %.2e\n",
						pars.updated_residuum_limit);
		fprintf(stderr, "5. min_improve = %lg\n", pars.min_improve);
		fprintf(stderr, "6. max_epochs_without_improve = %d\n",
						pars.max_epochs_without_improve);
		fprintf(stderr, "7. max_fail = %d\n", pars.max_fail);
		fprintf(stderr, "8. steps_per_epoch = %d\n", pars.steps_per_epoch);
		fprintf(stderr, "9. use_updated_residuum = %d\n",
						pars.use_updated_residuum);
		fprintf(stderr, "10. messages_on = %d\n", pars.messages_on);
		fprintf(stderr, "11. message_interval = %d\n", pars.message_interval);

		next = getchar();
		if (!isspace(next)) {
			int *fields[12] =
				{ 0, &pars.max_epochs, 0, 0, 0, 0, &pars.max_epochs_without_improve,
				&pars.max_fail, &pars.steps_per_epoch, &pars.use_updated_residuum,
				&pars.messages_on, &pars.message_interval
			};
			int value;

			ungetc(next, stdin);
			scanf("%d%*[ ]", &choice);
			switch (choice) {
				case 1:
				case 6:
				case 7:
				case 8:
				case 9:
				case 10:
				case 11:
					scanf("%d%*[ ]", &value);
					*fields[choice] = value;
					break;
				case 2:
					scanf("%s%*[ ]", pars.stop_file_name);
					break;
				case 3:
					scanf("%lg%*[ ]", &pars.residuum_goal);
					break;
				case 4:
					scanf("%lg%*[ ]", &pars.updated_residuum_limit);
					break;
				case 5:
					scanf("%lg%*[ ]", &pars.min_improve);
					break;
			}
			getchar();
		}
	} while (!isspace(next));

	pars.monitor_error = 1;
	pars.omega_ref = omega_ref;
	pars.omega_min_err =
		smart_malloc(texture_flat_length(new_N) * sizeof(double complex));

	fprintf(stderr,
					"Choose a weight policy (0 = flat, 1 = 1/n, 2 = 1/n^2 (default), 3 = 1/n^3): ");
	next = getchar();

	if (!isspace(next)) {
		ungetc(next, stdin);
		scanf("%d%*[ ]", &weight_policy);
		getchar();
	} else {
		weight_policy = 2;
	}

	fprintf(stderr, "Input a name for the data set (default=%s):",
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

	sprintf(output_file, "%s%s_N%d_N1%d_N2%d_newN%d_w%d", output_prefix, name,
					N, N1, N2, new_N, weight_policy);
	output = fopen(output_file, "w");

	fprintf(output, "# Name: %s\n", name);
	fprintf(output, "# N: %d, dim(J_N): %d\n", N, texture_flat_length(N));
	fprintf(output, "# N1: %d, N2: %d, N1*N2: %d\n", N1, N2, N1 * N2);
	fprintf(output, "# new_N: %d, dim(J_new_N): %d\n", new_N,
					texture_flat_length(new_N));
	fprintf(output, "# Weight policy: %d (%s)\n", weight_policy,
					weight_policy_descr[weight_policy]);

	fprintf(output, "#\n");

	fprintf(output, "# Solver parameters:\n");
	fprintf(output, "# max_epochs: %d\n", pars.max_epochs);
	fprintf(output, "# residuum_goal: %.2e, updated_residuum_limit: %.2e\n",
					pars.residuum_goal, pars.updated_residuum_limit);
	fprintf(output,
					"# min_improve: %lg, max_epochs_without_improve: %d, max_fail: %d\n",
					pars.min_improve, pars.max_epochs_without_improve, pars.max_fail);
	fprintf(output, "# steps_per_epoch: %d\n", pars.steps_per_epoch);
	fprintf(output, "# use_updated_residuum: %d\n", pars.use_updated_residuum);
}

void calculate_x()
{
	x_ref = smart_malloc(N1 * N2 * sizeof(double complex));
	texture_precompute(NFFT_MAX(N, new_N));
	texture_init(&plan, N, N1, N2, omega_ref, x_ref, h_phi, h_theta, r);
	texture_trafo(&plan);
}

void calculate_omega()
{
	omega = smart_malloc(texture_flat_length(new_N) * sizeof(double complex));
	x = smart_malloc(N1 * N2 * sizeof(double complex));
	texture_init(&plan2, new_N, N1, N2, omega, x, h_phi, h_theta, r);
	itexture_init_advanced(&iplan, &plan2, solver_flags(0, weight_policy));
	set_weights(&iplan, weight_policy);

	memcpy(iplan.y, x_ref,
				 texture_get_x_length(&plan2) * sizeof(double complex));

	texture_itrafo(&iplan, &pars);
}

void calculate_results()
{
	int l;
	double ref = l_2_norm(omega_ref, texture_flat_length(N));
	double sum = 0;

	error_harm = smart_malloc((N+1) * sizeof(double));
	for (l = 0; l <= new_N; l++) {
		int start = texture_flat_index(l, -l, -l);
		int end = texture_flat_index(l, l, l);
		int len = end - start + 1;
		error_harm[l] =
			l_2_dist(&(pars.omega_min_err[start]), &(omega_ref[start]), len);
	}
	for (l = new_N+1; l <= N; l++) {
		int start = texture_flat_index(l, -l, -l);
		int end = texture_flat_index(l, l, l);
		int len = end - start + 1;
		error_harm[l] =
			l_2_norm(&(omega_ref[start]), len);
	}
	for (l = 0; l <= N; l++) {
		error_harm[l] /= ref;
	}

	for (l = 0; l <= new_N; l++) {
		sum += error_harm[l] * error_harm[l];
	}
	sum = sqrt(sum);
	sum *= ref;
	sum /= l_2_norm(omega_ref, texture_flat_length(new_N));
	
	fprintf(stderr, "Error sum: %.4e, Minimal error: %.4e, Difference: %.4e\n",
					sum, pars.min_error, sum - pars.min_error);
	fprintf(stderr, "new error sum: %.4e\n", l_2_rel_dist(pars.omega_min_err, omega_ref, texture_flat_length(new_N)));
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
					"# For each degree l we give the relative l_2 error between the best frequency\n");
	fprintf(output,
					"# vector omega_min and the reference vector omega_ref projected to Harm_l, i.e.\n");
	fprintf(output,
					"# sqrt(sum_{m,n} (omega_min[l,m,n]-omega_ref[l,m,n])^2) / norm_2(omega_ref).\n");

	fprintf(output, "# %*s %*s\n", width1 - 2, "degree", width2, "error");
	for (l = 0; l <= N; l++) {
		fprintf(output, "%*d %*.2e\n", width1, l, width2, error_harm[l]);
	}
}

void cleanup()
{
	free(error_harm);

	fclose(output);

	texture_finalize(&plan2);
	free(omega);
	free(x);

	free(x_ref);
	texture_finalize(&plan);
	texture_forget();

	free(pars.omega_min_err);
	destroy_itexture_params(&pars);

	free(h_phi);
	free(h_theta);
	free(r);
	free(omega_ref);
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
/**
 * @}
 */
