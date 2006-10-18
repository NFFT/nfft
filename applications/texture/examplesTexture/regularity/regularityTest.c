#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include<nfft3_texture.h>
#include<util.h>
#include<texture_util.h>

const char *err_fun_id_descr[] = { "two norm relativ error" };
const char *res_fun_id_descr[] = { "two norm relativ error" };

typedef struct iteration_param_set_ {
	int err_fun_id;
	int res_fun_id;
	double res_delta;
	double err_delta;
	double min_improve;
	int iterations_without_check;
} iteration_param_set;

const char *starting_point_policy_descr[] = { "zero" };

typedef struct solver_param_set_ {
	int grid_id;
	char grid_file_h[100];
	char grid_file_r[100];
	int starting_point_policy;
	int solver_algorithm_id;
	int weight_policy;
	int omega_ref_policy;
} solver_param_set;

typedef struct test_case_param_set_0_ {
	int N1_start, N1_end, N1_incr;
	int N2_start, N2_end, N2_incr;
	int N_min_improve, max_iter_without_improve;
} test_case_param_set_0;

typedef struct test_case_param_set_1_ {
	int alpha;
	int N1_start, N1_end, N1_incr;
	int N_min_improve, max_iter_without_improve;
} test_case_param_set_1;

const char *method_descr[] = { "N1, N2 independent", "alpha * N1 = N2" };
const char *output_policy_descr[] = { "quiet", "human readable" };

typedef struct test_case_param_set_ {
	int method;
	test_case_param_set_0 c0;
	test_case_param_set_1 c1;
	int output_policy;
} test_case_param_set;

typedef struct param_set_ {
	iteration_param_set iteration;
	solver_param_set solver;
	test_case_param_set test_case;
} param_set;

void make_starting_point(itexture_plan * plan, int method)
{
	switch (method) {
		case 0:
			memset(plan->f_hat_iter, 0,
						 texture_get_omega_length(plan->mv) * sizeof(complex));
			break;
		default:
			printf("Illegal params.solver.starting_point_policy: %d\n", method);
			exit(-1);
	}
}

int is_regular_helper(int N1_new, int N2_new, int N,
											const iteration_param_set * it_par,
											const solver_param_set * sol_par, double *max_res,
											int clean_up)
{
	static int N1 = 0, N2 = 0;
	static double *h_phi = 0, *h_theta = 0, *r = 0;
	static complex *x = 0;

	if (!clean_up) {
		if (N <= 2) {
			return 1;
		} else {
			texture_plan plan, test_plan;
			itexture_plan iplan;
			complex *omega =
				(complex *) smart_malloc(texture_flat_length(N) * sizeof(complex));
			complex *omega_ref =
				(complex *) smart_malloc(texture_flat_length(N) * sizeof(complex));
			int success;

			if (N1 != N1_new || N2 != N2_new) {
				int h_phi_count, h_theta_count, r_phi_count, r_theta_count;
				grid_dim dims;

				if (N1_new > N1) {
					free(h_phi);
					free(h_theta);
					h_phi = (double *) smart_malloc(N1_new * sizeof(double));
					h_theta = (double *) smart_malloc(N1_new * sizeof(double));
				}
				if (N1_new * N2_new > N1 * N2) {
					free(r);
					free(x);
					r = (double *) smart_malloc(N1_new * N2_new * 2 * sizeof(double));
					x = (complex *) smart_malloc(N1_new * N2_new * sizeof(complex));
				}

				N1 = N1_new;
				N2 = N2_new;

				switch (sol_par->grid_id) {
					case 2:
					{
						dims.samples.N1 = N1;
						dims.samples.N2 = N2;
						calculate_grid(dims, h_phi, h_theta, r, sol_par->grid_id);
						break;
					}
					case 0:
					{
						split(N1, &h_phi_count, &h_theta_count);
						split(N2, &r_phi_count, &r_theta_count);
						dims.angles.h_phi_count = h_phi_count;
						dims.angles.h_theta_count = h_theta_count;
						dims.angles.r_phi_count = r_phi_count;
						dims.angles.r_theta_count = r_theta_count;
						calculate_grid(dims, h_phi, h_theta, r, sol_par->grid_id);
						break;
					}
					default:
						error("Illegal grid type!\n");
				}
			}

			init_omega(omega_ref, N, sol_par->omega_ref_policy);
			texture_precompute(N);
			texture_init(&plan, N, N1, N2, omega_ref, x, h_phi, h_theta, r);
			texture_init(&test_plan, N, N1, N2, omega, x, h_phi, h_theta, r);
			texture_trafo(&plan);
			texture_set_omega(&plan, omega);

			itexture_init_advanced(&iplan, &plan,
														 solver_flags(sol_par->solver_algorithm_id,
																					sol_par->weight_policy));

			make_starting_point(&iplan, sol_par->starting_point_policy);
			set_weights(&iplan, sol_par->weight_policy);
			memcpy(iplan.y, x, texture_get_x_length(&plan) * sizeof(complex));

			itexture_before_loop(&iplan);
			{
				double (*res_dist_arr[]) (const complex * vec, const complex * ref,
																	unsigned int length) = {
				l_2_rel_norm};
				double (*err_dist_arr[]) (const complex * vec, const complex * ref,
																	unsigned int length) = {
				l_2_rel_dist};
				double (*err_dist) (const complex * vec, const complex * ref,
														unsigned int length) =
					err_dist_arr[it_par->err_fun_id];
				double (*res_dist) (const complex * vec, const complex * ref,
														unsigned int length) =
					res_dist_arr[it_par->res_fun_id];
				double old_res, new_res;
#ifdef DEBUG_TRUE_RESIDUUM
				double true_res;
#endif
				int count = 0;

				new_res =
					res_dist(iplan.r_iter, iplan.y, texture_get_x_length(&plan));
				old_res = new_res;

				do {
					if (count % it_par->iterations_without_check == 0) {
						old_res = new_res;
					}

					itexture_loop_one_step(&iplan);
					count++;

					new_res =
						res_dist(iplan.r_iter, iplan.y, texture_get_x_length(&plan));
#ifdef DEBUG_RESIDUUM
					printf("residuum (%d): %lg\n", count, new_res);
					printf("error (%d):    %lg\n", count,
								 err_dist(iplan.f_hat_iter, omega_ref,
													texture_get_omega_length(&plan)));
					fflush(0);
#endif
#ifdef DEBUG_TRUE_RESIDUUM
					if (count % it_par->iterations_without_check == 0) {
						texture_set_omega(&test_plan, iplan.f_hat_iter);
						texture_trafo(&test_plan);
						true_res =
							err_dist(texture_get_x(&test_plan), iplan.y,
											 texture_get_x_length(&test_plan));
						fprintf(stderr, "residuum (true) (%d): %lg\n", count, true_res);
						fflush(0);
					}
#endif
				} while (new_res > it_par->res_delta
								 && !(count % it_par->iterations_without_check == 0 &&
											old_res > 0
											&& (old_res - new_res) / old_res <
											it_par->min_improve));

#ifdef DEBUG_CONVERGENCE_QUALITY
				if (new_res > it_par->res_delta) {
					char message[100];

					sprintf(message, "No convergence, residuum was %lg!", new_res);
					warning(message);
				}
#endif
				*max_res = NFFT_MAX(*max_res, new_res);

				success =
					(err_dist
					 (iplan.f_hat_iter, omega_ref,
						texture_get_omega_length(&plan)) <= it_par->err_delta);
			}

			itexture_finalize(&iplan);
			texture_finalize(&plan);
			// texture_finalize(&test_plan);
			free(omega);
			free(omega_ref);

			return success;
		}
	} else {
		free(h_phi);
		free(h_theta);
		free(r);
		free(x);
		return 0;
	}
}

int is_regular(int N1, int N2, int N, const iteration_param_set * it_par,
							 const solver_param_set * sol_par, double *max_res)
{
	return is_regular_helper(N1, N2, N, it_par, sol_par, max_res, 0);
}

int determine_max_N(int N1, int N2, int N_hint,
										const iteration_param_set * it_par,
										const solver_param_set * sol_par, double *max_res)
{
	int N = N_hint;
	if (is_regular(N1, N2, N, it_par, sol_par, max_res)) {
		do {
#ifdef DEBUG_N
			printf("regular: %d\n", N);
			fflush(0);
#endif
			N++;
		} while (is_regular(N1, N2, N, it_par, sol_par, max_res));
#ifdef DEBUG_N
		printf("not regular: %d\n", N);
		fflush(0);
#endif
		N--;
	} else {
		do {
#ifdef DEBUG_N
			printf("not regular: %d\n", N);
			fflush(0);
#endif
			N--;
		} while (!is_regular(N1, N2, N, it_par, sol_par, max_res));
#ifdef DEBUG_N
		printf("regular: %d\n", N);
		fflush(0);
#endif
	}

	return N;
}

void output_preliminaries(const param_set * params)
{
	const iteration_param_set *it_par = &params->iteration;
	const solver_param_set *sol_par = &params->solver;

	printf("# iteration parameters:\n");
	printf("# error function: %d (%s)\n", it_par->err_fun_id,
				 err_fun_id_descr[it_par->err_fun_id]);
	printf("# error_delta = %lg\n", it_par->err_delta);
	printf("# residuum function: %d (%s)\n", it_par->res_fun_id,
				 res_fun_id_descr[it_par->res_fun_id]);
	printf("# residuum_delta = %lg\n", it_par->res_delta);
	printf("# min_improve = %lg\n", it_par->min_improve);
	printf("# iterations_without_check = %d\n",
				 it_par->iterations_without_check);

	printf("# solver parameters:\n");
	printf("# grid_id: %d (%s)\n", sol_par->grid_id,
				 grid_descr[sol_par->grid_id]);
	printf("# grid_file_h: %s\n", sol_par->grid_file_h);
	printf("# grid_file_r: %s\n", sol_par->grid_file_r);
	printf("# starting_point_policy: %d (%s)\n",
				 sol_par->starting_point_policy,
				 starting_point_policy_descr[sol_par->starting_point_policy]);
	printf("# solver_algorithm_id: %d (%s)\n", sol_par->solver_algorithm_id,
				 solver_algo_descr[sol_par->solver_algorithm_id]);
	printf("# weight_policy: %d (%s)\n", sol_par->weight_policy,
				 weight_policy_descr[sol_par->weight_policy]);
	printf("# omega_ref_policy: %d (%s)\n", sol_par->omega_ref_policy,
				 omega_policy_descr[sol_par->omega_ref_policy]);
}

void output_N_curve_header(const param_set * params)
{
	if (params->test_case.output_policy != 0) {
		const test_case_param_set_0 *tc_par_0 = &params->test_case.c0;
		const test_case_param_set_1 *tc_par_1 = &params->test_case.c1;

		printf("# method: %d (%s)\n", params->test_case.method,
					 method_descr[params->test_case.method]);
		printf("# output_policy: %d (%s)\n", params->test_case.output_policy,
					 output_policy_descr[params->test_case.output_policy]);
		switch (params->test_case.method) {
			case 0:
				printf("# N1 = %d:%d:%d\n", tc_par_0->N1_start, tc_par_0->N1_incr,
							 tc_par_0->N1_end);
				printf("# N2 = %d:%d:%d\n", tc_par_0->N2_start, tc_par_0->N2_incr,
							 tc_par_0->N2_end);
				printf("# N_min_improve = %d\n", tc_par_0->N_min_improve);
				printf("# max_iter_without_improve = %d\n",
							 tc_par_0->max_iter_without_improve);
				break;
			case 1:
				printf("# N1 = %d:%d:%d\n", tc_par_1->N1_start, tc_par_1->N1_incr,
							 tc_par_1->N1_end);
				printf("# N2 = %d:%d:%d alpha = %d\n",
							 tc_par_1->N1_start * tc_par_1->alpha,
							 tc_par_1->N1_incr * tc_par_1->alpha,
							 tc_par_1->N1_end * tc_par_1->alpha, tc_par_1->alpha);
				printf("# N_min_improve = %d\n", tc_par_1->N_min_improve);
				printf("# max_iter_without_improve = %d\n",
							 tc_par_1->max_iter_without_improve);
				break;
			default:
				printf("Illegal params.test_case.method: %d\n",
							 params->test_case.method);
		}
	}
}

void output_N_curve_0_subheader(int N1, const param_set * params)
{
	printf("# N1 = %d\n", N1);
	fflush(0);
}

void output_N_curve_point(int x, int y, const param_set * params)
{
	printf("%d %d\n", x, y);
	fflush(0);
}

void determine_N_curve_0(const param_set * params)
{
#ifdef PARALLEL
	{
		int N1_num = (tc_par->N1_end - tc_par_N1_start) / tc_par->N1_incr + 1;
		int N2_num = (tc_par->N2_end - tc_par_N2_start) / tc_par->N2_incr + 1;
		int curve_size = N1_num * N2_num;
		int *N_val = (int *) smart_malloc(curve_size * sizeof(int));
		double *max_res = (double *) smart_malloc(curve_size * sizeof(double));
		int *calculated = (int *) smart_calloc(curve_size, sizeof(int));
		int N1_val = (int *) smart_malloc(curve_size * sizeof(int));
		int N2_val = (int *) smart_malloc(curve_size * sizeof(int));


		free(N_val);
		free(N1_val);
		free(N2_val);
		free(calculated);
		free(max_res);
	}
#else
	{
		const test_case_param_set_0 *tc_par = &params->test_case.c0;
		int N1, N2, old_N, new_N, iter_without_improve;
		double max_res = 0;

		output_N_curve_header(params);

		for (N1 = tc_par->N1_start; N1 <= tc_par->N1_end; N1 += tc_par->N1_incr) {
			output_N_curve_0_subheader(N1, params);

			new_N = -tc_par->N_min_improve;
			for (N2 = tc_par->N2_start, iter_without_improve = 0;
					 N2 <= tc_par->N2_end
					 && iter_without_improve <= tc_par->max_iter_without_improve;
					 N2 += tc_par->N2_incr) {
				old_N = new_N;
				new_N =
					determine_max_N(N1, N2, old_N, &params->iteration, &params->solver,
													&max_res);

				if (new_N - old_N < tc_par->N_min_improve) {
					iter_without_improve++;
				} else {
					iter_without_improve = 0;
				}

				output_N_curve_point(N2, new_N, params);
			}
		}
	}
#endif
}

void determine_N_curve_1(const param_set * params)
{
#ifdef PARALLEL
#else
	{
		const test_case_param_set_1 *tc_par = &params->test_case.c1;
		int N1, N2, old_N, new_N, iter_without_improve;
		double max_res = 0;

		output_N_curve_header(params);
		new_N = -tc_par->N_min_improve;
		for (N1 = tc_par->N1_start, iter_without_improve = 0;
				 N1 <= tc_par->N1_end
				 && iter_without_improve <= tc_par->max_iter_without_improve;
				 N1 += tc_par->N1_incr) {
			N2 = tc_par->alpha * N1;
			old_N = new_N;
			new_N =
				determine_max_N(N1, N2, old_N, &params->iteration, &params->solver,
												&max_res);

			if (new_N - old_N < tc_par->N_min_improve) {
				iter_without_improve++;
			} else {
				iter_without_improve = 0;
			}

			output_N_curve_point(N1, new_N, params);
		}
	}
#endif
}

void determine_N_curve(const param_set * params)
{
	const test_case_param_set *tc_par = &params->test_case;
	void (*determine_N_curve_funs[]) (const param_set * params) = {
	determine_N_curve_0, determine_N_curve_1};

	if (tc_par->method >= 0 && tc_par->method < 2) {
		determine_N_curve_funs[tc_par->method] (params);
	} else {
		printf("Illegal param: params.test_case.method = %d\n", tc_par->method);
		exit(-1);
	}
}

void init()
{
	unsigned short int seed1[] = { 1, 2, 3 };
	unsigned int seed2 = 1;

	seed48(seed1);
	srand(seed2);
}

void terminate()
{
	is_regular_helper(0, 0, 0, 0, 0, 0, 1);
	texture_forget();
}

void usage()
{
	// TODO
	printf("Some hints about the usage.");
}

void read_iteration_params(FILE * inp_file, iteration_param_set * param_set)
{
	fscanf(inp_file, "%d%d%lg%lg%lg%d", &param_set->err_fun_id,
				 &param_set->res_fun_id, &param_set->err_delta, &param_set->res_delta,
				 &param_set->min_improve, &param_set->iterations_without_check);
}

void read_solver_params(FILE * inp_file, solver_param_set * param_set)
{
	fscanf(inp_file, "%d\n%s\n%s\n%d%d%d%d", &param_set->grid_id,
				 param_set->grid_file_h, param_set->grid_file_r,
				 &param_set->starting_point_policy, &param_set->solver_algorithm_id,
				 &param_set->weight_policy, &param_set->omega_ref_policy);
}

void read_test_case_params(FILE * inp_file, test_case_param_set * param_set)
{
	fscanf(inp_file, "%d%d", &param_set->method, &param_set->output_policy);
	switch (param_set->method) {
		case 0:
			fscanf(inp_file, "%d%d%d%d%d%d%d%d", &param_set->c0.N1_start,
						 &param_set->c0.N1_end, &param_set->c0.N1_incr,
						 &param_set->c0.N2_start, &param_set->c0.N2_end,
						 &param_set->c0.N2_incr, &param_set->c0.N_min_improve,
						 &param_set->c0.max_iter_without_improve);
			break;
		case 1:
			fscanf(inp_file, "%d%d%d%d%d%d", &param_set->c1.alpha,
						 &param_set->c1.N1_start, &param_set->c1.N1_end,
						 &param_set->c1.N1_incr, &param_set->c1.N_min_improve,
						 &param_set->c1.max_iter_without_improve);
			break;
		default:
			printf("Illegal method for testcase!\n");
			exit(-1);
	}
}

/**
void make_grid_test_case()
{
	int i, j;
		solver_param_set sol_par;

		sol_par.grid_id = 0;

	{
		int N1 = 10;
		double h_phi[N1];
		double h_theta[N1];

		make_grid_h(h_phi, h_theta, N1, &sol_par);

		for (i = 0; i < N1; i++) {
			printf("%lg %lg\n", h_phi[i], h_theta[i]);
		}
		printf("\n");
	}

	{
		int N1 = 3, N2 = 10;
		double r[N1 * N2 * 2];

		make_grid_r(r, N1, N2, &sol_par);
		for (i = 0; i < N1; i++) {
			for (j = 0; j < N2; j++) {
				printf("%lg %lg\n", r[2 * (i * N2 + j)], r[2 * (i * N2 + j) + 1]);
			}
			printf("\n");
		}
	}

	exit(0);
}
*/

int main(int arglen, char *argv[])
{
	const char *test_cases = "test_cases.inp";
	const char *iteration_params = "iteration_parameters.inp";
	const char *solver_params = "solver_parameters.inp";

	// make_grid_test_case();

	if (arglen >= 2) {
		test_cases = argv[1];
		if (arglen >= 3) {
			solver_params = argv[2];
			if (arglen >= 4) {
				iteration_params = argv[3];
				if (arglen > 4) {
					printf("Too many parameters!");
					usage();
					return -1;
				}
			}
		}
	}

	init();

	{
		FILE *test_cases_file = fopen(test_cases, "r");
		FILE *iteration_params_file = fopen(iteration_params, "r");
		FILE *solver_params_file = fopen(solver_params, "r");
		param_set params;
		int test_cases_number;
		int count;

		read_iteration_params(iteration_params_file, &params.iteration);
		read_solver_params(solver_params_file, &params.solver);

		output_preliminaries(&params);

		fscanf(test_cases_file, "%d", &test_cases_number);
		for (count = 0; count < test_cases_number; count++) {
			read_test_case_params(test_cases_file, &params.test_case);
			determine_N_curve(&params);
		}

		fclose(test_cases_file);
		fclose(iteration_params_file);
		fclose(solver_params_file);
	}

	terminate();

	return 0;
}
