#include<nfft3.h>
#include<util.h>

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

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

const char *grid_id_descr[] = { "equidistant theta - equidistant phi" };
const char *starting_point_policy_descr[] = { "zero" };
const char *solver_algorithm_id_descr[] = { "CGNR" };
const char *weight_policy_descr[] = { "one" };
const char *omega_ref_policy_descr[] = { "rectangle" };

typedef struct solver_param_set_ {
	int grid_id;
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

void output_warning(const char *message)
{
	printf("Warning: %s\n", message);
	fflush(0);
}

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

inline double equidist(double start, double end, double i, double n,
											 int include_end)
{
	if (!include_end) {
		n++;
	}

	return start + (end - start) * ((double) i) / ((double) (n - 1));
}

void make_grid_h(double *h_phi, double *h_theta, int N1, int grid_id)
{
	switch (grid_id) {
		case 0:
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

			break;
		}
		default:
			printf("Illegal params.solver.grid_id: %d\n", grid_id);
			exit(-1);
	}
}

void make_grid_r(double *r, int N1, int N2, int grid_id)
{
	switch (grid_id) {
		case 0:
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
		default:
			printf("Illegal params.solver.grid_id: %d\n", grid_id);
			exit(-1);
	}
}

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

void make_weights(itexture_plan * plan, const solver_param_set * sol_par)
{
	switch (sol_par->weight_policy) {
		case 0:
			switch (sol_par->solver_algorithm_id) {
				case 0:
					return;
				default:
					printf("Illegal params.solver.solver_algorithm_id: %d\n",
								 sol_par->solver_algorithm_id);
					exit(-1);
			}
		default:
			printf("Illegal params.solver.weight_policy: %d\n",
						 sol_par->weight_policy);
			exit(-1);
	}
}

void make_omega_ref(complex * omega_ref, int length, int method)
{
	switch (method) {
		case 0:
		{
			int i;
			for (i = 0; i < length; i++) {
				omega_ref[i] = drand48() * rand() + I * drand48() * rand();
			}
			break;
		}
		default:
			printf("Illegal params.solver.omega_ref_policy: %d\n", method);
			exit(-1);
	}
}

int make_solver_flags(const solver_param_set * sol_par)
{
	const unsigned int algorithm_array[] = { CGNR };
	const unsigned int weight_array[] = { 0U };
	return algorithm_array[sol_par->
												 solver_algorithm_id] | weight_array[sol_par->
																														 weight_policy];
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

double dist_0(const complex * vec, const complex * ref, unsigned int length)
{
	double abs_dist = l_2_dist(vec, ref, length);
	double base = l_2_norm(ref, length);
	if (base != 0.0 || abs_dist != 0.0) {
		return abs_dist / base;
	} else {
		return 0.0;
	}
}

int is_regular_helper(int N1_new, int N2_new, int N,
											const iteration_param_set * it_par,
											const solver_param_set * sol_par, int clean_up)
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
				(complex *) malloc(texture_flat_length(N) * sizeof(complex));
			complex *omega_ref =
				(complex *) malloc(texture_flat_length(N) * sizeof(complex));
			int success;

			if (N1 != N1_new || N2 != N2_new) {
				if (N1_new > N1) {
					free(h_phi);
					free(h_theta);
					h_phi = (double *) malloc(N1_new * sizeof(double));
					h_theta = (double *) malloc(N1_new * sizeof(double));
				}
				if (N1_new * N2_new > N1 * N2) {
					free(r);
					free(x);
					r = (double *) malloc(N1_new * N2_new * 2 * sizeof(double));
					x = (complex *) malloc(N1_new * N2_new * sizeof(complex));
				}

				N1 = N1_new;
				N2 = N2_new;

				make_grid_h(h_phi, h_theta, N1, sol_par->grid_id);
				make_grid_r(r, N1, N2, sol_par->grid_id);
			}

			make_omega_ref(omega_ref, texture_flat_length(N),
										 sol_par->omega_ref_policy);
			texture_precompute(N);
			texture_init(&plan, N, N1, N2, omega_ref, x, h_phi, h_theta, r);
			texture_init(&test_plan, N, N1, N2, omega, x, h_phi, h_theta, r);
			texture_trafo(&plan);
			texture_set_omega(&plan, omega);

			itexture_init_advanced(&iplan, &plan, make_solver_flags(sol_par));

			make_starting_point(&iplan, sol_par->starting_point_policy);
			make_weights(&iplan, sol_par);
			memcpy(iplan.y, x, texture_get_x_length(&plan) * sizeof(complex));

			itexture_before_loop(&iplan);
			{
				double (*dist_arr[]) (const complex * vec, const complex * ref,
															unsigned int length) = {
				dist_0};
				double (*err_dist) (const complex * vec, const complex * ref,
														unsigned int length) =
					dist_arr[it_par->err_fun_id];
				double (*res_dist) (const complex * vec, const complex * ref,
														unsigned int length) =
					dist_arr[it_par->res_fun_id];
				double old_res, new_res;

				texture_set_omega(&test_plan, iplan.f_hat_iter);
				texture_trafo(&test_plan);
				new_res =
					res_dist(texture_get_x(&test_plan), iplan.y,
									 texture_get_x_length(&test_plan));

				do {
					int count;

					for (count = 0; count < it_par->iterations_without_check; count++) {
						itexture_loop_one_step(&iplan);
					}

					old_res = new_res;
					texture_set_omega(&test_plan, iplan.f_hat_iter);
					texture_trafo(&test_plan);
					new_res =
						res_dist(texture_get_x(&test_plan), iplan.y,
										 texture_get_x_length(&test_plan));
					//printf("residuum: %lg\n", new_res);
				} while (old_res > 0 && new_res > it_par->res_delta
								 && (old_res - new_res) / old_res >= it_par->min_improve);

				if (new_res > it_par->res_delta) {
					char message[100];

					sprintf(message, "No convergence, residuum was %lg!", new_res);
					output_warning(message);
				}

				success =
					(err_dist
					 (iplan.f_hat_iter, omega_ref,
						texture_get_omega_length(&plan)) <= it_par->err_delta);
			}

			itexture_finalize(&iplan);
			texture_finalize(&plan);
			texture_finalize(&test_plan);
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
							 const solver_param_set * sol_par)
{
	return is_regular_helper(N1, N2, N, it_par, sol_par, 0);
}

int determine_max_N(int N1, int N2, int N_hint,
										const iteration_param_set * it_par,
										const solver_param_set * sol_par)
{
	int N_max = MAX(1, N_hint);
	int N_min = N_max;
	if (is_regular(N1, N2, N_max, it_par, sol_par)) {
		do {
			printf("regular: %d\n", N_max);
			fflush(0);
			N_min = N_max;
			N_max *= 2;
		} while (is_regular(N1, N2, N_max, it_par, sol_par));
		printf("not regular: %d\n", N_max);
		fflush(0);
	} else {
		do {
			printf("not regular: %d\n", N_min);
			fflush(0);
			N_max = N_min;
			N_min /= 2;
		} while (!is_regular(N1, N2, N_min, it_par, sol_par));
		printf("regular: %d\n", N_min);
		fflush(0);
	}

	while (N_max - N_min > 1) {
		int N = (N_max + N_min) / 2;
		if (is_regular(N1, N2, N, it_par, sol_par)) {
			printf("regular: %d\n", N);
			fflush(0);
			N_min = N;
		} else {
			printf("not regular: %d\n", N);
			fflush(0);
			N_max = N;
		}
	}

	return MAX(2, N_min);
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
				 grid_id_descr[sol_par->grid_id]);
	printf("# starting_point_policy: %d (%s)\n",
				 sol_par->starting_point_policy,
				 starting_point_policy_descr[sol_par->starting_point_policy]);
	printf("# solver_algorithm_id: %d (%s)\n", sol_par->solver_algorithm_id,
				 solver_algorithm_id_descr[sol_par->solver_algorithm_id]);
	printf("# weight_policy: %d (%s)\n", sol_par->weight_policy,
				 weight_policy_descr[sol_par->weight_policy]);
	printf("# omega_ref_policy: %d (%s)\n", sol_par->omega_ref_policy,
				 omega_ref_policy_descr[sol_par->omega_ref_policy]);
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
	const test_case_param_set_0 *tc_par = &params->test_case.c0;
	int N1, N2, old_N, new_N, iter_without_improve;

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
				determine_max_N(N1, N2, old_N, &params->iteration, &params->solver);

			if (new_N - old_N < tc_par->N_min_improve) {
				iter_without_improve++;
			} else {
				iter_without_improve = 0;
			}

			output_N_curve_point(N2, new_N, params);
		}
	}
}

void determine_N_curve_1(const param_set * params)
{
	const test_case_param_set_1 *tc_par = &params->test_case.c1;
	int N1, N2, old_N, new_N, iter_without_improve;

	output_N_curve_header(params);
	new_N = -tc_par->N_min_improve;
	for (N1 = tc_par->N1_start, iter_without_improve = 0;
			 N1 <= tc_par->N1_end
			 && iter_without_improve <= tc_par->max_iter_without_improve;
			 N1 += tc_par->N1_incr) {
		N2 = tc_par->alpha * N1;
		old_N = new_N;
		new_N =
			determine_max_N(N1, N2, old_N, &params->iteration, &params->solver);

		if (new_N - old_N < tc_par->N_min_improve) {
			iter_without_improve++;
		} else {
			iter_without_improve = 0;
		}

		output_N_curve_point(N1, new_N, params);
	}
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
	is_regular_helper(0, 0, 0, 0, 0, 1);
	texture_forget();
}

void usage()
{
	//TODO
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
	fscanf(inp_file, "%d%d%d%d%d", &param_set->grid_id,
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

void make_grid_test_case()
{
	int i, j;

	{
		int N1 = 10;
		double h_phi[N1];
		double h_theta[N1];

		make_grid_h(h_phi, h_theta, N1, 0);

		for (i = 0; i < N1; i++) {
			printf("%lg %lg\n", h_phi[i], h_theta[i]);
		}
		printf("\n");
	}

	{
		int N1 = 3, N2 = 10;
		double r[N1 * N2 * 2];

		make_grid_r(r, N1, N2, 0);
		for (i = 0; i < N1; i++) {
			for (j = 0; j < N2; j++) {
				printf("%lg %lg\n", r[2 * (i * N2 + j)], r[2 * (i * N2 + j) + 1]);
			}
			printf("\n");
		}
	}

	exit(0);
}

int main(int arglen, char *argv[])
{
	const char *test_cases = "test_cases.inp";
	const char *iteration_params = "iteration_parameters.inp";
	const char *solver_params = "solver_parameters.inp";

	//make_grid_test_case();

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
