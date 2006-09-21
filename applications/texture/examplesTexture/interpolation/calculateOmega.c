#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <nfft3_texture.h>
#include <texture_util.h>

typedef struct properties_ {
	double res_delta;
	int iterations_without_check;
	double min_improve;
	int solver_algo;
	int weight_policy;
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
}

void read_properties(const char *propfile)
{
	FILE *f = fopen(propfile, "r");

	fscanf(f, "%lg%d%lg%d%d", &prop.res_delta, &prop.iterations_without_check,
				 &prop.min_improve, &prop.solver_algo, &prop.weight_policy);

	fclose(f);

	printf("# res_delta: %lg\n", prop.res_delta);
	printf("# iterations_without_check: %d\n", prop.iterations_without_check);
	printf("# min_improve: %lg\n", prop.min_improve);
	printf("# solver_algo: %d (%s)\n", prop.solver_algo,
				 solver_algo_descr[prop.solver_algo]);
	printf("# weight_policy: %d (%s)\n", prop.weight_policy,
				 weight_policy_descr[prop.weight_policy]);
	printf("#\n");
}

void calculate_omega()
{
	itexture_plan iplan;
	texture_plan plan, test_plan;
	unsigned int iflags;
	complex *test_omega;

	omega = (complex *) smart_malloc(texture_flat_length(N) * sizeof(complex));
	test_omega =
		(complex *) smart_malloc(texture_flat_length(N) * sizeof(complex));

	texture_precompute(N);
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);
	texture_init(&test_plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	iflags = solver_flags(prop.solver_algo, prop.weight_policy);

	itexture_init_advanced(&iplan, &plan, iflags);

	set_weights(&iplan, prop.weight_policy);

	memcpy(iplan.y, x, N1 * N2 * sizeof(complex));
	memset(iplan.f_hat_iter, 0, texture_flat_length(N) * sizeof(complex));
	itexture_before_loop(&iplan);
	{
		double new_res, old_res;
#ifdef DEBUG_TRUE_RESIDUUM
		double true_res;
#endif
		int count = 0;

		new_res =
			l_2_rel_norm(iplan.r_iter, iplan.y, texture_get_x_length(&plan));
		old_res = new_res;

		do {
			if (count % prop.iterations_without_check == 0) {
				old_res = new_res;
			}

			itexture_loop_one_step(&iplan);
			count++;

			new_res =
				l_2_rel_norm(iplan.r_iter, iplan.y, texture_get_x_length(&plan));

#ifdef DEBUG_RESIDUUM
			fprintf(stderr, "residuum (%d):        %lg\n", count, new_res);
			fflush(0);
#endif
#ifdef DEBUG_TRUE_RESIDUUM
			if (count % prop.iterations_without_check == 0) {
				texture_set_omega(&test_plan, iplan.f_hat_iter);
				texture_trafo(&test_plan);
				true_res =
					l_2_rel_dist(texture_get_x(&test_plan), iplan.y,
											 texture_get_x_length(&test_plan));
				fprintf(stderr, "residuum (true) (%d): %lg\n", count, true_res);
				fflush(0);
			}
#endif
		}
		while (new_res > prop.res_delta
					 && !(count % prop.iterations_without_check == 0
								&& old_res > 0
								&& (old_res - new_res) / old_res < prop.min_improve));

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


	free(test_omega);
	texture_forget();
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
	const char *propfile_name = "propfile_omega";
	const char *sample_file = "samples.in";
	const char *grid_h_file = "grid_h.in";
	const char *grid_r_file = "grid_r.in";
	FILE *f1, *f2, *f3;

	if (argc > 1) {
		sample_file = argv[1];
	}
	if (argc > 2) {
		grid_h_file = argv[2];
	}
	if (argc > 3) {
		grid_r_file = argv[3];
	}
	if (argc > 4) {
		propfile_name = argv[4];
	}

	if (argc <= 5) {
		init();
		read_N();

		printf("Omega\n");

		read_properties(propfile_name);

		f1 = fopen(grid_h_file, "r");
		f2 = fopen(grid_r_file, "r");
		f3 = fopen(sample_file, "r");
		read_samples(&N1, &N2, &h_phi, &h_theta, &r, &x, f1, f2, f3, stdout);
		fclose(f1);
		fclose(f2);
		fclose(f3);
		fflush(0);

		calculate_omega();

		write_omega(N, omega, stdout);

		cleanup();
	} else {
		usage();
	}

	return 0;
}
