#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"nfft3_texture.h"
#include"util.h"
#include"texture_util.h"

/** @defgroup texture_regularity2 Texture: Regularity 2
 *
 * The program creates some random frequencies in omega_ref and evaluates them
 * for some given grid. From the samples it tries to reconstruct the 
 * frequencies using the solver.
 * This is repeated several times.
 * The output is some statistic about the size of the residuum and the distance
 * between the original frequencies omega_ref and the frequencies omega 
 * calculated by the solver.
 * 
 * @section CMA Command Line Arguments
 * There is only one argument: The prefix of the file to which output is 
 * printed. (default: prefix = "output/out3")
 *
 * @section Inp Input
 * The requested input is explained by messages on stderr.
 * - The grid has equidistant angles if it is not given by a file (f).
 * - grid_size determines the number of different @f$\theta@f$s and @f$\phi@f$s
 *   of the pole figures. The number of nodes is equal (e) or 16 times higher.
 *   (d).
 *
 * @section ProcOutp Processing and Output
 * The outputfile is of the form 
 * "<prefix>.omega<omega_type>.weight<weight_type>.<h_file>.<r_file>.<N>.<cases>".
 * Firstly, the program outputs all parameters to the outputfile.
 * 
 * The procedure afterwards is repeated for cases times.
 * The program creates some frequencies omega_ref of bandwidth N.
 * If and how much they are decreasing is determined by omega_type.
 * The samples x_ref are the result of a texture transform on omega_ref.
 * Then, the program uses the solver to reconstruct omega_ref from x_ref.
 * The starting point is 0.
 * Let omega be the frequencies of the current iteration, x the corresponding
 * samples, 
 * @f$\|x-x_{ref}\|_2 / \|x_{ref}\|_2 = res@f$ the current true relative 
 * residuum and @f$res_{min}@f$ the minimal true relative residuum that has 
 * occured until this iteration.
 * The solver stops if the updated squared residuum 
 * @f$\|x-x_{ref}\|_2^2 = dot\_r\_iter \leq 1e-20@f$,
 * if the true relative residuum 
 * @f$res \leq 1e-20@f$,
 * or if @f$(res_{min} - res) / res_{min} < min\_improve@f$ for more than 
 * max_iter_without_improve times.
 * Finally, @f$res_{min}@f$, the number min_iter of the iteration at which 
 * @f$res_{min}@f$ occured and the relative error 
 * @f$\|omega\_min - omega\_ref\|_2 / \|omega\_ref\|_2@f$ are saved.
 * The frequences omega_min are the frequencies of iteration min_iter.
 *
 * Finally, the residuums errors and iterations to achieve them is printed to
 * the output file in a histogram like manner.
 * Details are described by the comments of the output file.
 * 
 * @author Matthias Schmalz
 * @ingroup texture_examples
 */

// constants
const double min_improve = 0.01;
const int max_iter_without_improve = 10;
// unsigned short int dseed[3] = { 1, 2, 3 };
unsigned int iseed = 0;

// parameters influenced by the command line arguments

char *output_file_prefix = "output/out3";

// parameters given by the user
int N;
int grid_size;
char grid_type = 'e';
// char solver_type = 'n';
char h_file_name[100] = "dubna1_h";
char r_file_name[100] = "dubna1_r";
int weight_type = 1;
int omega_type = 3;
int cases = 100;

void read_params()
{
	int use_default;

	fprintf(stderr, "Use default values for most of the options (0 = no)? ");
	scanf("%d", &use_default);

	fprintf(stderr, "N: ");
	scanf("%d", &N);

	if (!use_default) {
		fprintf(stderr, "grid_type\n(e = all parameters are equal,\n");
		fprintf(stderr,
						"d = dense, i.e. about 16 times more nodes than pole figures,\n");
		fprintf(stderr, "f = file): ");
		scanf(" %c", &grid_type);

		if (grid_type == 'f') {
			fprintf(stderr, "file with pole figures: ");
			scanf(" %s", h_file_name);
			fprintf(stderr, "file with nodes: ");
			scanf(" %s", r_file_name);
		}
	}

	if (grid_type != 'f') {
		fprintf(stderr, "grid_size: ");
		scanf("%d", &grid_size);
	}

	if (!use_default) {
		// fprintf(stderr, "solver_type (n = normal, r = with regularization):
		// ");
		// scanf("%c", &solver_type);

		fprintf(stderr, "omega_type (0 = flat, 1 = 1/n, 2 = one, 3 = 1/n^2): ");
		scanf("%d", &omega_type);

		fprintf(stderr,
						"weight_type (0 = flat, 1 = linear, 2 = quadratic, 3 = cubic): ");
		scanf("%d", &weight_type);

		fprintf(stderr, "number of testcases: ");
		scanf("%d", &cases);
	}
}

// dependent parameters calculated from the user parameters
grid_dim gdim;
int N1, N2;
char output_path[500];
double *h_phi, *h_theta, *r;

void calculate_dependent_params()
{
	if (grid_type != 'f') {
		gdim.angles.h_phi_count = grid_size;
		gdim.angles.h_theta_count = grid_size;
		gdim.angles.r_phi_count = grid_size;
		gdim.angles.r_theta_count = grid_size;
		if (grid_type == 'd') {
			gdim.angles.r_phi_count *= 4;
			gdim.angles.r_theta_count *= 4;
		}

		N1 = gdim.angles.h_phi_count * (gdim.angles.h_theta_count - 2) + 2;
		N2 = gdim.angles.r_phi_count * (gdim.angles.r_theta_count - 2) + 2;

		h_phi = smart_malloc(N1 * sizeof(double));
		h_theta = smart_malloc(N1 * sizeof(double));
		r = smart_malloc(N1 * N2 * 2 * sizeof(double));

		calculate_grid(gdim, h_phi, h_theta, r, 0);

		sprintf(output_path, "%s.omega%d.weight%d.%c.%02d.%02d.%04d",
						output_file_prefix, omega_type, weight_type, grid_type, N,
						grid_size, cases);
	} else {
		FILE *h_file = fopen(h_file_name, "r");
		FILE *r_file = fopen(r_file_name, "r");

		read_grid(&N1, &N2, &h_phi, &h_theta, &r, h_file, r_file, stderr);

		sprintf(output_path, "%s.omega%d.weight%d.%s.%s.%02d.%04d",
						output_file_prefix, omega_type, weight_type, h_file_name,
						r_file_name, N, cases);

		fclose(h_file);
		fclose(r_file);
	}
}

texture_plan plan;
itexture_plan iplan;
complex *omega, *omega_ref, *omega_min, *x;
FILE *out;
double *residuums, *errors;
int *iterations;

void init()
{
	int sol_alg = ((N1 * N2) >= texture_flat_length(N)) ? 0 : 1;
	unsigned int iflags = solver_flags(sol_alg, weight_type);

	omega = smart_malloc(texture_flat_length(N) * sizeof(complex));
	omega_ref = smart_malloc(texture_flat_length(N) * sizeof(complex));
	omega_min = smart_malloc(texture_flat_length(N) * sizeof(complex));
	x = smart_malloc(N1 * N2 * sizeof(complex));

	residuums = smart_malloc(cases * sizeof(double));
	errors = smart_malloc(cases * sizeof(double));
	iterations = smart_malloc(cases * sizeof(int));

	texture_precompute(N);
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	itexture_init_advanced(&iplan, &plan, iflags);
	set_weights(&iplan, weight_type);

	out = fopen(output_path, "w");

	srand(iseed);
	// seed48(dseed);
}

void output_params()
{
	fprintf(out, "# omega: %d, weights: %d, testcases: %d\n", omega_type,
					weight_type, cases);
	fprintf(out, "# dim(J_N): %d, N: %d\n", texture_flat_length(N), N);

	if (grid_type == 'f') {
		fprintf(out, "# h: %s, r: %s\n", h_file_name, r_file_name);
	}

	fprintf(out, "# N1*N2: %d\n", N1 * N2);
	fprintf(out, "# N1: %d (phi: %d, theta: %d)\n", N1, gdim.angles.h_phi_count,
					gdim.angles.h_theta_count);
	fprintf(out, "# N2: %d (phi: %d, theta: %d)\n", N2, gdim.angles.r_phi_count,
					gdim.angles.r_theta_count);
	fflush(0);
}

void calculate_results()
{
	int count;

	for (count = 0; count < cases; count++) {
		double min_res = 1, res;
		int iter = 0;
		int iter_without_improve = 0;
		int min_iter = 0;

		init_omega(omega_ref, N, omega_type);
		texture_set_omega(&plan, omega_ref);
		texture_trafo(&plan);
		texture_set_omega(&plan, omega);

		memcpy(iplan.y, x, N1 * N2 * sizeof(complex));
		memset(iplan.f_hat_iter, 0, texture_flat_length(N) * sizeof(complex));
		itexture_before_loop(&iplan);

		memcpy(omega_min, iplan.f_hat_iter,
					 texture_flat_length(N) * sizeof(complex));
		do {
			itexture_loop_one_step(&iplan);
			iter++;

			texture_set_omega(&plan, iplan.f_hat_iter);
			texture_trafo(&plan);
			texture_set_omega(&plan, omega);
			res = l_2_rel_dist(x, iplan.y, N1 * N2);
			// fprintf(stderr, "residuum: %lg\n", res);

			if ((min_res - res) / min_res < min_improve) {
				iter_without_improve++;
			} else {
				iter_without_improve = 0;
			}

			if (res < min_res) {
				min_res = res;
				memcpy(omega_min, iplan.f_hat_iter,
							 texture_flat_length(N) * sizeof(complex));
				min_iter = iter;
			}

		} while (min_res > 1e-20 && iplan.dot_r_iter > 1e-20
						 && iter_without_improve <= max_iter_without_improve);

		errors[count] =
			l_2_rel_dist(omega_min, omega_ref, texture_flat_length(N));
		residuums[count] = min_res;
		iterations[count] = min_iter;

		fprintf(stderr, "case %d complete!\n", count + 1);
		fprintf(stderr, "res: %lg, err: %lg, iter: %d\n", residuums[count],
						errors[count], min_iter);

	}
}

inline int getBin(double x)
{
	return NFFT_MAX(ceil(log10(x)), -16);
}

void output_results()
{
	double minres = residuums[0], maxres = residuums[0];
	double minerr = 1000000000, maxerr = 0;
	int minbin, maxbin;
	int maxcountbin = -1000, maxcount = -1;
	int i;

	for (i = 0; i < cases; i++) {
		if (residuums[i] < minres) {
			minres = residuums[i];
		} else if (residuums[i] > maxres) {
			maxres = residuums[i];
		}
	}

	minbin = getBin(minres);
	maxbin = getBin(maxres);

	fprintf(out, "# upper bound of residuum\n");
	fprintf(out, "# number of cases in %%\n");
	fprintf(out, "# mean of iterations\n");
	for (i = minbin; i <= maxbin; i++) {
		fprintf(out, "10^%03d ", i);
	}
	fprintf(out, "\n");

	for (i = minbin; i <= maxbin; i++) {
		int count = 0;
		int j;

		for (j = 0; j < cases; j++) {
			if (i == getBin(residuums[j])) {
				count++;
			}
		}

		if (count > maxcount) {
			maxcount = count;
			maxcountbin = i;
		}

		fprintf(out, "%6.0f ", ((float) count) / cases * 100);
	}
	fprintf(out, "\n");

	for (i = minbin; i <= maxbin; i++) {
		int count = 0;
		int iter = 0;
		int j;

		for (j = 0; j < cases; j++) {
			if (i == getBin(residuums[j])) {
				count++;
				iter += iterations[j];
			}
		}

		fprintf(out, "%6.0f ", ((float) iter) / count);
	}
	fprintf(out, "\n\n");

	fprintf(out,
					"# For the class 10^%03d we give the distribution of errors:\n",
					maxcountbin);

	for (i = 0; i < cases; i++) {
		if (getBin(residuums[i]) == maxcountbin) {
			if (errors[i] < minerr) {
				minerr = errors[i];
			}
			if (errors[i] > maxerr) {
				maxerr = errors[i];
			}
		}
	}

	minbin = getBin(minerr);
	maxbin = getBin(maxerr);

	fprintf(out, "# upper bound of error\n");
	fprintf(out, "# number of cases in %%\n");
	for (i = minbin; i <= maxbin; i++) {
		fprintf(out, "10^%03d ", i);
	}
	fprintf(out, "\n");

	for (i = minbin; i <= maxbin; i++) {
		int count = 0;
		int j;

		for (j = 0; j < cases; j++) {
			if (getBin(residuums[j]) == maxcountbin && getBin(errors[j]) == i) {
				count++;
			}
		}

		fprintf(out, "%6.0f ", ((float) count) / maxcount * 100);
	}
	fprintf(out, "\n");
}

void cleanup()
{
	fclose(out);

	itexture_finalize(&iplan);
	texture_finalize(&plan);
	texture_forget();

	free(omega);
	free(omega_ref);
	free(omega_min);
	free(x);
	free(h_phi);
	free(h_theta);
	free(r);
	free(residuums);
	free(errors);
	free(iterations);
}

int main(int argc, char *argv[])
{
	if (argc >= 2) {
		output_file_prefix = argv[1];
	}

	read_params();
	calculate_dependent_params();
	init();
	output_params();
	calculate_results();
	output_results();
	cleanup();

	return 0;
}
