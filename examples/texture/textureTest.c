#include <nfft3.h>
#include <util.h>
#include <nfsft.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void print_error(complex* vec1, complex* vec2, int n) {
	double infNorm = 0, twoNorm = 0;
	int i;

	for(i = 0; i < n; i++) {
		double absDiff = cabs(vec1[i] - vec2[i]); 
		if(absDiff > infNorm) {
			infNorm = absDiff;
		}

		twoNorm += absDiff * absDiff;
	}
	twoNorm = sqrt(twoNorm);

	printf("error in twoNorm: %g in infNorm: %g\n", twoNorm, infNorm);
}

inline complex expi(double phi) {
	return cos(phi) + I * sin(phi);
}

inline int equal(complex x, complex y, double delta) {
	return cabs(x - y) <= delta;
}

inline double equidist(double min_value, double max_value, int i, int n, 
		int incl) {
	if (!incl) {
		n++;
	}
	return min_value 
		+ ((double) (i)  * (max_value - min_value)) / (double) (n - 1);
}

inline double two_norm_dist(const complex *vec1, const complex *vec2, 
		unsigned int length) {
	double norm = 0;
	unsigned int i;

	for(i = 0; i < length; i++) {
		double x = cabs(vec1[i] - vec2[i]);

		norm += x*x;
	}
	norm = sqrt(norm);

	return norm;
}

inline double mse(const complex *vec1, const complex *vec2,
		unsigned int length) {
	return two_norm_dist(vec1, vec2, length) / length;
}

inline void initialize_angles(double *h_phi, double *h_theta, double *r, 
		int h_phi_count, int h_theta_count, int N2, int r_phi_count, 
		int r_theta_count) {
	int k = 0;
	int o = 0;
	int s, t, j;
	
	for(s = 0; s < h_phi_count; s++) {
		for(t = 0; t < h_theta_count; t++) {
			int i = s*h_theta_count + t;
			
			h_phi[i] = equidist(0, TEXTURE_MAX_ANGLE, s, h_phi_count, 0);
			h_theta[i] = equidist(0, TEXTURE_MAX_ANGLE / 2, t, h_theta_count, 1);

			for(j = 0; j < N2; j++) {
				r[2*(i*N2 + j)] = equidist(0, TEXTURE_MAX_ANGLE, k, r_phi_count, 0);
				r[2*(i*N2 + j) + 1] = 
					equidist(0, TEXTURE_MAX_ANGLE / 2, o, r_theta_count, 1);
			
				o++;
				if (o == r_theta_count) {
					o = 0;
					k = (k+1) % r_phi_count;
				}
			}
		}
	}
}

complex spherical_harmonic(int k, int n, double phi, double theta, int init) {
	static double p;
	static double p_old;
	static double old_theta;
	static int k0;
	static int n0;

	if (init || 
			abs(n) != n0 || k < k0 || 
			!equal(theta, old_theta, 1E-15)) {
		k0 = abs(n);
		old_theta = theta;
		p_old = 0;
		p = 1;

		for(n0 = 0; n0 < abs(n); n0++) {
			p *= sqrt((double) (2*n0 + 1) / (double) (2*n0 + 2)) 
						* sin(theta*2*PI/TEXTURE_MAX_ANGLE);
		}
	} 

	for(; k0 < k; k0++) {
		double p_new = 
			(double) (2*k0 + 1) / sqrt((double) ((k0 - n0 + 1) * (k0 + n0 + 1))) 
				* cos(theta*2*PI/TEXTURE_MAX_ANGLE) * p 
			- sqrt((double) ((k0 - n0) * (k0 + n0)) / 
						 (double) ((k0 - n0 + 1) * (k0 + n0 +	1))) 
				* p_old;
		p_old = p;
		p = p_new;
	}

	return p * expi(n*phi*2*PI/TEXTURE_MAX_ANGLE);
}

/*
void simple_solver_test() {
	texture_plan my_plan;
	itexture_plan my_iplan;
	int N = 2;
	int N1 = (texture_flat_index(N, N, N) + 1)/2 + 3;
	int N2 = 1;
	int i;
	complex *f_hat_ref; 

	printf("simple solver test:\n");
	
	texture_init(&my_plan, N, N1, N2);
	itexture_init(&my_iplan, &my_plan);

	f_hat_ref = (complex*) malloc(sizeof(complex) * my_plan.N_total);
	for(i = 0; i < my_plan.N_total; i++) {
		f_hat_ref[i] = i + 1;
		my_plan.f_hat[i] = f_hat_ref[i];
	}

	texture_trafo(&my_plan);
	vpr_complex(my_plan.f, my_plan.M_total, "f - the samples");
	vpr_complex(my_plan.f_hat, my_plan.N_total, 
			"f_hat - the corresponding frequencies");

	// destroying f_hat
	texture_adjoint(&my_plan);

	cp_complex(my_iplan.y, my_plan.f, my_plan.M_total);
	cp_complex(my_iplan.f_hat_iter, my_plan.f_hat, my_plan.N_total);

	itexture_before_loop(&my_iplan);
	for(i = 0; i < 4; i++)	{
		vpr_complex(my_iplan.f_hat_iter, my_plan.N_total, "f_hat_iter - the guess");

		SWAP_complex(my_iplan.f_hat_iter, my_plan.f_hat);
		texture_trafo(&my_plan);
		print_error(my_iplan.y, my_plan.f, my_plan.M_total);
		SWAP_complex(my_iplan.f_hat_iter, my_plan.f_hat);

		itexture_loop_one_step(&my_iplan);
	}

	free(f_hat_ref);
	itexture_finalize(&my_iplan);
	texture_finalize(&my_plan);
}
*/

void spherical_harmonic_test(const char *inp) {
	int N1, N2, N;
	int i, j, k, n0;
	char err_prefix[100];
	FILE *inp_file = fopen(inp, "r");

	sprintf(err_prefix, "spherical_harmonic_test failed (%s):\n", inp);
	fscanf(inp_file, "%d%d%d", &N1, &N2, &N);

	for(i = 1; i <= N1; i++) {
		double theta;
		
		fscanf(inp_file, "%lE", &theta);

		for(n0 = 0; n0 <= N; n0++) {
			for(k = n0; k <= N; k++) {
				for(j = 1; j <= N2; j++) {
					int sign;
					for (sign = -1; sign <= 1; sign += 2) {
						double phi, re_out, im_out;
						complex out;
						complex res;
						int n = n0 * sign;
						
						fscanf(inp_file, "%lE %lE %lE", &phi, &re_out, &im_out);
						out = re_out + I*im_out;
						
						res = spherical_harmonic(k, n, phi, theta, 0);
						if (!equal(out, res, 1E-13)) {
							printf("%s%d %d %g %g %g%+gi %g%+gi %g\n", 
								err_prefix,	k, n, phi, theta, creal(out), cimag(out), 
								creal(res), cimag(res),	cabs(out - res));
							return;
						}
					}
				}
			}
		}
	}
	fclose(inp_file);
}

//TODO Parameter
void unit_vector_test(const char *inp) {
	FILE *inp_file = fopen(inp, "r");
	int N, h_phi_count, h_theta_count, N1, N2, r_phi_count, r_theta_count;
	texture_plan plan;
	complex *omega, *x;
	double *r, *h_phi, *h_theta;
	int l, m0, n0;
	char err_prefix[100]; 
		
	sprintf(err_prefix, "unit_vector_test failed (%s):\n", inp);

	fscanf(inp_file, "%d%d%d%d%d%d", 
			&N, &h_phi_count, &h_theta_count, &N2, &r_phi_count, &r_theta_count);
	N1 = h_phi_count * h_theta_count;

	texture_precompute(N);
	
	omega = (complex*) calloc(texture_flat_length(N), sizeof(complex));
	x = (complex*) calloc(N1 * N2, sizeof(complex));
	h_phi = (double*) malloc(N1 * sizeof(double));
	h_theta = (double*) malloc(N1 * sizeof(double));
	r = (double*) malloc(N1 * N2 * 2 * sizeof(double));

	initialize_angles(h_phi, h_theta, r, h_phi_count, h_theta_count, N2, 
			r_phi_count, r_theta_count);
	
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	texture_set_r(&plan, r);
	texture_set_h_phi(&plan, h_phi);
	texture_set_h_theta(&plan, h_theta);

	for(m0 = 0; m0 <= N; m0++) {
		for(n0 = 0; n0 <= N; n0++) {
			for(l = MAX(m0, n0); l <= N; l++) {
				int sign;
				int m_sign[] = {-1, 1, -1, 1};
				int n_sign[] = {-1, -1, 1, 1};
				for(sign = 0; sign < 4; sign++) {
					int i, j;
					int m = m_sign[sign] * m0, n = n_sign[sign] * n0;
					
					omega[texture_flat_index(l, m, n)] = 1;
					texture_set_omega(&plan, omega);

					texture_trafo(&plan);

					for(i = 0; i < N1; i++) {
						for(j = 0; j < N2; j++) {
							complex out = 
								conj(
									spherical_harmonic(l, n, texture_get_h_phi(&plan)[i], 
										texture_get_h_theta(&plan)[i], 0))
								* spherical_harmonic(l, m, texture_get_r(&plan)[2*(i*N2 + j)], 
										texture_get_r(&plan)[2*(i*N2 + j) + 1], 0);
							complex res = texture_get_x(&plan)[i*N2 + j];

							if(!equal(out, res, 1E-10)) {
								printf(
									"%sl=%d m=%d n=%d h_phi=%g h_theta=%g r_phi=%g r_theta=%g\n", 
									err_prefix, l, m, n, texture_get_h_phi(&plan)[i], 
									texture_get_h_theta(&plan)[i], 
									texture_get_r(&plan)[2*(i*N2 + j)], 
									texture_get_r(&plan)[2*(i*N2 + j) + 1]);
								printf("result=%g%+gi reference=%g%+gi difference=%g\n",
									creal(res), cimag(res), creal(out), cimag(out), 
									cabs(out - res));
								return;
							}
						}
					}
					
					omega[texture_flat_index(l, m, n)] = 0;
					texture_set_omega(&plan, omega);
				}
			}
		}
	}

	texture_finalize(&plan);
	texture_forget();
	free(omega);
	free(x);
	free(r);
	free(h_phi);
	free(h_theta);
	fclose(inp_file);
}

void nfsft_test(const char *inp) {
	int N, theta_count, phi_count, N1;
	double threshold, tolerance;
	FILE *inp_file = fopen(inp, "r");
	nfsft_plan_old plan;
	int l, m, i, s, t;
	complex **f_hat;
	complex *f;
	double *x;
	char err_prefix[100];

	sprintf(err_prefix, "nfsft_test failed (%s):\n", inp);
	
	fscanf(inp_file, "%d%d%d%lg%lg", 
			&N, &phi_count, &theta_count, &threshold, &tolerance);
	N1 = phi_count * theta_count;

	nfsft_precompute_old(N, threshold, 0U);

	x = (double*) malloc(sizeof(double) * N1 * 2);
	i = 0;
	for(s = 0; s < phi_count; s++) {
		for(t = 0; t < theta_count; t++) {
			x[2*i] = -equidist(0, 1.0, s, phi_count, 0) + 1;
			x[2*i+1] = equidist(0, 0.5, t, theta_count, 1) + 2;
			i++;
		}
	}
	
	f = (complex*) malloc(sizeof(complex) * N1);
	
	f_hat = (complex**) malloc(sizeof(complex*) * (2*N + 1));
	for(m = 0; m <= 2*N; m++) {
		f_hat[m] = (complex*) malloc(sizeof(complex) * (next_power_of_2(N) + 1));
	}

	for(m = 0; m <= N; m++) {
		for(l = m; l <= N; l++) {
			int mm;
			//int ll;
			

			for(mm = 0; mm <= 2*N; mm++) {
				memset(f_hat[mm], 0, sizeof(complex) * (next_power_of_2(N) + 1));
			}
			
			f_hat[m+N][l] = 1;
			
			i = 0;
			for(s = 0; s < phi_count; s++) {
				for(t = 0; t < theta_count; t++) {
					double phi = -equidist(0, 1.0, s, phi_count, 0) + 1;
					double theta = equidist(0, 0.5, t, theta_count, 1) + 2;
					
					if(fabs(x[2*i] - phi) > 1E-15
						|| fabs(x[2*i+1] - theta) > 1E-15) {
						printf("%sl=%d m=%d i=%d\n",
							err_prefix, l, m, i);
						printf("plan corrupted:\n");
						printf("expected: phi=%lf theta=%lf\n",
							phi, theta);
						printf("given: phi=%lf theta=%lf\n",
							x[2*i], x[2*i+1]);
						printf("diff: phi=%lf theta=%lf\n",
							fabs(x[2*i]-phi), fabs(x[2*i+1]-theta));
					}
					i++;
				}
			}
			plan = nfsft_init_old(N, N1, f_hat, x, f, 0U);
			nfsft_trafo_old(plan);

/*			if(cabs(f_hat[m+N][l] - 1.0) > 1E-15)
			{
				printf("%sl=%d m=%d\n",
					err_prefix, l, m);
				printf("plan corrupted: expected=%lg%+lg received=%lg%+lg\n",
					1.0, 0.0, creal(f_hat[m+N][l]), cimag(f_hat[m+N][l]));
				return;
			}
					
			f_hat[m+N][l] = 0;
			for(mm = -N; mm <= N; mm++) {
				for(ll = abs(mm); ll <= N; ll++) {
					if(cabs(f_hat[mm+N][ll]) > 1E-15) {
						printf("%sl=%d m=%d\n",
							err_prefix, l, m);
						printf("plan corrupted: ll=%d mm=%d expected=0.0 given=%lg%+lg\n", 
							ll, mm, creal(f_hat[mm+N][ll]), cimag(f_hat[mm+N][ll]));
						return;
					}
				}
			}*/
			
			for(i = 0; i < N1; i++) {
				complex out = f[i];
				complex ref = spherical_harmonic(l, m, -x[2*i] * TEXTURE_MAX_ANGLE, 
												x[2*i+1] * TEXTURE_MAX_ANGLE, 0);
				if(!equal(ref, out, tolerance)) {
					printf("%sl=%d m=%d i=%d phi=%lg theta=%lg\n",
						err_prefix, l, m, i, x[2*i], x[2*i+1]);
					printf("out=%lg%+lgi ref=%lg%+lgi diff=%lg\n",
						creal(out), cimag(out), creal(ref), cimag(ref), cabs(ref - out));
					return;
				}
			}

//			printf("finish: %d %d\n", l, m);

			nfsft_finalize_old(plan);
		}
	}
	
	free(x);
	free(f);
	for(m = 0; m <= 2*N; m++) {
		free(f_hat[m]);
	}
	free(f_hat);
	
	fclose(inp_file);

	nfsft_forget_old();
}

void linearity_test(const char *inp) {
	FILE *inp_file = fopen(inp, "r");
	int N, h_phi_count, h_theta_count, N1, N2, r_phi_count, r_theta_count, w, it;
	double delta;
	char err_prefix[100];
	complex *omega, *x;
	complex *x_ref;
	double *h_phi, *h_theta, *r;
	texture_plan plan;
	int i, j, count, k;
	unsigned short int seed[] = {1, 2, 3};

	sprintf(err_prefix, "linearity_test failed (%s):\n", inp);

	fscanf(inp_file, "%d%d%d%d%d%d%d%d%lg", &N, &w, &h_phi_count, &h_theta_count, 
		&N2, &r_phi_count, &r_theta_count, &it, &delta);
	N1 = h_phi_count * h_theta_count;

	texture_precompute(N);
	seed48(seed);

	omega = (complex*) malloc(texture_flat_length(N) * sizeof(complex));
	x = (complex*) malloc(N1 * N2 * sizeof(complex));
	x_ref = (complex*) malloc(N1 * N2 * sizeof(complex));
	
	h_phi = (double*) malloc(N1 * sizeof(double));
	h_theta = (double*) malloc(N1 * sizeof(double));
	r = (double*) malloc(N1 * N2 * 2 * sizeof(double));
	initialize_angles(h_phi, h_theta, r, h_phi_count, h_theta_count, N2,
			r_phi_count, r_theta_count);
	
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	for(count = 0; count < it; count++) {
		double diff;
		
		memset(x_ref, 0, N1 * N2 * sizeof(complex));
		memset(omega, 0, texture_flat_length(N) * sizeof(complex));

		for(k = 0; k < w; k++) {
			int l, m, n;
			complex offset;
			
			l = rand() % (N + 1);
			m = (rand() % (2*l + 1)) - l;
			n = (rand() % (2*l + 1)) - l;
			offset = drand48() * rand() + I * drand48() * rand();

			omega[texture_flat_index(l, m, n)] += offset;

			for(i = 0; i < N1; i++) {
				for(j = 0; j < N2; j++) {
					x_ref[i*N2 + j] += offset 
						* conj(spherical_harmonic(l, n, texture_get_h_phi(&plan)[i], 
									 texture_get_h_theta(&plan)[i], 0))
						* spherical_harmonic(l, m, texture_get_r(&plan)[2*(i*N2 + j)], 
								texture_get_r(&plan)[2*(i*N2 + j) + 1], 0);
				}
			}
		}
		texture_set_omega(&plan, omega);

		texture_trafo(&plan);

		diff = mse(texture_get_x(&plan), x_ref, texture_get_x_length(&plan));
		if(diff	> delta) {
			printf("%scount=%d diff=%lg\n", err_prefix, count, diff);
		}	
	}
	
	texture_finalize(&plan);
	
	free(omega);
	free(x);
	free(x_ref);
	
	free(h_phi);
	free(h_theta);
	free(r);

	texture_forget();
}

void init() {
	spherical_harmonic(0, 0, 0, 0, 1);
}

void usage() {
	printf("A test for the texture transform.\n");
	printf("textureTest [Options]\n");
	printf("Options:\n");
	printf("-val : use small input values for usage with valgrind\n");
}

int main(int arglen, char *argv[]) {
	
	if (arglen == 2 && !strcmp(argv[1], "-val")) {
		// usage with valgrind
		
		printf("Test for usage with valgrind.\n");
		printf("If no further output is produced, everything seems to be okay.\n");

		init();
		
		//spherical_harmonic_test("spherical_harmonic_test.inp");

		//nfsft_test("nfsft_moderate_test.inp");

		//nfsft_test("nfsft_small_test.inp"); //reveals a bug in nfsft
		
		//unit_vector_test("unit_vector_moderate_test.inp");

		linearity_test("linearity_moderate_test.inp");
		
		//simple_solver_test();
	} else if (arglen == 1) {
		// default usage
		
		printf("Test with default parameters.\n");
		printf("If no further output is produced, everything seems to be okay.\n");

		init();
		
		//spherical_harmonic_test("spherical_harmonic_test.inp");

		//nfsft_test("nfsft_test.inp");

		//nfsft_test("nfsft_small_test.inp"); //reveals a bug in nfsft
		
		nfsft_test("nfsft_test2.inp");
		
		//unit_vector_test("unit_vector_test.inp");

		//linearity_test("linearity_test.inp");
		
		//simple_solver_test();
	} else {
		usage();
	}

	return 0;
}
