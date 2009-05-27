#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <nfft3util.h>
#include <nfft3_texture.h>
#include <texture_nfft3util.h>

/** @defgroup texture_correctness Texture: Correctness
 * This program was designed to check the correctness of the implementation
 * of the texture transforms.
 * 
 * @author Matthias Schmalz
 * @ingroup texture_examples
 * @{
 */

/** Returns if @f$ \| vec - ref \|_2 \leq delta \cdot \|ref\|_2 @f$
 *
 * @par length - the length of vec and ref
 */
inline int equal_two_norm_rel(const double _Complex * vec,
															const double _Complex * ref, unsigned int length,
															double delta)
{
	return l_2_dist(vec, ref, length) <= delta * l_2_norm(ref, length);
}

/** Initializes the nodes.
 * The nodes are taken from a
 * cartesian rectangle of h_phi_count / r_phi_count equidistant angles for the 
 * latitudes and 
 * h_theta_count / r_theta_count
 * equidistant angles for the longitudes.
 * 
 * If @f$ h\_phi\_count \cdot h\_theta\_count \cdot N2 
 * \leq r\_phi\_count \cdot r\_theta\_count @f$, i.e. if the rectangle is too
 * small, the nodes will be extended periodically.
 *
 * @par h_phi, h_theta store the pole figures.
 * @par r stores the nodes of the pole figures.
 * @par h_phi_count, h_theta_count determines the number of different 
 * @par                            latitudes / longitudes that will be taken for
 * @par                            the pole figures.
 * @par r_phi_count, r_theta_count determines the number of different
 * @par                            latitudes / longitudes that will be taken for
 * @par                            the nodes in the pole figures.
 * @par N2 - the number of nodes per pole figure
 *
 * @pre The size of h_phi and h_theta has to be 
 * @f$ h\_phi\_count \cdot h\_theta\_count @f$.
 * The size of r has to be N2.
 */
inline void initialize_angles(double *h_phi, double *h_theta, double *r,
															int h_phi_count, int h_theta_count, int N2,
															int r_phi_count, int r_theta_count)
{
	int k = 0;
	int o = 0;
	int s, t, j;

	for (s = 0; s < h_phi_count; s++) {
		for (t = 0; t < h_theta_count; t++) {
			int i = s * h_theta_count + t;

			h_phi[i] = equidist(0, TEXTURE_MAX_ANGLE, s, h_phi_count, 0);
			h_theta[i] = equidist(0, TEXTURE_MAX_ANGLE / 2, t, h_theta_count, 1);

			for (j = 0; j < N2; j++) {
				r[2 * (i * N2 + j)] = equidist(-TEXTURE_MAX_ANGLE / 2,
																			 TEXTURE_MAX_ANGLE / 2, k, r_phi_count,
																			 0);
				r[2 * (i * N2 + j) + 1] =
					equidist(0, TEXTURE_MAX_ANGLE / 2, o, r_theta_count, 1);

				o++;
				if (o == r_theta_count) {
					o = 0;
					k = (k + 1) % r_phi_count;
				}
			}
		}
	}
}

/** A test for the inverse texture transform.
 * Runs the inverse transformation for a given number of iterations.
 * Afterwards, if the residuum is not sufficently small, prints an error 
 * message containing the relative error (see ::l_2_dist) "diff" and the
 * absolute error goal "ref".
 *
 * The test parameters are read from a parameter file.
 * The parameters in the file are (in this order):
 * - the bandwidth N
 * - the number of longitudes of the pole figures
 * - the number of latitudes of the pole figures
 * - the number of nodes per pole figure
 * - the number of longitudes of nodes
 * - the number of latitudes of nodes
 * - the maximum number of iterations
 * - the relative error goal
 *
 * @par inp - the name of the parameter file.
 */
void simple_solver_test(const char *inp)
{
	FILE *inp_file = fopen(inp, "r");
	texture_plan plan;
	itexture_plan iplan;
	int N, N1, N2, r_phi_count, r_theta_count, h_phi_count, h_theta_count,
		max_iter;
	double delta;
	int i;
	double _Complex *omega, *x, *omega_ref;
	double *h_phi, *h_theta, *r;
	char err_prefix[100];
	// unsigned short int seed[] = { 1, 2, 3 };

	// Output user information.
	printf("*** simple_solver_test (%s)\n", inp);
	sprintf(err_prefix, "simple_solver_test failed (%s):\n", inp);

	// Read parameters.
	fscanf(inp_file, "%d%d%d%d%d%d%d%lg", &N, &h_phi_count, &h_theta_count,
				 &N2, &r_phi_count, &r_theta_count, &max_iter, &delta);
	N1 = h_phi_count * h_theta_count;
	// seed48(seed);
	srand(0);

	// Prepare input parameters.
	texture_precompute(N);
	omega =
		(double _Complex *) smart_malloc(texture_flat_length(N) *
																		sizeof(double _Complex));
	x = (double _Complex *) smart_malloc(N1 * N2 * sizeof(double _Complex));
	omega_ref =
		(double _Complex *) smart_malloc(texture_flat_length(N) *
																		sizeof(double _Complex));
	h_phi = (double *) smart_malloc(N1 * sizeof(double));
	h_theta = (double *) smart_malloc(N2 * sizeof(double));
	r = (double *) smart_malloc(N1 * N2 * 2 * sizeof(double));

	initialize_angles(h_phi, h_theta, r, h_phi_count, h_theta_count, N2,
										r_phi_count, r_theta_count);
	for (i = 0; i < texture_flat_length(N); i++) {
		omega_ref[i] = crand();
	}

	// Prepare plans.
	texture_init(&plan, N, N1, N2, omega_ref, x, h_phi, h_theta, r);
	texture_trafo(&plan);
	texture_set_omega(&plan, omega);

	itexture_init(&iplan, &plan);

	memset(iplan.f_hat_iter, 0,
				 texture_get_omega_length(&plan) * sizeof(double _Complex));
	memcpy(iplan.y, x, texture_get_x_length(&plan) * sizeof(double _Complex));

	// Run the inverse transformation.
	itexture_before_loop(&iplan);
	for (i = 0;
			 i < max_iter &&
			 !equal_two_norm_rel(iplan.f_hat_iter, omega_ref,
													 texture_get_omega_length(&plan), delta); i++) {
		itexture_loop_one_step(&iplan);
	}

	// If the residuum is not small, print an error message.
	if (!equal_two_norm_rel(iplan.f_hat_iter, omega_ref,
													texture_get_omega_length(&plan), delta)) {
		printf("%sdiff=%lg ref=%lg\n",
					 err_prefix,
					 l_2_dist(iplan.f_hat_iter, omega_ref,
										texture_get_omega_length(&plan)),
					 delta * l_2_norm(omega_ref, texture_get_omega_length(&plan)));
	}
	// Clean up.
	itexture_finalize(&iplan);
	texture_finalize(&plan);

	texture_forget();

	free(omega);
	free(x);
	free(omega_ref);
	free(h_phi);
	free(h_theta);
	free(r);

	fclose(inp_file);
}

/** A test for ::spherical_harmonic.
 * Computes the spherical harmonic with a given set of parameters and compares
 * them with given reference values.
 * In case of error the execution will be interrupted and there will be an 
 * error message containing degree, order,
 * longitude, latitude, reference value, computed value and absolute difference
 * (in this order).
 *
 * The parameter file contains (in this order):
 * - the number of pole figures N1
 * - the number of nodes per pole figure N2
 * - the bandwidth N
 * - the maximum absolute error
 * - N1 times:
 *   - the latitude
 *   - N2 times:
 *     - the longitude
 *     - the real part of the reference value
 *     - the imaginary part of the reference value
 *
 * @par inp - the name of the parameter file
 */
void spherical_harmonic_test(const char *inp)
{
	int N1, N2, N;
	double delta;
	int i, j, k, n0;
	char err_prefix[100];
	FILE *inp_file = fopen(inp, "r");

	// Print user information.
	sprintf(err_prefix, "spherical_harmonic_test failed (%s):\n", inp);
	printf("*** spherical_harmonic_test (%s)\n", inp);

	// Read parameters from the input file.
	fscanf(inp_file, "%d%d%d%lf", &N1, &N2, &N, &delta);

	for (i = 1; i <= N1; i++) {
		double theta;

		fscanf(inp_file, "%lE", &theta);

		for (n0 = 0; n0 <= N; n0++) {
			for (k = n0; k <= N; k++) {
				for (j = 1; j <= N2; j++) {
					int sign;
					for (sign = -1; sign <= 1; sign += 2) {
						double phi, re_out, im_out;
						double _Complex out;
						double _Complex res;
						int n = n0 * sign;

						fscanf(inp_file, "%lE %lE %lE", &phi, &re_out, &im_out);
						out = re_out + _Complex_I * im_out;

						res = spherical_harmonic(k, n, phi, theta);
						if (!equal(out, res, delta)) {
							printf("%s%d %d %g %g %g%+gi %g%+gi %g\n",
										 err_prefix, k, n, phi, theta, creal(out), cimag(out),
										 creal(res), cimag(res), cabs(out - res));
							return;
						}
					}
				}
			}
		}
	}
	fclose(inp_file);
}

/** A test of the texture transform.
 * The texture transform will be executed on unit vectors.
 * The output is compared with the results produced by ::spherical_harmonic.
 * In case of error the execution will be interrupted and an error message 
 * will be printed containing
 * - the indices l, m, n of the unit vector,
 * - the pole figure and node of the erroneos intensity,
 * - the intensity calculated by the texture transform,
 * - the intensity calculated by ::spherical_harmonic and
 * - the absolute difference.
 * 
 * The parameter file contains (in this order):
 * - the bandwidth N,
 * - the number of different longitudes and latitudes of the pole figures,
 * - the number of nodes,
 * - the number of different longitudes and latitudes of the nodes and
 * - the maximum difference between intensities that is not an error.
 * 
 * @par inp - the filename of the parameter file.
 */
void unit_vector_test(const char *inp)
{
	FILE *inp_file = fopen(inp, "r");
	int N, h_phi_count, h_theta_count, N1, N2, r_phi_count, r_theta_count;
	double delta;
	texture_plan plan;
	double _Complex *omega, *x;
	double *r, *h_phi, *h_theta;
	int l, m0, n0;
	char err_prefix[100];

	sprintf(err_prefix, "unit_vector_test failed (%s):\n", inp);
	printf("*** unit_vector_test (%s)\n", inp);

	fscanf(inp_file, "%d%d%d%d%d%d%lf",
				 &N, &h_phi_count, &h_theta_count, &N2, &r_phi_count, &r_theta_count,
				 &delta);
	N1 = h_phi_count * h_theta_count;

	texture_precompute(N);

	omega =
		(double _Complex *) smart_calloc(texture_flat_length(N),
																		sizeof(double _Complex));
	x = (double _Complex *) smart_malloc(N1 * N2 * sizeof(double _Complex));
	h_phi = (double *) smart_malloc(N1 * sizeof(double));
	h_theta = (double *) smart_malloc(N1 * sizeof(double));
	r = (double *) smart_malloc(N1 * N2 * 2 * sizeof(double));

	initialize_angles(h_phi, h_theta, r, h_phi_count, h_theta_count, N2,
										r_phi_count, r_theta_count);

	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	for (m0 = 0; m0 <= N; m0++) {
		for (n0 = 0; n0 <= N; n0++) {
			for (l = NFFT_MAX(m0, n0); l <= N; l++) {
				int sign;
				int m_sign[] = { -1, 1, -1, 1 };
				int n_sign[] = { -1, -1, 1, 1 };
				for (sign = 0; sign < 4; sign++) {
					int i, j;
					int m = m_sign[sign] * m0, n = n_sign[sign] * n0;

					omega[texture_flat_index(l, m, n)] = 1;
					texture_set_omega(&plan, omega);

					texture_trafo(&plan);

					for (i = 0; i < N1; i++) {
						for (j = 0; j < N2; j++) {
							double _Complex out =
								conj(spherical_harmonic(l, n, texture_get_h_phi(&plan)[i],
																				texture_get_h_theta(&plan)[i]))
								* spherical_harmonic(l, m,
																		 texture_get_r(&plan)[2 * (i * N2 + j)],
																		 texture_get_r(&plan)[2 * (i * N2 + j) +
																													1]);
							double _Complex res = texture_get_x(&plan)[i * N2 + j];

							if (!equal(out, res, delta)) {
								printf
									("%sl=%d m=%d n=%d h_phi=%g h_theta=%g r_phi=%g r_theta=%g\n",
									 err_prefix, l, m, n, texture_get_h_phi(&plan)[i],
									 texture_get_h_theta(&plan)[i],
									 texture_get_r(&plan)[2 * (i * N2 + j)],
									 texture_get_r(&plan)[2 * (i * N2 + j) + 1]);
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

/** A test of the nfsft.
 * The nfsft will be executed on unit vectors.
 * The output is compared with the results produced by ::spherical_harmonic.
 * In case of an erroneos sample value the execution will be interrupted and 
 * an error message 
 * will be printed containing
 * - the indices l, m of the unit vector,
 * - the index i of the node
 * - the node itself,
 * - the intensity calculated by the nfsft,
 * - the intensity calculated by ::spherical_harmonic and
 * - the absolute difference.
 *
 * If the nfsft destroys the nodes in the plan, an error message will be printed
 * containing
 * - the indices l, m of the unit vector,
 * - the index i of the node,
 * - the expected node,
 * - the node stored in the plan and
 * - the absolute differences of angles.
 *
 * In both cases the execution will be interrupted.
 * 
 * The parameter file contains (in this order):
 * - the bandwidth N,
 * - the number of different longitudes and latitudes of the nodes,
 * - the threshold and
 * - the maximum difference between samples that is not an error.
 * 
 * @par inp - the filename of the parameter file.
 */
void nfsft_test(const char *inp)
{
	int N, theta_count, phi_count, N1;
	double threshold, tolerance;
	FILE *inp_file = fopen(inp, "r");
	nfsft_plan plan;
	int l, m, i, s, t;
	double _Complex *f_hat;
	double _Complex *f;
	double *x;
	char err_prefix[100];

	sprintf(err_prefix, "nfsft_test failed (%s):\n", inp);
	printf("*** nfsft_test (%s)\n", inp);

	fscanf(inp_file, "%d%d%d%lg%lg",
				 &N, &phi_count, &theta_count, &threshold, &tolerance);
	N1 = phi_count * theta_count;

	nfsft_precompute(N, threshold, 0U, 0U);

	x = (double *) smart_malloc(sizeof(double) * N1 * 2);
	i = 0;
	for (s = 0; s < phi_count; s++) {
		for (t = 0; t < theta_count; t++) {
			x[2 * i] = equidist(-0.5, 0.5, s, phi_count, 0);
			x[2 * i + 1] = equidist(0, 0.5, t, theta_count, 1);
			i++;
		}
	}

	f = (double _Complex *) smart_malloc(sizeof(double _Complex) * N1);

	f_hat = smart_malloc(sizeof(double _Complex) * NFSFT_F_HAT_SIZE(N));

	nfsft_init_guru(&plan, N, N1, TEXTURE_DEF_NFSFT_INIT_FLAGS,
									TEXTURE_DEF_NFFT_INIT_FLAGS, TEXTURE_DEF_NFFT_CUTOFF);

	for (m = 0; m <= N; m++) {
		for (l = m; l <= N; l++) {
			// int mm;
			// int ll;

			memset(f_hat, 0, sizeof(double _Complex) * NFSFT_F_HAT_SIZE(N));

			f_hat[NFSFT_INDEX(l, m, &plan)] = 1;

			i = 0;
			for (s = 0; s < phi_count; s++) {
				for (t = 0; t < theta_count; t++) {
					double phi = equidist(-0.5, 0.5, s, phi_count, 0);
					double theta = equidist(0, 0.5, t, theta_count, 1);

					if (fabs(x[2 * i] - phi) > 1E-15
							|| fabs(x[2 * i + 1] - theta) > 1E-15) {
						printf("%sl=%d m=%d i=%d\n", err_prefix, l, m, i);
						printf("plan corrupted:\n");
						printf("expected: phi=%lf theta=%lf\n", phi, theta);
						printf("given: phi=%lf theta=%lf\n", x[2 * i], x[2 * i + 1]);
						printf("diff: phi=%lf theta=%lf\n",
									 fabs(x[2 * i] - phi), fabs(x[2 * i + 1] - theta));
						return;
					}
					i++;
				}
			}
			plan.f_hat = f_hat;
			plan.x = x;
			plan.f = f;
			nfsft_precompute_x(&plan);
			nfsft_trafo(&plan);

			/* if(cabs(f_hat[m+N][l] - 1.0) > 1E-15) { printf("%sl=%d m=%d\n",
			   err_prefix, l, m); printf("plan corrupted: expected=%lg%+lg
			   received=%lg%+lg\n", 1.0, 0.0, creal(f_hat[m+N][l]),
			   cimag(f_hat[m+N][l])); }

			   f_hat[m+N][l] = 0; for(mm = -N; mm <= N; mm++) { for(ll = abs(mm);
			   ll <= N; ll++) { if(cabs(f_hat[mm+N][ll]) > 1E-15) { printf("%sl=%d
			   m=%d\n", err_prefix, l, m); printf("plan corrupted: ll=%d mm=%d
			   expected=0.0 given=%lg%+lg\n", ll, mm, creal(f_hat[mm+N][ll]),
			   cimag(f_hat[mm+N][ll])); } } } */

			for (i = 0; i < N1; i++) {
				double _Complex out = f[i];
				double _Complex ref =
					spherical_harmonic(l, m, x[2 * i] * TEXTURE_MAX_ANGLE,
														 x[2 * i + 1] * TEXTURE_MAX_ANGLE);
				if (!equal(ref, out, tolerance)) {
					printf("%sl=%d m=%d i=%d phi=%lg theta=%lg\n",
								 err_prefix, l, m, i, x[2 * i], x[2 * i + 1]);
					printf("out=%lg%+lgi ref=%lg%+lgi diff=%lg\n",
								 creal(out), cimag(out), creal(ref), cimag(ref),
								 cabs(ref - out));
					return;
				}
			}

			// printf("finish: %d %d\n", l, m);

		}
	}

	nfsft_finalize(&plan);
	free(x);
	free(f);
	free(f_hat);

	fclose(inp_file);

	nfsft_forget();
}

/** A test of the adjoint texture transform.
 * The adjoint texture transform will be executed on unit vectors.
 * The output is compared with the results produced by ::spherical_harmonic.
 * In case of error the execution will be interrupted and an error message 
 * will be printed containing
 * - the indices l, m, n of the erroneous frequency,
 * - the pole figure and node of the intensity beeing one,
 * - the frequency calculated by the texture transform,
 * - the frequency calculated by ::spherical_harmonic and
 * - the absolute difference.
 * 
 * The parameter file contains (in this order):
 * - the bandwidth N,
 * - the number of different longitudes and latitudes of the pole figures,
 * - the number of nodes,
 * - the number of different longitudes and latitudes of the nodes and
 * - the maximum difference between frequencies that is not an error.
 * 
 * @par inp - the filename of the parameter file.
 */
void unit_vector_adjoint_test(const char *inp)
{
	FILE *inp_file = fopen(inp, "r");
	int N, h_phi_count, h_theta_count, N1, N2, r_phi_count, r_theta_count;
	double delta;
	texture_plan plan;
	double _Complex *omega, *x;
	double *r, *h_phi, *h_theta;
	int i, j;
	char err_prefix[100];

	sprintf(err_prefix, "unit_vector_adjoint_test failed (%s):\n", inp);
	printf("*** unit_vector_adjoint_test (%s)\n", inp);

	fscanf(inp_file, "%d%d%d%d%d%d%lf",
				 &N, &h_phi_count, &h_theta_count, &N2, &r_phi_count, &r_theta_count,
				 &delta);
	N1 = h_phi_count * h_theta_count;

	texture_precompute(N);

	omega =
		(double _Complex *) smart_malloc(texture_flat_length(N) *
																		sizeof(double _Complex));
	x = (double _Complex *) smart_calloc(N1 * N2, sizeof(double _Complex));
	h_phi = (double *) smart_malloc(N1 * sizeof(double));
	h_theta = (double *) smart_malloc(N1 * sizeof(double));
	r = (double *) smart_malloc(N1 * N2 * 2 * sizeof(double));

	initialize_angles(h_phi, h_theta, r, h_phi_count, h_theta_count, N2,
										r_phi_count, r_theta_count);

	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			int l, m0, n0;

			x[i * N2 + j] = 1;
			texture_set_x(&plan, x);

			texture_adjoint(&plan);

			for (m0 = 0; m0 <= N; m0++) {
				for (n0 = 0; n0 <= N; n0++) {
					for (l = NFFT_MAX(m0, n0); l <= N; l++) {
						int sign;
						int m_sign[] = { 1, 1, -1, -1 };
						int n_sign[] = { 1, -1, 1, -1 };

						for (sign = 0; sign < 4; sign++) {
							int m = m0 * m_sign[sign];
							int n = n0 * n_sign[sign];

							double _Complex out =
								spherical_harmonic(l, n, texture_get_h_phi(&plan)[i],
																	 texture_get_h_theta(&plan)[i])
								*
								conj(spherical_harmonic
										 (l, m, texture_get_r(&plan)[2 * (i * N2 + j)],
											texture_get_r(&plan)[2 * (i * N2 + j) + 1]));
							double _Complex res
								= texture_get_omega(&plan)[texture_flat_index(l, m, n)];

							if (!equal(out, res, delta)) {
								printf
									("%sl=%d m=%d n=%d h_phi=%g h_theta=%g r_phi=%g r_theta=%g\n",
									 err_prefix, l, m, n, texture_get_h_phi(&plan)[i],
									 texture_get_h_theta(&plan)[i],
									 texture_get_r(&plan)[2 * (i * N2 + j)],
									 texture_get_r(&plan)[2 * (i * N2 + j) + 1]);
								printf("result=%g%+gi reference=%g%+gi difference=%g\n",
											 creal(res), cimag(res), creal(out), cimag(out),
											 cabs(out - res));
								return;
							}
						}
					}
				}
			}

			x[i * N2 + j] = 0;
			texture_set_x(&plan, x);
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

/** A test of the texture transform.
 * The texture transform will be executed on various sparse vectors.
 * The output is compared with the results produced by ::spherical_harmonic.
 * In case of error the execution will be interrupted and an error message 
 * will be printed containing
 * - the index of the test,
 * - the error (see ::l_2_dist) and
 * - the maximum absolute error in this case.
 * 
 * The parameter file contains (in this order):
 * - the bandwidth N,
 * - the number of non-zero components in the transformed vectors,
 * - the number of different longitudes and latitudes of the pole figures,
 * - the number of nodes,
 * - the number of different longitudes and latitudes of the nodes,
 * - the number of tests ands
 * - the maximum relative error. (see ::l_2_rel_dist)
 * 
 * @par inp - the filename of the parameter file.
 */
void linearity_test(const char *inp)
{
	FILE *inp_file = fopen(inp, "r");
	int N, h_phi_count, h_theta_count, N1, N2, r_phi_count, r_theta_count, w,
		it;
	double delta;
	char err_prefix[100];
	double _Complex *omega, *x;
	double _Complex *x_ref;
	double *h_phi, *h_theta, *r;
	texture_plan plan;
	int i, j, count, k;
	// unsigned short int seed[] = { 1, 2, 3 };

	sprintf(err_prefix, "linearity_test failed (%s):\n", inp);
	printf("*** linearity_test (%s)\n", inp);

	fscanf(inp_file, "%d%d%d%d%d%d%d%d%lg", &N, &w, &h_phi_count,
				 &h_theta_count, &N2, &r_phi_count, &r_theta_count, &it, &delta);
	N1 = h_phi_count * h_theta_count;

	texture_precompute(N);
	// seed48(seed);
	srand(0);

	omega =
		(double _Complex *) smart_malloc(texture_flat_length(N) *
																		sizeof(double _Complex));
	x = (double _Complex *) smart_malloc(N1 * N2 * sizeof(double _Complex));
	x_ref = (double _Complex *) smart_malloc(N1 * N2 * sizeof(double _Complex));

	h_phi = (double *) smart_malloc(N1 * sizeof(double));
	h_theta = (double *) smart_malloc(N1 * sizeof(double));
	r = (double *) smart_malloc(N1 * N2 * 2 * sizeof(double));
	initialize_angles(h_phi, h_theta, r, h_phi_count, h_theta_count, N2,
										r_phi_count, r_theta_count);

	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	for (count = 0; count < it; count++) {
		memset(omega, 0, texture_flat_length(N) * sizeof(double _Complex));
		memset(x_ref, 0, N1 * N2 * sizeof(double _Complex));

		for (k = 0; k < w; k++) {
			int l, m, n;
			double _Complex offset;

			l = rand() % (N + 1);
			m = (rand() % (2 * l + 1)) - l;
			n = (rand() % (2 * l + 1)) - l;
			offset = crand();

			omega[texture_flat_index(l, m, n)] += offset;

			for (i = 0; i < N1; i++) {
				for (j = 0; j < N2; j++) {
					x_ref[i * N2 + j] += offset
						* conj(spherical_harmonic(l, n, texture_get_h_phi(&plan)[i],
																			texture_get_h_theta(&plan)[i]))
						* spherical_harmonic(l, m, texture_get_r(&plan)[2 * (i * N2 + j)],
																 texture_get_r(&plan)[2 * (i * N2 + j) + 1]
						);
				}
			}
		}
		texture_set_omega(&plan, omega);

		texture_trafo(&plan);

		if (!equal_two_norm_rel(texture_get_x(&plan), x_ref,
														texture_get_x_length(&plan), delta)) {
			printf("%scount=%d diff=%lg delta=%lg\n", err_prefix, count,
						 l_2_dist(texture_get_x(&plan), x_ref,
											texture_get_x_length(&plan)),
						 delta * l_2_norm(x_ref, texture_get_x_length(&plan)));
			return;
		}
	}

	texture_finalize(&plan);

	free(omega);
	free(x);
	free(x_ref);

	free(h_phi);
	free(h_theta);
	free(r);

	fclose(inp_file);

	texture_forget();
}

/** A test of the adjoint texture transform.
 * The adjoint texture transform will be executed on various sparse vectors.
 * The output is compared with the results produced by ::spherical_harmonic.
 * In case of error the execution will be interrupted and an error message 
 * will be printed containing
 * - the index of the test,
 * - the error (see ::l_2_dist) and
 * - the maximum absolute error in this case.
 * 
 * The parameter file contains (in this order):
 * - the bandwidth N,
 * - the number of non-zero components in the transformed vectors,
 * - the number of different longitudes and latitudes of the pole figures,
 * - the number of nodes,
 * - the number of different longitudes and latitudes of the nodes,
 * - the number of tests ands
 * - the maximum relative error. (see ::l_2_rel_dist)
 * 
 * @par inp - the filename of the parameter file.
 */
void linearity_adjoint_test(const char *inp)
{
	FILE *inp_file = fopen(inp, "r");
	int N, h_phi_count, h_theta_count, N1, N2, r_phi_count, r_theta_count, w,
		it;
	double delta;
	char err_prefix[100];
	double _Complex *omega, *x;
	double _Complex *omega_ref;
	double *h_phi, *h_theta, *r;
	texture_plan plan;
	int count, k;
	// unsigned short int seed[] = { 1, 2, 3 };

	sprintf(err_prefix, "linearity_adjoint_test failed (%s):\n", inp);
	printf("*** linearity_adjoint_test (%s)\n", inp);

	fscanf(inp_file, "%d%d%d%d%d%d%d%d%lg", &N, &w, &h_phi_count,
				 &h_theta_count, &N2, &r_phi_count, &r_theta_count, &it, &delta);
	N1 = h_phi_count * h_theta_count;

	texture_precompute(N);
	// seed48(seed);
	srand(0);

	omega =
		(double _Complex *) smart_malloc(texture_flat_length(N) *
																		sizeof(double _Complex));
	x = (double _Complex *) smart_malloc(N1 * N2 * sizeof(double _Complex));
	omega_ref =
		(double _Complex *) smart_malloc(texture_flat_length(N) *
																		sizeof(double _Complex));

	h_phi = (double *) smart_malloc(N1 * sizeof(double));
	h_theta = (double *) smart_malloc(N1 * sizeof(double));
	r = (double *) smart_malloc(N1 * N2 * 2 * sizeof(double));
	initialize_angles(h_phi, h_theta, r, h_phi_count, h_theta_count, N2,
										r_phi_count, r_theta_count);

	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	for (count = 0; count < it; count++) {
		memset(omega_ref, 0, texture_flat_length(N) * sizeof(double _Complex));
		memset(x, 0, N1 * N2 * sizeof(double _Complex));

		for (k = 0; k < w; k++) {
			int i, j;
			int l, m, n;
			double _Complex offset;

			i = rand() % N1;
			j = rand() % N2;
			offset = crand();

			x[i * N2 + j] += offset;

			for (l = 0; l <= N; l++) {
				for (m = -l; m <= l; m++) {
					for (n = -l; n <= l; n++) {
						omega_ref[texture_flat_index(l, m, n)] += offset
							* spherical_harmonic(l, n, texture_get_h_phi(&plan)[i],
																	 texture_get_h_theta(&plan)[i])
							*
							conj(spherical_harmonic
									 (l, m, texture_get_r(&plan)[2 * (i * N2 + j)],
										texture_get_r(&plan)[2 * (i * N2 + j) + 1]));
					}
				}
			}
		}

		texture_set_x(&plan, x);

		texture_adjoint(&plan);

		if (!equal_two_norm_rel(texture_get_omega(&plan), omega_ref,
														texture_get_omega_length(&plan), delta)) {
			printf("%scount=%d diff=%lg delta=%lg\n", err_prefix, count,
						 l_2_dist(texture_get_omega(&plan), omega_ref,
											texture_get_omega_length(&plan)),
						 delta * l_2_norm(omega_ref, texture_get_omega_length(&plan)));
			return;
		}
	}

	texture_finalize(&plan);

	free(omega);
	free(x);
	free(omega_ref);

	free(h_phi);
	free(h_theta);
	free(r);

	fclose(inp_file);

	texture_forget();
}

/** A test for ::texture_precompute.
 * The function will be called with small values.
 */
void precompute_extreme_values_test()
{
	printf("*** precompute_extreme_values_test\n");

	texture_precompute(2);
	texture_forget();

	texture_precompute(1);
	texture_forget();

	texture_precompute(0);
	texture_forget();
}

/** A test for ::texture_trafo.
 * The function will be called with plans having small values for N, N1 or N2.
 */
void texture_trafo_extreme_values_test(const char *inp)
{
	FILE *inp_file = fopen(inp, "r");
	texture_plan plan;
	double _Complex *omega;
	double _Complex *x;
	double *h_phi, *h_theta, *r;
	int N, N1, N2;

	printf("*** texture_trafo_extreme_values_test (%s)\n", inp);

	fscanf(inp_file, "%d%d%d", &N, &N1, &N2);

	omega = smart_calloc(texture_flat_length(N), sizeof(double _Complex));
	x = smart_malloc(N1 * N2 * sizeof(double _Complex));
	h_phi = smart_calloc(N1, sizeof(double));
	h_theta = smart_calloc(N1, sizeof(double));
	r = smart_calloc(N1 * N2 * 2, sizeof(double));

	texture_precompute(N);

	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	texture_trafo(&plan);
	texture_adjoint(&plan);

	texture_finalize(&plan);

	texture_forget();
	free(omega);
	free(x);
	free(h_phi);
	free(h_theta);
	free(r);
	fclose(inp_file);
}

/** Print a message explaining how to run the tests.
 */
void usage()
{
	printf("A test for the texture transform.\n");
	printf("textureTest [Options]\n");
	printf("Options:\n");
	printf("-val : use small input values for usage with valgrind\n");
}

/** The testsuite.
 * See usage or run the programm with "-help" for more information about 
 * command line arguments.
 */
int main(int arglen, char *argv[])
{

	if (arglen == 2 && !strcmp(argv[1], "-val")) {
		// usage with valgrind

		printf("*** Test for usage with valgrind.\n");
		printf("*** If some output not preceded by *** is produced, ");
		printf("there is some error.\n");

		spherical_harmonic_test("spherical_harmonic_test.inp");

		nfsft_test("nfsft_moderate_test.inp");

		nfsft_test("nfsft_small_test.inp");

		unit_vector_test("unit_vector_moderate_test.inp");

		unit_vector_adjoint_test("unit_vector_adjoint_moderate_test.inp");

		linearity_test("linearity_moderate_test.inp");

		linearity_adjoint_test("linearity_adjoint_moderate_test.inp");

		precompute_extreme_values_test();

		texture_trafo_extreme_values_test
			("texture_trafo_extreme_values_test_1.inp");

		texture_trafo_extreme_values_test
			("texture_trafo_extreme_values_test_2.inp");

		texture_trafo_extreme_values_test
			("texture_trafo_extreme_values_test_3.inp");

		texture_trafo_extreme_values_test
			("texture_trafo_extreme_values_test_4.inp");

		simple_solver_test("simple_solver_moderate_test.inp");
	} else if (arglen == 1) {
		// default usage

		printf("*** Test with default parameters.\n");
		printf("*** If some output not preceded by *** is produced, ");
		printf("there is some error.\n");

		spherical_harmonic_test("spherical_harmonic_test.inp");

		nfsft_test("nfsft_test.inp");

		unit_vector_test("unit_vector_test.inp");

		unit_vector_adjoint_test("unit_vector_adjoint_test.inp");

		linearity_test("linearity_test.inp");

		linearity_adjoint_test("linearity_adjoint_test.inp");

		simple_solver_test("simple_solver_test.inp");
	} else {
		usage();
	}

	return 0;
}

/**
 * @}
 */
