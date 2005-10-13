#include<complex.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include<nfft3.h>

int N, N1, N2;
int type;
int iterations;

complex *omega, *x;
double *h_phi, *h_theta, *r;
texture_plan plan;
itexture_plan iplan;

void usage()
{
	printf("'performance N N1 N2 type iterations' carries out some ");
	printf("texture_trafos.\n");
	printf("It has to be run with exactly 5 parameters as above.\n");
	printf("N N1 N2 are the parameters for the texture_plan");
	printf("type determines the type of transformation that will be run.\n");
	printf("(0 = trafo, 1 = adjoint, 2 = inverse)\n");
	printf("iterations gives the number of transformations.\n");
}

void parse_input(int arglen, char *argv[])
{
	if (arglen == 6) {
		N = atoi(argv[1]);
		N1 = atoi(argv[2]);
		N2 = atoi(argv[3]);
		type = atoi(argv[4]);
		iterations = atoi(argv[5]);
	} else {
		usage();
		exit(-1);
	}
}

void init()
{
	int i;
	unsigned short int seed[] = { 1, 2, 3 };

	seed48(seed);

	omega = (complex *) malloc(texture_flat_length(N) * sizeof(complex));
	x = (complex *) malloc(N1 * N2 * sizeof(complex));
	h_phi = (double *) malloc(N1 * sizeof(double));
	h_theta = (double *) malloc(N1 * sizeof(double));
	r = (double *) malloc(N1 * N2 * 2 * sizeof(double));

	for (i = 0; i < texture_flat_length(N); i++) {
		omega[i] = rand() * drand48() + I * rand() * drand48();
	}
	for (i = 0; i < N1 * N2; i++) {
		x[i] = rand() * drand48() + I * rand() * drand48();
	}
	for (i = 0; i < N1; i++) {
		h_phi[i] = (drand48() - 0.5) * TEXTURE_MAX_ANGLE;
		h_theta[i] = (drand48() / 2) * TEXTURE_MAX_ANGLE;
	}
	for (i = 0; i < N1 * N2; i++) {
		r[2 * i] = (drand48() - 0.5) * TEXTURE_MAX_ANGLE;
		r[2 * i + 1] = (drand48() / 2) * TEXTURE_MAX_ANGLE;
	}

	texture_precompute(N);
	
	texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

	itexture_init(&iplan, &plan);
	memcpy(iplan.y, texture_get_x(&plan), N1 * N2 * sizeof(complex));
	memset(iplan.f_hat_iter, 0, texture_flat_length(N) * sizeof(complex));
	itexture_before_loop(&iplan);
}

void run_test()
{
	int count;

	switch (type) {
		case 0:
			for (count = 0; count < iterations; count++) {
				texture_trafo(&plan);
			}
			break;
		case 1:
			for (count = 0; count < iterations; count++) {
				texture_adjoint(&plan);
			}
			break;
		case 2:
			for (count = 0; count < iterations; count++) {
				itexture_loop_one_step(&iplan);
			}
			break;
	}
}

void cleanup()
{
	itexture_finalize(&iplan);
	texture_finalize(&plan);

	texture_forget();

	free(omega);
	free(x);
	free(h_phi);
	free(h_theta);
	free(r);
}

int main(int arglen, char *argv[])
{
	parse_input(arglen, argv);
	init();
	run_test();
	cleanup();

	return 0;
}
