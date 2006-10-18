#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>

#include<util.h>
#include<nfft3.h>
#include"texture_util.h"

int main(int argc, char *argv[])
{
	FILE *f1 = fopen(argv[1], "r");
	FILE *f2 = fopen(argv[2], "r");
	char line1[MAX_LINE];
	char line2[MAX_LINE];

	double *h_phi, *h_theta, *r;
	double *h_phi_2, *h_theta_2, *r_2;
	complex *x, *omega;
	complex *x_2, *omega_2;
	int N, N1, N2;
	int N_2, N1_2, N2_2;

	fgets(line1, MAX_LINE, f1);
	fgets(line2, MAX_LINE, f2);
	fclose(f1);
	fclose(f2);
	f1 = fopen(argv[1], "r");
	f2 = fopen(argv[2], "r");
	if (!strcmp(line1, "Polefigures\n") && !strcmp(line2, "Polefigures\n")) {
		read_h(&N1, &h_phi, &h_theta, f1, stderr);
		read_h(&N1_2, &h_phi_2, &h_theta_2, f2, stderr);
		if (N1 == N1_2) {
			printf("%lg\n", nfft_error_l_2_double(h_phi, h_phi_2, N1));
			printf("%lg\n", nfft_error_l_2_double(h_theta, h_theta, N1));
		} else {
			printf("different N1\n");
		}
	} else if (!strcmp(line1, "Nodes\n") && !strcmp(line2, "Nodes\n")) {
		read_r(1, &N2, &r, f1, stderr);
		read_r(1, &N2_2, &r_2, f2, stderr);
		if (N2 == N2_2) {
			printf("%lg\n", nfft_error_l_2_double(r, r_2, N2 * 2));
		} else {
			printf("different N2\n");
		}
	} else if (!strcmp(line1, "Omega\n") && !strcmp(line2, "Omega\n")) {
		read_omega(&N, &omega, f1, stderr);
		read_omega(&N_2, &omega_2, f2, stderr);
		if (N == N_2) {
			printf("%lg\n",
						 nfft_error_l_2_complex(omega, omega_2, texture_flat_length(N)));
		} else {
			printf("different N\n");
		}
	} else if (!strcmp(line1, "Samples\n") && !strcmp(line2, "Samples\n")) {
		read_x(&N1, &N2, &x, f1, stderr);
		read_x(&N1_2, &N2_2, &x_2, f2, stderr);
		if (N1 == N1_2 && N2 == N2_2) {
			printf("%lg\n", nfft_error_l_2_complex(x, x_2, N1 * N2));
		} else {
			printf("different N1 / N2\n");
		}
	} else {
	}

	fclose(f1);
	fclose(f2);

	return 0;
}
