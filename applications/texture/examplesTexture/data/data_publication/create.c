#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<string.h>

#include<texture_util.h>

#include"data_util.h"

int N1, N2;
double *h_phi, *h_theta, *r;
double _Complex *x;

int main()
{
	FILE *fid = fopen("grid_size", "r");
	int x1, x2, x3;
	int i, j;
	double factor;

	fscanf(fid, "%d%d", &N1, &N2);
	fclose(fid);

	h_phi = smart_malloc(N1 * sizeof(double));
	h_theta = smart_malloc(N1 * sizeof(double));
	r = smart_malloc(2 * N1 * N2 * sizeof(double));
	x = smart_malloc(N1 * N2 * sizeof(double _Complex));

	fid = fopen("h_files", "r");
	for (i = 0; i < N1; i++) {
		double norm;
		fscanf(fid, "%*[^0-9]%1d%1d%1d", &x1, &x2, &x3);
		norm = sqrt(x1 * x1 + x2 * x2 + x3 * x3);

		h_theta[i] = normalise_theta(acos(x3 / norm));
		h_phi[i] = normalise_phi(atan2(x2, x1));
	}
	fclose(fid);

	fid = fopen("factor", "r");
	fscanf(fid, "%lg", &factor);
	fclose(fid);

	fid = fopen("h_files", "r");
	for (i = 0; i < N1; i++) {
		FILE *pfid;
		char pfname[20];
		char pfpath[40];

		fscanf(fid, "%s\n", pfname);
		sprintf(pfpath, "data/%s", pfname);
		pfid = fopen(pfpath, "r");

		for (j = 0; j < N2; j++) {
			double dummy;
			fscanf(pfid, "%lg%lg%lg", &(r[2 * (i * N2 + j) + 1]),
						 &(r[2 * (i * N2 + j)]), &dummy);
			x[i * N2 + j] = dummy;
		}
		fclose(pfid);
		// memcpy(&(x[i * N2 + N2 / 2]), &(x[i * N2]), N2 / 2 * sizeof(double _Complex));
	}
	fclose(fid);
	// expand_r(N2, r);
	normalise_r(N2, r);

	fid = fopen("h_file", "w");
	fprintf(fid, "Polefigures\n");
	fprintf(fid, "# Data set for the publication.\n");
	write_h(N1, h_phi, h_theta, fid);
	fclose(fid);

	// The result of the following piece of code is that all nodesets are
	// equal.
	/* 
	   for (i = 0; i < N1; i++) { char name[20];

	   sprintf (name, "r_file%d", i);

	   fid = fopen (name, "w"); write_r (N2/2, &(r[2 * i * N2]), fid); fclose
	   (fid); } */

	fid = fopen("r_file", "w");
	fprintf(fid, "Nodes\n");
	fprintf(fid, "# Data set for the publication.\n");
	// fprintf(fid, "# The grid is mirrored.\n");
	write_r(N2, r, fid);
	fclose(fid);

	for (j = 0; j < N2; j++) {
		r[2 * j] *= -1;
	}

	fid = fopen("r_file_flipped", "w");
	fprintf(fid, "Nodes\n");
	fprintf(fid, "# Data set for the publication.\n");
	fprintf(fid, "# \\phis (\\rhos) are flipped.\n");
	// fprintf(fid, "# The grid is mirrored.\n");
	write_r(N2, r, fid);
	fclose(fid);

	fid = fopen("samples", "w");
	fprintf(fid, "Samples\n");
	fprintf(fid, "# Data set for the publication.\n");
	// fprintf(fid, "# Samples are mirrored.\n");
	write_x(N1, N2, x, fid);
	fclose(fid);

	for (i = 0; i < N1 * N2; i++) {
		x[i] *= factor;
	}

	fid = fopen("samples_unnormed", "w");
	fprintf(fid, "Samples\n");
	fprintf(fid, "# Data set for the publication.\n");
	fprintf(fid, "# Integral on the sphere should be %.2e.\n", factor);
	// fprintf(fid, "# Samples are mirrored.\n");
	write_x(N1, N2, x, fid);
	fclose(fid);

	free(h_phi);
	free(h_theta);
	free(r);
	free(x);

	return 0;
}
