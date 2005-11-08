#include<stdio.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include<stdlib.h>
#include<complex.h>

#define MAX_LINE 100

int N1, N2;
complex *x;

void read_samples(const char *file)
{
	int i, j;
	FILE *f = fopen(file, "r");
	char line[MAX_LINE];
	char blank;

	fgets(line, MAX_LINE, f);
	if (strcmp(line, "Samples\n")) {
		fprintf(stderr, "Invalid sample file!\n");
		fflush(0);
		exit(-1);
	}

	printf("# From the sample file:\n");

	fgets(line, MAX_LINE, f);
	while (line[0] == '#') {
		printf("%s", line);
		fgets(line, MAX_LINE, f);
	}

	if (line[0] != '\n') {
		fprintf(stderr, "Invalid header in sample file!\n");
		fflush(0);
		exit(-1);
	}

	fscanf(f, "%d%d", &N1, &N2);

	printf("# N1=%d N2=%d\n", N1, N2);
	printf("#\n");

	x = (complex *) malloc(N1 * N2 * sizeof(complex));
	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			double r, im;
			if (!fscanf(f, "%lg", &r)) {
				fprintf(stderr, "Parse error in sample file!\n");
				fflush(0);
				exit(-1);
			} else {
				fscanf(f, " + %lgi", &im);
				x[i * N2 + j] = r + I * im;
			}
		}
	}

	do {
		blank = fgetc(f);
	} while (isspace(blank));

	if (!feof(f)) {
		fprintf(stderr, "Parse error at the end of sample file!\n");
		fflush(0);
		exit(-1);
	}

	fclose(f);
}

void cleanup()
{
	free(x);
}

void test_abs()
{
	int i, j;
	double maxabs = -1;
	int i_max = -1, j_max = -1;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			double absval = cabs(x[i * N2 + j]);
			if (absval > maxabs) {
				maxabs = absval;
				i_max = i;
				j_max = j;
			}
		}
	}

	printf("maximum absolute value: %lg at i=%d j=%d\n", maxabs, i_max, j_max);
}

void test_real()
{
	int i, j;
	double maxim = -1;
	int i_max = -1, j_max = -1;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			double im = fabs(cimag(x[i * N2 + j]));
			if (im > maxim) {
				maxim = im;
				i_max = i;
				j_max = j;
			}
		}
	}

	printf("maximum imaginary part: %lg at i=%d j=%d\n", maxim, i_max, j_max);
}

void test_positiv()
{
	int i, j;
	double maxneg = -1;
	int i_max = -1, j_max = -1;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			double neg = (creal(x[i * N2 + j]) < 0) ? -creal(x[i * N2 + j]) : 0;
			if (neg > maxneg) {
				maxneg = neg;
				i_max = i;
				j_max = j;
			}
		}
	}

	printf("maximum negative part: %lg at i=%d j=%d\n", maxneg, i_max, j_max);
}

int main(int argc, char *argv[])
{
	const char *sample_file = "samples.out";

	if (argc > 1) {
		sample_file = argv[1];
	}

	if (argc <= 2) {
		read_samples(sample_file);
		test_abs();
		test_real();
		test_positiv();
		cleanup();
	}
	return 0;
}
