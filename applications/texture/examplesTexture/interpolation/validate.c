#include<stdio.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include<stdlib.h>
#include<complex.h>

#include <texture_util.h>

int N1, N2;
complex *x;

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

void usage()
{
	// TODO
}

int main(int argc, char *argv[])
{
	const char *sample_file = "samples.out";
	FILE *f;

	if (argc > 1) {
		sample_file = argv[1];
	}

	if (argc <= 2) {
		// read_samples(sample_file);
		f = fopen(sample_file, "r");
		read_x(&N1, &N2, &x, f, stdout);
		fclose(f);

		test_abs();

		test_real();

		test_positiv();

		cleanup();
	} else {
		usage();
	}

	return 0;
}
