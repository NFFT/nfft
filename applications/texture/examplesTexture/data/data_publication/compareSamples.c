#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include <texture_util.h>

int N1, N2, N1_ref, N2_ref;
complex *x, *x_ref;

int main(int argic, char *argv[])
{
	FILE *f1, *f2;
	double min_k, max_k;
	int i;

	f1 = fopen(argv[1], "r");
	f2 = fopen(argv[2], "r");
	read_x(&N1, &N2, &x, f1, stderr);
	read_x(&N1_ref, &N2_ref, &x_ref, f2, stderr);
	assert(N1 == N1_ref && N2 == N2_ref);

	fprintf(stderr, "relative l_2 distance: %.2e\n",
					l_2_rel_dist(x, x_ref, N1 * N2));

	min_k = max_k = cabs(x_ref[0] / x[0]);
	for (i = 0; i < N1 * N2; i++) {
		double quot = cabs(x_ref[i] / x[i]);
		if (quot < min_k) {
			min_k = quot;
		}
		if (quot > max_k) {
			max_k = quot;
		}
	}

	fprintf(stderr, "x_ref[i] / x[i] is between %.2e and %.2e\n", min_k, max_k);

	free(x);
	free(x_ref);
	fclose(f1);
	fclose(f2);
	return 0;
}
