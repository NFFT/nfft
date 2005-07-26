#ifndef TRANSFORM_H
#define TRANSFORM_H

#include<complex.h>
#include<nfft3.h>

/**
 * The structure for the transform plan.
 */
struct transform_plan_s {
	int N;
	int N1;
	int N2;
	double* h;
	double** r;
	complex* omega;
	complex* x;
};

/**
 * Datatype for the transform.
 */
typedef struct transform_plan_s *transform_plan;

/**
 * Konvert a non-flat index to a flat index.
 * \arg l the frequence
 * \arg m ranges from -l to l
 * \arg n ranges from -l to l
 */
int flatIndex(int l, int m, int n);

/**
 * Carries out the transform.
 */
void transform(transform_plan);

#endif
