#ifndef TRANSFORM_H
#define TRANSFORM_H

#include<complex.h>
#include<nfft3.h>

/**
 * The structure for the transform plan.
 */
struct texture_plan {

	MACRO_MV_PLAN(complex);

	int N;
	int N1;
	int N2;
	double* h;
	double** r;
	complex* omega;
	complex* x;
};

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
void texture_trafo(texture_plan *ths);

/**
 * Simple initialisation of a plan.
 */
void texture_init(texture_plan *ths);

/*
 * Initialisation of a plan.
 */
void texture_init_advanced(texture_plan *ths);

/**
 * The adjoint version of the transform.
 */
void texture_adjoint(texture_plan *ths);

#endif
