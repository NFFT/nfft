#ifndef TRANSFORM_H
#define TRANSFORM_H

#include<complex.h>
#include<nfft3.h>

/**
 * The structure for the transform plan.
 */
typedef struct texture_plan_ {

	MACRO_MV_PLAN(complex);

	int N;
	int N1;
	int N2;
	double *h_phi;
	double *h_theta;
	double **r_phi;
	double **r_theta;
} texture_plan;

/**
 * Convert a non-flat index to a flat index.
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
 * Use texture_finalize to free allocated memory.
 */
void texture_init(texture_plan *ths, int N, int N1, int N2);

/**
 * Advanced initialisation of a plan.
 */
void texture_init_advanced(texture_plan *ths, int N, int N1, int N2);

/**
 * Frees all memory allocated by texture_init or texture_init_advanced.
 */
void texture_finalize(texture_plan *ths);

/**
 * The adjoint version of the transform.
 */
void texture_adjoint(texture_plan *ths);

#endif
