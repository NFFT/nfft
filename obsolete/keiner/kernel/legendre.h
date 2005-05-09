/**
 * Header file for functions related to associated Legendre functions/polynomials
 */
#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "api.h"

/** \defgroup nfsft_legendre NFSFT: Functions related to associated Legendre functions/polynomials */

/**
 *
 * 
 * \ingroup nfsft_legendre
 */
void alpha_al_all(double *alpha, int N);

/**
 *
 * 
 * \ingroup nfsft_legendre
 */
void beta_al_all(double *beta, int N);

/**
 *
 * 
 * \ingroup nfsft_legendre
 */
void gamma_al_all(double *gamma, int N);

/**
 *
 * 
 * \ingroup nfsft_legendre
 */
void gamma_al_m1_all(double *gamma, int N);

/**
 * Computes three-term recurrence coefficients alpha for associated Legendre 
 * functions. 
 * 
 * \ingroup nfsft_legendre
 */ 
double alpha_al (int k, int n);

/**
 * Computes three-term recurrence coefficients beta for associated Legendre 
 * functions. 
 * 
 * \ingroup nfsft_legendre
 */ 
double beta_al (int k, int n);

/**
 * Computes three-term recurrence coefficients gamma for associated Legendre 
 * functions. 
 * 
 * \ingroup nfsft_legendre
 */ 
double gamma_al (int k, int n);

/** 
 * Initializes three-term recurrence coefficients of associated Legendre 
 * functions. 
 * 
 * \ingroup nfsft_legendre
 */
void init_al_coeff(double *alpha, double *beta, double *gamma,
  int l, int m, int n);

/**
 * Evaluates associated Legendre functions using the Clenshaw-algorithm. 
 * 
 * \ingroup nfsft_legendre
 */
void eval_al(double *x, double *y, int size, int m, double *alpha, 
  double *beta, double *gamma);

/**
 * Evaluates associated Legendre functions using the Clenshaw-algorithm. 
 * 
 * \ingroup nfsft_legendre
 */
bool eval_al_thresh(double *x, double *y, int size, int m, double *alpha, 
  double *beta, double *gamma, double threshold);
#endif
