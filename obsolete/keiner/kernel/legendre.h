#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "api.h"

void alpha_al_all(double *alpha, int N);
void beta_al_all(double *alpha, int N);
void gamma_al_all(double *alpha, int N);
void gamma_al_m1_all(double *alpha, int N);

/**
 * Computes three-term recurrence coefficients alpha for associated Legendre 
 * functions. 
 */ 
double alpha_al (int k, int n);

/**
 * Computes three-term recurrence coefficients beta for associated Legendre 
 * functions. 
 */ 
double beta_al (int k, int n);

/**
 * Computes three-term recurrence coefficients gamma for associated Legendre 
 * functions. 
 */ 
double gamma_al (int k, int n);

/** 
 * Initializes three-term recurrence coefficients of associated Legendre 
 * functions. 
 */
void init_al_coeff(double *alpha, double *beta, double *gamma,
  int l, int m, int n);

/**
 * Evaluates associated Legendre functions using the Clenshaw-algorithm. 
 */
void eval_al(double *x, double *y, int size, int m, double *alpha, 
  double *beta, double *gamma);

/**
 * Evaluates associated Legendre functions using the Clenshaw-algorithm. 
 */
bool eval_al_thresh(double *x, double *y, int size, int m, double *alpha, 
  double *beta, double *gamma, double threshold);
#endif
