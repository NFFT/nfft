/*
 * Copyright (c) 2002, 2012 Jens Keiner, Daniel Potts, Stefan Kunis
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* $Id$ */

/**
 * \file nfsft.h
 *
 * \brief Header file for functions related to associated Legendre
 *        functions/polynomials
 *
 * \author Jens Keiner
 */
#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "api.h"

/** \addtogroup nfsft
 * \{
 */

void alpha_al_row(double *alpha, int N, int n);
void beta_al_row(double *beta, int N, int n);
void gamma_al_row(double *gamma, int N, int n);

/**
 * Compute three-term-recurrence coefficients \f$\alpha_{k-1}^n\f$ of associated
 * Legendre functions for \f$k,n = 0,1,\ldots,N\f$.
 *
 * \arg alpha A pointer to an array of doubles of size \f$(N+1)^2\f$ where the
 *   coefficients will be stored such that alpha[n+(N+1)+k] =
 *   \f$\alpha_{k-1}^n\f$.
 * \arg N The upper bound \f$N\f$.
 */
void alpha_al_all(double *alpha, int N);

/**
 * Compute three-term-recurrence coefficients \f$\beta_{k-1}^n\f$ of associated
 * Legendre functions for \f$k,n = 0,1,\ldots,N\f$.
 *
 * \arg beta A pointer to an array of doubles of size \f$(N+1)^2\f$ where the
 *   coefficients will be stored such that beta[n+(N+1)+k] =
 *   \f$\beta_{k-1}^n\f$.
 * \arg N The upper bound \f$N\f$.
 */
void beta_al_all(double *beta, int N);

/**
 * Compute three-term-recurrence coefficients \f$\gamma_{k-1}^n\f$ of associated
 * Legendre functions for \f$k,n = 0,1,\ldots,N\f$.
 *
 * \arg beta A pointer to an array of doubles of size \f$(N+1)^2\f$ where the
 *   coefficients will be stored such that gamma[n+(N+1)+k] =
 *   \f$\gamma_{k-1}^n\f$.
 * \arg N The upper bound \f$N\f$.
 */
void gamma_al_all(double *gamma, int N);

/**
 * Evaluates an associated Legendre polynomials \f$P_k^n(x,c)\f$ using the
 * Clenshaw-algorithm.
 *
 * \arg x A pointer to an array of nodes where the function is to be evaluated
 * \arg y A pointer to an array where the function values are returned
 * \arg size The length of x and y
 * \arg k The index \f$k\f$
 * \arg alpha A pointer to an array containing the recurrence coefficients
 *   \f$\alpha_c^n,\ldots,\alpha_{c+k}^n\f$
 * \arg beta A pointer to an array containing the recurrence coefficients
 *   \f$\beta_c^n,\ldots,\beta_{c+k}^n\f$
 * \arg gamma A pointer to an array containing the recurrence coefficients
 *   \f$\gamma_c^n,\ldots,\gamma_{c+k}^n\f$
 */
void eval_al(double *x, double *y, int size, int k, double *alpha,
  double *beta, double *gamma);

/**
 * Evaluates an associated Legendre polynomials \f$P_k^n(x,c)\f$ using the
 * Clenshaw-algorithm if it no exceeds a given threshold.
 *
 * \arg x A pointer to an array of nodes where the function is to be evaluated
 * \arg y A pointer to an array where the function values are returned
 * \arg size The length of x and y
 * \arg k The index \f$k\f$
 * \arg alpha A pointer to an array containing the recurrence coefficients
 *   \f$\alpha_c^n,\ldots,\alpha_{c+k}^n\f$
 * \arg beta A pointer to an array containing the recurrence coefficients
 *   \f$\beta_c^n,\ldots,\beta_{c+k}^n\f$
 * \arg gamma A pointer to an array containing the recurrence coefficients
 *   \f$\gamma_c^n,\ldots,\gamma_{c+k}^n\f$
 * \arg threshold The threshold
 */
int eval_al_thresh(double *x, double *y, int size, int k, double *alpha,
  double *beta, double *gamma, double threshold);
/* \} */
#endif
