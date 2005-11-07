#ifndef DPT_H_
#define DPT_H_

#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>

#define DPT_NO_STABILIZATION (1U << 0)
#define DPT_BANDWIDTH_WINDOW (1U << 1)

typedef struct dpt_step_
{
  bool stable;                            /**< Indicates if the values 
                                               contained represent a fast or 
                                               a slow stabilized step.       */
  double **a11,**a12,**a21,**a22;         /**< The matrix components         */
} dpt_step;

typedef struct dpt_data_
{
  dpt_step **steps;                       /**< The cascade summation steps   */
} dpt_data;

typedef struct dpt_set_
{
  int flags;                              /**< The flags                     */
  int M;                                  /**< The number of DPT transforms  */
  int N;                                  /**< The transform length. Must be 
                                               a power of two.               */
  int t;                                  /**< The exponent of N             */
  dpt_data *dpt;                          /**< The DPT transform data        */
  double **xcvecs;                        /**< Array of pointers to arrays 
                                               containing the Chebyshev 
                                               nodes                         */
  double *xc;                             /**< Array for Chebychev-nodes.    */  
} dpt_set_s;

typedef dpt_set_s *dpt_set;

dpt_set dpt_init(const int M, const int t, const int flags);

void dpt_precompute(dpt_set set, const int m, double const* alpha, 
                    double const* beta, double const* gamma, int k_start,
                    double threshold);

void dpt_trafo(dpt_set set, const int m, complex *x);

void dpt_transposed(dpt_set set, const int m, complex *x);

void dpt_finalize(dpt_set set); 

#endif /*DPT_H_*/

