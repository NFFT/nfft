/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: simple_test.c 3509 2010-05-25 19:00:59Z keiner $ */

/* Standard headers. */
#include <complex.h>
#include "nfft3.h"
#include "infft.h"

/* Testcase delegate. */
typedef struct testcase_delegate_s testcase_delegate_t;

typedef void (*setup_t)(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, C **f_hat, C **f);
typedef void (*destroy_t)(const testcase_delegate_t *ego_, R *x, C *f_hat, C *f);

struct testcase_delegate_s
{
  setup_t setup;
  destroy_t destroy;
};

typedef struct testcase_delegate_file_s
{
  setup_t setup;
  destroy_t destroy;
  const char *filename;
} testcase_delegate_file_t;

void X(setup_file)(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, C **f_hat, C **f);
void X(destroy_file)(const testcase_delegate_t *ego_, R *x, C *f_hat, C *f);

typedef struct testcase_delegate_online_s
{
  setup_t setup;
  destroy_t destroy;
  const int d;
  const int N;
  const int M;
} testcase_delegate_online_t;

void X(setup_online)(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, C **f_hat, C **f);
void X(destroy_online)(const testcase_delegate_t *ego_, R *x, C *f_hat, C *f);

/* Init delegate. */
typedef struct init_delegate_s init_delegate_t;
typedef void (*init_t)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);

struct init_delegate_s
{
  const char *name;
  init_t init;
  const int m;
  const unsigned nfft_flags;
  const unsigned fftw_flags;
};

/* Trafo delegate. */
typedef void (*trafo_t)(X(plan) *p);
typedef double (*cost_t)(X(plan) *p);
typedef const char* (*check_t)(X(plan) *p);
typedef R (*acc_t)(X(plan) *p);

typedef struct trafo_delegate_s
{
  const char *name;
  trafo_t trafo;
  check_t check;
  cost_t cost;
  acc_t acc;

} trafo_delegate_t;

double X(trafo_direct_cost)(X(plan) *p);

R X(err_trafo)(X(plan) *p);
R X(err_trafo_direct)(X(plan) *p);

/* Check single test case.*/
int X(check_single)(const testcase_delegate_t *testcase,
    init_delegate_t *init_delegate, trafo_delegate_t *trafo_delegate);

/* Check multiple test cases.*/
void X(check_many)(const size_t nf, const size_t ni, const size_t nt,
  const testcase_delegate_t **testcases, init_delegate_t **initializers,
  trafo_delegate_t **trafos);

/* Size of array. */
#define SIZE(x) sizeof(x)/sizeof(x[0])

/* Initializers. */
void X(init_1d_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
void X(init_2d_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
void X(init_3d_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
void X(init_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
void X(init_advanced_pre_psi_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);

#define DEFAULT_NFFT_FLAGS MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT | FFT_OUT_OF_PLACE
#define DEFAULT_FFTW_FLAGS FFTW_ESTIMATE | FFTW_DESTROY_INPUT

init_delegate_t init_1d;
init_delegate_t init_2d;
init_delegate_t init_3d;
init_delegate_t init;
init_delegate_t init_advanced_pre_psi;
init_delegate_t init_advanced_pre_full_psi;
init_delegate_t init_advanced_pre_lin_psi;
#if defined(GAUSSIAN)
init_delegate_t init_advanced_pre_fg_psi;
#endif

trafo_delegate_t trafo_direct;
trafo_delegate_t trafo;
