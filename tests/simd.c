/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* The SIMD kernels answer to the plain C ones.
 *
 * Every instruction set this build carries *and* this host supports runs the
 * same transform on the same inputs, and the result has to match the scalar
 * one. The two paths sum the same terms in a different order, so they agree to
 * rounding rather than exactly -- but nothing beyond rounding, which is what
 * separates a working vector kernel from one whose lane handling is wrong.
 *
 * On a host that offers nothing but the scalar path (or in a build configured
 * with --disable-simd, or in long-double precision) the comparison loop has a
 * single pass and the test degenerates to a self-check. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "simd.h" /* tests/simd.h; include/simd.h comes via infft.h */

/* Relative agreement demanded between two instruction sets. Both sums run over
 * the same (2m+2)^d window terms, so the gap is a handful of ulps; the margin
 * is wide enough not to depend on the association a compiler picks. */
#define SIMD_TOL (K(1024.0) * EPSILON)

/* Deterministic inputs: the runs being compared must differ in the instruction
 * set and in nothing else, so the data cannot come from rand(). */
static unsigned int lcg_state;

static void lcg_seed(unsigned int s)
{
  lcg_state = s;
}

/** Next pseudo random value in [-1/2, 1/2). */
static R lcg_next(void)
{
  lcg_state = lcg_state * 1103515245u + 12345u;
  return (R)((lcg_state >> 8) & 0xFFFFu) / K(65536.0) - K(0.5);
}

static void fill_nodes(R *x, const INT len)
{
  INT j;

  for (j = 0; j < len; j++)
    x[j] = lcg_next();
}

static void fill_complex(C *v, const INT len)
{
  INT j;

  for (j = 0; j < len; j++)
    v[j] = lcg_next() + II * lcg_next();
}

/** Runs one forward and one adjoint transform for every instruction set
 *  available here and checks each against the scalar result. */
static void compare_isa(const int d, const int *N, const int M,
    const unsigned flags)
{
  const int isa_saved = Y(simd_isa)();
  int n[3], NN = 1, t, isa, runs = 0;
  C *f_ref = NULL, *f_hat_ref = NULL;

  for (t = 0; t < d; t++)
  {
    INT n2;
    INT e;

    Y(next_power_of_2_exp)((INT)(N[t]), &n2, &e);
    n[t] = (int)(2 * n2);
    NN *= N[t];
  }

  f_ref = (C*)Y(malloc)((size_t)(M) * sizeof(C));
  f_hat_ref = (C*)Y(malloc)((size_t)(NN) * sizeof(C));

  for (isa = NFFT_SIMD_SCALAR; isa <= NFFT_SIMD_MAX; isa++)
  {
    X(plan) p;

    if (!Y(simd_isa_available)(isa))
      continue;

    /* Kernels are bound when the plan is initialised, so this has to come
     * first. */
    Y(simd_force_isa)(isa);
    CU_ASSERT(Y(simd_isa)() == isa);

    X(init_guru)(&p, d, (int*)N, M, n, 6, flags, FFTW_ESTIMATE
        | FFTW_DESTROY_INPUT);

    lcg_seed(4711u);
    fill_nodes(p.x, (INT)(M) * (INT)(d));
    fill_complex(p.f_hat, (INT)NN);

    X(precompute_one_psi)(&p);

    X(trafo)(&p);

    if (runs == 0)
      memcpy(f_ref, p.f, (size_t)(M) * sizeof(C));
    else
    {
      const R err = Y(error_l_infty_complex)(f_ref, p.f, (INT)M);
      if (!(err <= SIMD_TOL))
        fprintf(stderr, "simd: trafo d=%d %s vs scalar: err = " __FE__ "\n", d,
            Y(simd_isa_name)(isa), err);
      CU_ASSERT(err <= SIMD_TOL);
    }

    lcg_seed(112211u);
    fill_complex(p.f, (INT)M);

    X(adjoint)(&p);

    if (runs == 0)
      memcpy(f_hat_ref, p.f_hat, (size_t)(NN) * sizeof(C));
    else
    {
      const R err = Y(error_l_infty_complex)(f_hat_ref, p.f_hat, (INT)NN);
      if (!(err <= SIMD_TOL))
        fprintf(stderr, "simd: adjoint d=%d %s vs scalar: err = " __FE__ "\n",
            d, Y(simd_isa_name)(isa), err);
      CU_ASSERT(err <= SIMD_TOL);
    }

    X(finalize)(&p);
    runs++;
  }

  CU_ASSERT(runs >= 1);

  Y(free)(f_hat_ref);
  Y(free)(f_ref);

  Y(simd_force_isa)(isa_saved);
}

/* The flag sets worth separating: PRE_PSI hands the kernels a precomputed
 * window run, without it the run is built per node, and the blockwise adjoint
 * is a different kernel again. */
#define SIMD_FLAGS_BASE (PRE_PHI_HUT | MALLOC_X | MALLOC_F_HAT | MALLOC_F \
    | FFTW_INIT | FFT_OUT_OF_PLACE)

static void check_simd_d(const int d, const int *N, const int M)
{
  compare_isa(d, N, M, SIMD_FLAGS_BASE | PRE_PSI);
  compare_isa(d, N, M, SIMD_FLAGS_BASE);
#ifdef _OPENMP
  compare_isa(d, N, M, SIMD_FLAGS_BASE | PRE_PSI
      | NFFT_OMP_BLOCKWISE_ADJOINT | NFFT_SORT_NODES);
#endif
}

void X(check_simd_detection)(void)
{
  const int isa = Y(simd_isa)();

  CU_ASSERT(isa >= NFFT_SIMD_SCALAR);
  CU_ASSERT(isa <= NFFT_SIMD_MAX);
  CU_ASSERT(Y(simd_isa_available)(isa));
  CU_ASSERT(Y(simd_isa_available)(NFFT_SIMD_SCALAR));
  CU_ASSERT(strcmp(Y(simd_isa_name)(isa), "unknown") != 0);
  CU_ASSERT(strcmp(Y(simd_isa_name)(NFFT_SIMD_SCALAR), "scalar") == 0);

  /* An instruction set this build does not carry can never be selected. */
  Y(simd_force_isa)(NFFT_SIMD_MAX + 1);
  CU_ASSERT(Y(simd_isa)() == NFFT_SIMD_SCALAR);
  Y(simd_force_isa)(isa);
  CU_ASSERT(Y(simd_isa)() == isa);
}

void X(check_simd_1d)(void)
{
  const int N[1] = { 64 };
  check_simd_d(1, N, 97);
}

void X(check_simd_2d)(void)
{
  const int N[2] = { 24, 20 };
  check_simd_d(2, N, 61);
}

void X(check_simd_3d)(void)
{
  const int N[3] = { 12, 10, 14 };
  check_simd_d(3, N, 43);
}
