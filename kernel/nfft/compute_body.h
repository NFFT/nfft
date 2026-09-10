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

/*! \file compute_body.h
 *  \brief The per-node convolution kernels of the 1D, 2D and 3D NFFT, written
 *  once and compiled once per SIMD variant.
 *
 *  No include guard: nfft.c includes this after simd_run.h has defined the
 *  primitives for one instruction set, and does so once per variant it builds.
 *  Names come out as nfft_trafo_1d_compute_scalar, ..._avx2 and so on;
 *  nfft.c picks between them at run time.
 *
 *  Every kernel here walks the same shape. A node touches a 2m+2 point window
 *  run per dimension, wrapped around the periodic oversampled grid; along the
 *  last dimension the grid points it hits are contiguous, so that innermost
 *  loop is exactly the real-weight/complex-data dot product (forward) or
 *  scaled addition (adjoint) that simd_run.h vectorises. The outer dimensions
 *  only pick the row, and their window factors are hoisted out of the inner
 *  loop instead of being re-multiplied per grid point.
 */

/** One dimension's window run, expressed as the (at most two) contiguous
 *  pieces it occupies on a periodic grid of n points. The 2m+2 points wrap
 *  around the end of the grid at most once, since the transforms require
 *  2m+2 <= n. */
NFFT_SIMD_TGT static inline void NFFT_SIMD_ID(nfft_run_init)(nfft_run *r,
    const R x, const INT n, const INT m)
{
  INT u, o;

  uo2(&u, &o, x, n, m);

  r->start[0] = u;

  if (u < o)
  {
    r->len[0] = 2 * m + 2;
    r->pieces = 1;
  }
  else
  {
    /* o is u + 2m+1 - n here, so the first piece reaches the end of the grid
     * and the second one restarts at index 0. */
    r->len[0] = 2 * m + 1 - o;
    r->start[1] = 0;
    r->len[1] = o + 1;
    r->pieces = 2;
  }
}

/** sum_l psi[l] * g[run point l]. */
NFFT_SIMD_TGT static inline C NFFT_SIMD_ID(nfft_run_cdot)(const C *g,
    const nfft_run *r, const R *psi)
{
  C s = NFFT_SIMD_ID(nfft_simd_cdot)(psi, g + r->start[0], r->len[0]);

  if (r->pieces == 2)
    s += NFFT_SIMD_ID(nfft_simd_cdot)(psi + r->len[0], g + r->start[1],
        r->len[1]);

  return s;
}

/** g[run point l] += psi[l] * f. */
NFFT_SIMD_TGT static inline void NFFT_SIMD_ID(nfft_run_caxpy)(C *g,
    const nfft_run *r, const R *psi, const C f)
{
  NFFT_SIMD_ID(nfft_simd_caxpy)(g + r->start[0], psi, f, r->len[0]);

  if (r->pieces == 2)
    NFFT_SIMD_ID(nfft_simd_caxpy)(g + r->start[1], psi + r->len[0], f,
        r->len[1]);
}

/* ------------------------------------------------------------------ 1D --- */

NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_trafo_1d_compute)(C *fj,
    const C *g, const R *psij_const, const R *xj, const INT n, const INT m)
{
  nfft_run r;

  NFFT_SIMD_ID(nfft_run_init)(&r, *xj, n, m);

  *fj = NFFT_SIMD_ID(nfft_run_cdot)(g, &r, psij_const);
}

#ifndef _OPENMP
NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_adjoint_1d_compute_serial)(
    const C *fj, C *g, const R *psij_const, const R *xj, const INT n,
    const INT m)
{
  nfft_run r;

  NFFT_SIMD_ID(nfft_run_init)(&r, *xj, n, m);

  NFFT_SIMD_ID(nfft_run_caxpy)(g, &r, psij_const, *fj);
}
#endif

#ifdef _OPENMP
/** The part of one node's run in the first dimension that falls into the
 *  index range [u, o] the calling thread owns, empty when o < u. */
NFFT_SIMD_TGT static inline void NFFT_SIMD_ID(nfft_block_caxpy)(C *g,
    const R *psi, const C f, const INT u, const INT o, const INT offset_psij)
{
  if (o >= u)
    NFFT_SIMD_ID(nfft_simd_caxpy)(g + u, psi + offset_psij, f, o - u + 1);
}

/**
 * Adjoint NFFT for one-dimensional case updating only a specified range of
 * vector g.
 *
 * \arg f input coefficient f[j]
 * \arg g output vector g
 * \arg psij_const vector of window function values
 * \arg xj node x[j]
 * \arg n FFTW length (number oversampled Fourier coefficients)
 * \arg m window length
 * \arg my_u0 lowest index the current thread writes to in g
 * \arg my_o0 highest index the current thread writes to in g
 *
 * \author Toni Volkmer
 */
NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_adjoint_1d_compute_omp_blockwise)(
    const C f, C *g, const R *psij_const, const R *xj, const INT n, const INT m,
    const INT my_u0, const INT my_o0)
{
  INT ar_u, ar_o;

  uo2(&ar_u, &ar_o, *xj, n, m);

  if (ar_u < ar_o)
  {
    const INT u = MAX(my_u0, ar_u);
    const INT o = MIN(my_o0, ar_o);

    NFFT_SIMD_ID(nfft_block_caxpy)(g, psij_const, f, u, o, u - ar_u);
  }
  else
  {
    const INT u = MAX(my_u0, ar_u);

    NFFT_SIMD_ID(nfft_block_caxpy)(g, psij_const, f, u, my_o0, u - ar_u);
    NFFT_SIMD_ID(nfft_block_caxpy)(g, psij_const, f, my_u0,
        MIN(my_o0, ar_o), u - ar_u + my_u0 - ar_u + n);
  }
}
#endif

/* ------------------------------------------------------------------ 2D --- */

NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_trafo_2d_compute)(C *fj,
    const C *g, const R *psij_const0, const R *psij_const1, const R *xj0,
    const R *xj1, const INT n0, const INT n1, const INT m)
{
  nfft_run r0, r1;
  const R *psij0 = psij_const0;
  C s = K(0.0);
  INT p0, l0;

  NFFT_SIMD_ID(nfft_run_init)(&r0, *xj0, n0, m);
  NFFT_SIMD_ID(nfft_run_init)(&r1, *xj1, n1, m);

  for (p0 = 0; p0 < r0.pieces; p0++)
    for (l0 = 0; l0 < r0.len[p0]; l0++, psij0++)
      s += (*psij0) * NFFT_SIMD_ID(nfft_run_cdot)(g + (r0.start[p0] + l0) * n1,
          &r1, psij_const1);

  *fj = s;
}

#ifndef _OPENMP
NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_adjoint_2d_compute_serial)(
    const C *fj, C *g, const R *psij_const0, const R *psij_const1,
    const R *xj0, const R *xj1, const INT n0, const INT n1, const INT m)
{
  nfft_run r0, r1;
  const R *psij0 = psij_const0;
  INT p0, l0;

  NFFT_SIMD_ID(nfft_run_init)(&r0, *xj0, n0, m);
  NFFT_SIMD_ID(nfft_run_init)(&r1, *xj1, n1, m);

  for (p0 = 0; p0 < r0.pieces; p0++)
    for (l0 = 0; l0 < r0.len[p0]; l0++, psij0++)
      NFFT_SIMD_ID(nfft_run_caxpy)(g + (r0.start[p0] + l0) * n1, &r1,
          psij_const1, (*psij0) * (*fj));
}
#endif

#ifdef _OPENMP
/** The rows u0 .. o0 of one node's 2D contribution, empty when o0 < u0. */
NFFT_SIMD_TGT static inline void NFFT_SIMD_ID(nfft_adjoint_2d_block)(C *g,
    const R *psij_const0, const R *psij_const1, const nfft_run *r1, const C f,
    const INT n1, const INT u0, const INT o0, const INT offset_psij)
{
  INT l0;

  for (l0 = 0; l0 <= o0 - u0; l0++)
    NFFT_SIMD_ID(nfft_run_caxpy)(g + (u0 + l0) * n1, r1, psij_const1,
        psij_const0[offset_psij + l0] * f);
}

/**
 * Adjoint NFFT for two-dimensional case updating only a specified range of
 * vector g.
 *
 * \author Toni Volkmer
 */
NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_adjoint_2d_compute_omp_blockwise)(
    const C f, C *g, const R *psij_const0, const R *psij_const1, const R *xj0,
    const R *xj1, const INT n0, const INT n1, const INT m, const INT my_u0,
    const INT my_o0)
{
  INT ar_u0, ar_o0;
  nfft_run r1;

  uo2(&ar_u0, &ar_o0, *xj0, n0, m);
  NFFT_SIMD_ID(nfft_run_init)(&r1, *xj1, n1, m);

  if (ar_u0 < ar_o0)
  {
    const INT u0 = MAX(my_u0, ar_u0);
    const INT o0 = MIN(my_o0, ar_o0);

    NFFT_SIMD_ID(nfft_adjoint_2d_block)(g, psij_const0, psij_const1, &r1, f,
        n1, u0, o0, u0 - ar_u0);
  }
  else
  {
    const INT u0 = MAX(my_u0, ar_u0);
    const INT offset_psij = u0 - ar_u0;

    NFFT_SIMD_ID(nfft_adjoint_2d_block)(g, psij_const0, psij_const1, &r1, f,
        n1, u0, my_o0, offset_psij);
    NFFT_SIMD_ID(nfft_adjoint_2d_block)(g, psij_const0, psij_const1, &r1, f,
        n1, my_u0, MIN(my_o0, ar_o0), offset_psij + my_u0 - ar_u0 + n0);
  }
}
#endif

/* ------------------------------------------------------------------ 3D --- */

NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_trafo_3d_compute)(C *fj,
    const C *g, const R *psij_const0, const R *psij_const1,
    const R *psij_const2, const R *xj0, const R *xj1, const R *xj2,
    const INT n0, const INT n1, const INT n2, const INT m)
{
  nfft_run r0, r1, r2;
  const R *psij0 = psij_const0;
  C s = K(0.0);
  INT p0, l0, p1, l1;

  NFFT_SIMD_ID(nfft_run_init)(&r0, *xj0, n0, m);
  NFFT_SIMD_ID(nfft_run_init)(&r1, *xj1, n1, m);
  NFFT_SIMD_ID(nfft_run_init)(&r2, *xj2, n2, m);

  for (p0 = 0; p0 < r0.pieces; p0++)
    for (l0 = 0; l0 < r0.len[p0]; l0++, psij0++)
    {
      const R *psij1 = psij_const1;
      const INT i0 = (r0.start[p0] + l0) * n1;
      C s1 = K(0.0);

      for (p1 = 0; p1 < r1.pieces; p1++)
        for (l1 = 0; l1 < r1.len[p1]; l1++, psij1++)
          s1 += (*psij1) * NFFT_SIMD_ID(nfft_run_cdot)(
              g + (i0 + r1.start[p1] + l1) * n2, &r2, psij_const2);

      s += (*psij0) * s1;
    }

  *fj = s;
}

#ifndef _OPENMP
NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_adjoint_3d_compute_serial)(
    const C *fj, C *g, const R *psij_const0, const R *psij_const1,
    const R *psij_const2, const R *xj0, const R *xj1, const R *xj2,
    const INT n0, const INT n1, const INT n2, const INT m)
{
  nfft_run r0, r1, r2;
  const R *psij0 = psij_const0;
  INT p0, l0, p1, l1;

  NFFT_SIMD_ID(nfft_run_init)(&r0, *xj0, n0, m);
  NFFT_SIMD_ID(nfft_run_init)(&r1, *xj1, n1, m);
  NFFT_SIMD_ID(nfft_run_init)(&r2, *xj2, n2, m);

  for (p0 = 0; p0 < r0.pieces; p0++)
    for (l0 = 0; l0 < r0.len[p0]; l0++, psij0++)
    {
      const R *psij1 = psij_const1;
      const INT i0 = (r0.start[p0] + l0) * n1;
      const C v0 = (*psij0) * (*fj);

      for (p1 = 0; p1 < r1.pieces; p1++)
        for (l1 = 0; l1 < r1.len[p1]; l1++, psij1++)
          NFFT_SIMD_ID(nfft_run_caxpy)(g + (i0 + r1.start[p1] + l1) * n2, &r2,
              psij_const2, (*psij1) * v0);
    }
}
#endif

#ifdef _OPENMP
/** The rows u0 .. o0 of one node's 3D contribution, empty when o0 < u0. */
NFFT_SIMD_TGT static inline void NFFT_SIMD_ID(nfft_adjoint_3d_block)(C *g,
    const R *psij_const0, const R *psij_const1, const R *psij_const2,
    const nfft_run *r1, const nfft_run *r2, const C f, const INT n1,
    const INT n2, const INT u0, const INT o0, const INT offset_psij)
{
  INT l0, p1, l1;

  for (l0 = 0; l0 <= o0 - u0; l0++)
  {
    const R *psij1 = psij_const1;
    const INT i0 = (u0 + l0) * n1;
    const C v0 = psij_const0[offset_psij + l0] * f;

    for (p1 = 0; p1 < r1->pieces; p1++)
      for (l1 = 0; l1 < r1->len[p1]; l1++, psij1++)
        NFFT_SIMD_ID(nfft_run_caxpy)(g + (i0 + r1->start[p1] + l1) * n2, r2,
            psij_const2, (*psij1) * v0);
  }
}

/**
 * Adjoint NFFT for three-dimensional case updating only a specified range of
 * vector g.
 *
 * \author Toni Volkmer
 */
NFFT_SIMD_TGT static void NFFT_SIMD_ID(nfft_adjoint_3d_compute_omp_blockwise)(
    const C f, C *g, const R *psij_const0, const R *psij_const1,
    const R *psij_const2, const R *xj0, const R *xj1, const R *xj2,
    const INT n0, const INT n1, const INT n2, const INT m, const INT my_u0,
    const INT my_o0)
{
  INT ar_u0, ar_o0;
  nfft_run r1, r2;

  uo2(&ar_u0, &ar_o0, *xj0, n0, m);
  NFFT_SIMD_ID(nfft_run_init)(&r1, *xj1, n1, m);
  NFFT_SIMD_ID(nfft_run_init)(&r2, *xj2, n2, m);

  if (ar_u0 < ar_o0)
  {
    const INT u0 = MAX(my_u0, ar_u0);
    const INT o0 = MIN(my_o0, ar_o0);

    NFFT_SIMD_ID(nfft_adjoint_3d_block)(g, psij_const0, psij_const1,
        psij_const2, &r1, &r2, f, n1, n2, u0, o0, u0 - ar_u0);
  }
  else
  {
    const INT u0 = MAX(my_u0, ar_u0);
    const INT offset_psij = u0 - ar_u0;

    NFFT_SIMD_ID(nfft_adjoint_3d_block)(g, psij_const0, psij_const1,
        psij_const2, &r1, &r2, f, n1, n2, u0, my_o0, offset_psij);
    NFFT_SIMD_ID(nfft_adjoint_3d_block)(g, psij_const0, psij_const1,
        psij_const2, &r1, &r2, f, n1, n2, my_u0, MIN(my_o0, ar_o0),
        offset_psij + my_u0 - ar_u0 + n0);
  }
}
#endif
