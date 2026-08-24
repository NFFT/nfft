/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* Choosing the window cut-off m for the one-dimensional Kaiser-Bessel NFFT.
 *
 * Two effects oppose each other, both set by sigma = n/N alone. Writing
 * u = 1 - 1/sigma:
 *
 *   b = 2*pi*(1 - 1/(2*sigma))   Kaiser-Bessel shape parameter
 *   D = 2*pi*sqrt(u)             truncation decays as exp(-D*m)
 *   A = b - D                    deconvolution amplifies roundoff as exp(A*m)
 *
 * The window is cut to 2m+2 samples, and the deconvolution divides by its
 * Fourier transform, smallest at the band edge. Since b^2 - (pi/sigma)^2 = D^2
 * we have b > D, so A > 0 at every finite sigma and the error always falls,
 * bottoms out, then climbs again. Where A is large the minimum sits near
 * eps^(D/b); by sigma = 2 the amplification is weak enough that the floor is
 * back to a small multiple of eps in every precision. The cost is steep below
 * that, and steepest for the widest mantissa: measured floors at sigma = 3/2
 * are 21*eps in double but 6400*eps in long double. Hence the 5/4 floor
 * enforced below, where double already bottoms out near 1000*eps.
 *
 * Only the prefactors are fitted, over sigma in [5/4, 4], N in
 * {64, 256, 1024}, m in [1, 40] and all three precisions, then raised until
 * the formula dominates every measurement. It is an upper bound, so it errs
 * towards a larger m and towards calling a just-reachable goal unreachable.
 * See .scratch/sigma-m-study/ for the driver, the data and the fit.
 */

#include "nfft3.h"
#include "infft.h"

/* Prefactors of the two branches, from the fit. */
typedef struct {
  R a, p, r; /* truncation: a * u^p * m^r * exp(-D*m), u = 1 - 1/sigma */
  R c, alpha, q; /* roundoff: c * eps * u^q * exp(alpha*A*m) */
} tune_coeff;

static const tune_coeff tune_forward = {K(3.0431),  K(0.902205), K(-0.0183106),
                                        K(68.9787), K(0.967705), K(1.67263)};

static const tune_coeff tune_adjoint = {K(2.03698), K(0.234342), K(0.401585),
                                        K(33.5633), K(0.994013), K(1.15639)};

/* Largest cut-off ever returned, and the top of the fitted range. The error
 * curve bottoms out well below this for every sigma and precision. */
#define TUNE_M_MAX 40

static R tune_error(const tune_coeff *k, R sigma, R eps, int m)
{
  const R u = K(1.0) - K(1.0) / sigma;
  const R d = K2PI * SQRT(u);
  const R b = K2PI * (K(1.0) - K(1.0) / (K(2.0) * sigma));
  const R mm = (R)m;
  const R truncation = k->a * POW(u, k->p) * POW(mm, k->r) * EXP(-d * mm);
  const R roundoff = k->c * eps * POW(u, k->q) * EXP(k->alpha * (b - d) * mm);

  return truncation + roundoff;
}

/* Scan m over 1..m_max. Sets *best_m/*best_e to the argmin and returns the
 * smallest m at or below goal, or 0 when none is. */
static int tune_scan(const tune_coeff *k, R sigma, R eps, int m_max, R goal,
                     int *best_m, R *best_e)
{
  int i, hit = 0;

  *best_m = 1;
  *best_e = tune_error(k, sigma, eps, 1);
  for (i = 1; i <= m_max; i++) {
    const R e = tune_error(k, sigma, eps, i);
    if (e < *best_e) {
      *best_e = e;
      *best_m = i;
    }
    if (!hit && e <= goal)
      hit = i;
  }
  return hit;
}

/* Largest cut-off usable on an oversampled grid of size n: the window spans
 * 2m+2 samples, and the fit stops at TUNE_M_MAX. */
static int tune_m_max(NFFT_INT n)
{
  int m_max = (int)(n / 2 - 1);
  return m_max > TUNE_M_MAX ? TUNE_M_MAX : m_max;
}

int Y(tune)(NFFT_INT N, NFFT_INT n, int adjoint, R goal, int *m, R *attained)
{
  const tune_coeff *k = adjoint ? &tune_adjoint : &tune_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  int m_max, best_m, hit;
  R best_e;

  if (m == 0 || N < (NFFT_INT)1 || n <= N || !(goal > K(0.0)))
    return -1;

  /* sigma = n/N must be at least 5/4, tested in exact integer arithmetic.
   * Below that the deconvolution amplifies roundoff nearly as fast as the
   * window truncation decays, the attainable accuracy collapses, and the
   * formula here is outside the range it was fitted and checked on. */
  if ((NFFT_INT)4 * n < (NFFT_INT)5 * N)
    return -1;

  m_max = tune_m_max(n);
  if (m_max < 1)
    return -1;

  hit = tune_scan(k, (R)n / (R)N, eps, m_max, goal, &best_m, &best_e);

  if (hit) {
    *m = hit;
    if (attained)
      *attained = tune_error(k, (R)n / (R)N, eps, hit);
    return 1;
  }

  *m = best_m;
  if (attained)
    *attained = best_e;
  return 0;
}

/* Best accuracy any m can reach on this geometry. */
static R tune_floor(const tune_coeff *k, NFFT_INT N, NFFT_INT n, R eps)
{
  int best_m;
  R best_e;
  const int m_max = tune_m_max(n);

  if (m_max < 1)
    return K(-1.0);
  tune_scan(k, (R)n / (R)N, eps, m_max, K(0.0), &best_m, &best_e);
  return best_e;
}

/* Bottom and top of the searchable band: even n with sigma in [5/4, 4]. */
static NFFT_INT tune_band_lo(NFFT_INT N)
{
  NFFT_INT lo = (((NFFT_INT)5 * N + (NFFT_INT)3) / (NFFT_INT)4 + (NFFT_INT)1)
                / (NFFT_INT)2 * (NFFT_INT)2;
  if (lo <= N)
    lo = (N + (NFFT_INT)2) / (NFFT_INT)2 * (NFFT_INT)2;
  return lo;
}

/* Smallest even n in the band whose floor reaches goal, or 0 if none does.
 * The floor falls as sigma grows, so bisect on "this n reaches the goal",
 * keeping hi always reaching it. */
static NFFT_INT tune_smallest_n(const tune_coeff *k, NFFT_INT N, R goal, R eps)
{
  NFFT_INT lo = tune_band_lo(N);
  NFFT_INT hi = (NFFT_INT)4 * N;

  if (hi < lo)
    hi = lo;
  if (tune_floor(k, N, hi, eps) > goal)
    return 0;

  while (hi - lo > (NFFT_INT)2) {
    const NFFT_INT mid =
         (lo + (hi - lo) / (NFFT_INT)2) / (NFFT_INT)2 * (NFFT_INT)2;
    if (mid <= lo || mid >= hi)
      break;
    if (tune_floor(k, N, mid, eps) <= goal)
      hi = mid;
    else
      lo = mid;
  }
  if (tune_floor(k, N, lo, eps) <= goal)
    hi = lo;
  return hi;
}

int Y(tune_sigma)(NFFT_INT N, int adjoint, R goal, NFFT_INT *n, R *attained)
{
  const tune_coeff *k = adjoint ? &tune_adjoint : &tune_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  NFFT_INT hit;

  if (n == 0 || N < (NFFT_INT)1 || !(goal > K(0.0)))
    return -1;

  hit = tune_smallest_n(k, N, goal, eps);
  if (hit == 0) {
    const NFFT_INT top = (NFFT_INT)4 * N;
    *n = top;
    if (attained)
      *attained = tune_floor(k, N, top, eps);
    return 0;
  }

  *n = hit;
  if (attained)
    *attained = tune_floor(k, N, hit, eps);
  return 1;
}

/* Smallest even 5-smooth number (2^a*3^b*5^c, a >= 1) that is at least kk.
 * Enumerates every 3^b*5^c up to the bound and rounds each up by powers of
 * two, so it touches each candidate once: O(log^3 kk), no upward scan. */
static NFFT_INT tune_next_smooth(NFFT_INT kk)
{
  const NFFT_INT limit = (NFFT_INT)1 << (8 * (int)sizeof(NFFT_INT) - 3);
  NFFT_INT best = 0, p5, p35, v;

  if (kk < (NFFT_INT)2)
    return (NFFT_INT)2;

  /* The answer is below 2*kk, so its odd part 3^b*5^c is at most kk. */
  for (p5 = 1; p5 <= kk && p5 < limit; p5 *= 5)
    for (p35 = p5; p35 <= kk && p35 < limit; p35 *= 3) {
      v = p35 * (NFFT_INT)2; /* a >= 1 keeps n even */
      while (v < kk && v < limit)
        v *= (NFFT_INT)2;
      if (v >= kk && (best == 0 || v < best))
        best = v;
    }
  return best;
}

int Y(tune_plan)(NFFT_INT N, int adjoint, R goal, NFFT_INT *n, int *m,
                 R *attained)
{
  const tune_coeff *k = adjoint ? &tune_adjoint : &tune_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  R cap, a;
  NFFT_INT n0, nn;
  int tries;

  if (n == 0 || m == 0 || N < (NFFT_INT)1 || !(goal > K(0.0)))
    return -1;

  /* Cap the goal at what the band can actually deliver: the floor at the
   * widest oversampling on offer. Asking for less than this is asking for
   * accuracy no cut-off and no grid in range can reach. */
  cap = tune_floor(k, N, (NFFT_INT)4 * N, eps);
  a = goal > cap ? goal : cap;

  n0 = tune_smallest_n(k, N, a, eps);
  if (n0 == 0)
    n0 = (NFFT_INT)4 * N; /* a == cap: the top of the band is the answer */

  /* Round up to a size the FFT likes. That only raises sigma, which never
   * hurts -- but the error is not quite monotone in sigma at small m, so the
   * cut-off is re-derived at the size actually chosen and the answer checked
   * rather than assumed. */
  nn = tune_next_smooth(n0);
  if (nn == 0)
    nn = n0;

  for (tries = 0; tries < 8; tries++) {
    const int m_max = tune_m_max(nn);
    int best_m, hit;
    R best_e;

    if (m_max >= 1) {
      hit = tune_scan(k, (R)nn / (R)N, eps, m_max, a, &best_m, &best_e);
      if (hit) {
        const R got = tune_error(k, (R)nn / (R)N, eps, hit);
        *n = nn;
        *m = hit;
        if (attained)
          *attained = got;
        return got <= goal ? 1 : 0;
      }
    }
    nn = tune_next_smooth(nn + (NFFT_INT)1);
    if (nn == 0)
      break;
  }

  /* No smooth size cleared it: fall back to the unrounded n, whose floor was
   * established above. */
  {
    int best_m;
    R best_e;
    const int m_max = tune_m_max(n0);
    if (m_max < 1)
      return -1;
    tune_scan(k, (R)n0 / (R)N, eps, m_max, K(0.0), &best_m, &best_e);
    *n = n0;
    *m = best_m;
    if (attained)
      *attained = best_e;
    return best_e <= goal ? 1 : 0;
  }
}
