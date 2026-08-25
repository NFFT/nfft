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

/* Choosing the oversampling n and the window cut-off m for the
 * one-dimensional Kaiser-Bessel NFFT.
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
 * Both branches also scale with the geometry, because the error measure is
 * relative:
 *
 *   forward   max_j |f_j - s_j| / ||fhat||_1
 *   adjoint   max_k |fhat_k - s_k| / ||f||_1
 *
 * The forward numerator is a max over M nodes of a sum of N terms whose
 * phases cancel, over a denominator that grows like N, so it falls roughly
 * like N^-1/2 and rises slowly with M. The adjoint swaps the two. Roundoff
 * accumulates instead of cancelling, so its exponents differ from the
 * truncation branch's. Hence the N and M powers below: they are the
 * measure's own normalisation, not curve fitting, and a model without them
 * misses the goal outright once M leaves the band it was fitted in.
 *
 * Only the prefactors are fitted, over sigma in [5/4, 8], N in [32, 1024],
 * M/N in [1/4, 8], m in [1, 32] and all three precisions, then raised until
 * the formula dominates every measurement. It is an upper bound, so it errs
 * towards a larger m and towards calling a just-reachable goal unreachable.
 *
 * The fit measures the error on input whose real and imaginary parts are
 * uniform on [0, 1), as Y(vrand_unit_complex) draws it. Centred input
 * measures a forward error up to 2.6x smaller, so an envelope fitted to that
 * does not hold here. See .scratch/sigma-m-study/ for the drivers, the data
 * and the fit.
 */

#include "nfft3.h"
#include "infft.h"

/* Prefactors of the two branches, from the fit in dfit.py. */
typedef struct {
  /* truncation: a * u^p * m^r * exp(-D*m) * N^tn * M^tm, u = 1 - 1/sigma */
  R a, p, r, tn, tm;
  /* roundoff: c * eps * u^q * exp(alpha*A*m) * N^un * M^um */
  R c, alpha, q, un, um;
} tune_coeff;

/* Largest cut-off ever returned, and the top of the fitted range. The error
 * curve bottoms out well below this for every sigma and precision. */
#define TUNE_M_MAX 32

static R tune_error(const tune_coeff *k, R sigma, R eps, int m, NFFT_INT N,
                    NFFT_INT M)
{
  const R u = K(1.0) - K(1.0) / sigma;
  const R d = K2PI * SQRT(u);
  const R b = K2PI * (K(1.0) - K(1.0) / (K(2.0) * sigma));
  const R mm = (R)m;
  const R nn = (R)N, mn = (R)M;
  const R truncation = k->a * POW(u, k->p) * POW(mm, k->r) * EXP(-d * mm)
                       * POW(nn, k->tn) * POW(mn, k->tm);
  const R roundoff = k->c * eps * POW(u, k->q) * EXP(k->alpha * (b - d) * mm)
                     * POW(nn, k->un) * POW(mn, k->um);

  return truncation + roundoff;
}

/* Scan m over 1..m_max. Sets *best_m/*best_e to the argmin and returns the
 * smallest m at or below goal, or 0 when none is. */
static int tune_scan(const tune_coeff *k, R sigma, R eps, int m_max, R goal,
                     NFFT_INT N, NFFT_INT M, int *best_m, R *best_e)
{
  int i, hit = 0;

  *best_m = 1;
  *best_e = tune_error(k, sigma, eps, 1, N, M);
  for (i = 1; i <= m_max; i++) {
    const R e = tune_error(k, sigma, eps, i, N, M);
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

/* Best accuracy any m can reach on this geometry. */
static R tune_floor(const tune_coeff *k, NFFT_INT N, NFFT_INT n, NFFT_INT M,
                    R eps)
{
  int best_m;
  R best_e;
  const int m_max = tune_m_max(n);

  if (m_max < 1)
    return K(-1.0);
  tune_scan(k, (R)n / (R)N, eps, m_max, K(0.0), N, M, &best_m, &best_e);
  return best_e;
}

/* Relative price of one node-convolution sample against one FFT butterfly.
 *
 * An operation count, not a measurement: a complex FFT of size n takes about
 * 5*n*log2(n) real operations, and each of the M*(2m+2) window samples
 * multiplies a complex grid value by a real window value and accumulates,
 * four real operations. A weight fitted to timings would bind the library to
 * one machine's cache and its scatter/gather asymmetry. It would also buy
 * nothing: over the 27600 timings in .scratch/sigma-m-study/costfit.c, 4/5
 * orders 92.6 % of candidate pairs as measured, the best weight there 92.7 %,
 * and the 2/5 the planner's own pcost assumes 90.4 %.
 */
#define TUNE_NODE_WEIGHT (K(4.0) / K(5.0))

/* FFT term of the cost, in the same units. */
static R tune_fft_cost(NFFT_INT n)
{
  const R nn = (R)n;
  return nn * LOG(nn) / LOG(K(2.0));
}

/* Predicted run time of one transform, in arbitrary units. */
static R tune_cost(NFFT_INT n, int m, NFFT_INT M)
{
  return tune_fft_cost(n) + TUNE_NODE_WEIGHT * (R)M * (R)(2 * m + 2);
}

/* The dyadic ladder, n = 2^j * next_power_of_2(N).
 *
 * next_power_of_2(2*N) is 2*next_power_of_2(N) for every N, so the legacy
 * grid is rung 1. With t = next_power_of_2(N)/N in (1, 2] the rungs
 * oversample by sigma = t, 2t and 4t, three disjoint bands, and rung 0 is
 * legal only when t >= 5/4. One coefficient row per rung, selected by j and
 * not by testing sigma, since rung 0 and rung 1 meet at sigma = 2 when t = 1.
 *
 * A row has to dominate its own band only, where tune_forward/tune_adjoint
 * must dominate all of [5/4, 4] at once. Rung 2's alpha is pinned at 1 rather
 * than fitted: A = b - D is 0.056 there falling to 0.012, so exp(alpha*A*m)
 * is flat over m <= 32 and the fit reads noise.
 *
 * Every rung is a power of two, which is why no tie-break appears below.
 *
 * Fitted by .scratch/sigma-m-study/dfit.py, same envelope discipline as the
 * tables at the top of this file.
 */
static const tune_coeff tune_dyadic_forward[3] = {
    {K(60.6176), K(0.928375), K(-0.370109), K(-0.43395), K(0.0568873),
     K(24.4189), K(0.936522), K(0.193903), K(-0.17226), K(0.121651)},
    {K(27.3498), K(2.17547), K(-0.309852), K(-0.217142), K(0.0868104),
     K(21.6553), K(0.593618), K(1.9439), K(0.0856313), K(0.151626)},
    {K(10.2261), K(2.68804), K(-0.212348), K(-0.12615), K(0.0901836),
     K(13.0956), K(1), K(9.04299), K(0.347203), K(0.110038)},
};

static const tune_coeff tune_dyadic_adjoint[3] = {
    {K(22.4855), K(0.208243), K(0.464053), K(0.0432581), K(-0.501897),
     K(173.261), K(1.0011), K(0.310334), K(-0.174676), K(0.00122196)},
    {K(30.2286), K(0.411002), K(0.419399), K(0.0467219), K(-0.503082),
     K(51.0493), K(0.798611), K(0.254312), K(0.330917), K(-0.208085)},
    {K(22.7536), K(0.437151), K(0.368359), K(0.0715423), K(-0.489816),
     K(64.1882), K(1), K(12.1863), K(0.84782), K(-0.411437)},
};

/* Size of rung j, or 0 when that rung is not usable at this bandwidth. */
static NFFT_INT tune_dyadic_size(NFFT_INT N, int j)
{
  NFFT_INT nj;

  if (j < 0 || j > 2)
    return 0;
  nj = Y(next_power_of_2)(N) << j;
  if (nj <= N)
    return 0; /* rung 0 when N is itself a power of two */
  if ((NFFT_INT)4 * nj < (NFFT_INT)5 * N)
    return 0; /* rung 0 when N sits less than a quarter above one */
  if (tune_m_max(nj) < 1)
    return 0;
  return nj;
}

int Y(tune_dyadic_at)(NFFT_INT N, NFFT_INT M, int adjoint, R goal, int j,
                      NFFT_INT *n, int *m, R *attained)
{
  const tune_coeff *band = adjoint ? tune_dyadic_adjoint : tune_dyadic_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  NFFT_INT nj;
  int best_m, hit;
  R best_e;

  if (n == 0 || m == 0 || N < (NFFT_INT)1 || M < (NFFT_INT)1
      || !(goal > K(0.0)))
    return -1;

  nj = tune_dyadic_size(N, j);
  if (nj == (NFFT_INT)0)
    return -1;

  hit = tune_scan(&band[j], (R)nj / (R)N, eps, tune_m_max(nj), goal, N, M,
                  &best_m, &best_e);
  *n = nj;
  if (hit) {
    *m = hit;
    if (attained)
      *attained = tune_error(&band[j], (R)nj / (R)N, eps, hit, N, M);
    return 1;
  }
  *m = best_m;
  if (attained)
    *attained = best_e;
  return 0;
}

int Y(tune_plan_dyadic)(NFFT_INT N, NFFT_INT M, int adjoint, R goal,
                        NFFT_INT *n, int *m, R *attained)
{
  const tune_coeff *band = adjoint ? tune_dyadic_adjoint : tune_dyadic_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  NFFT_INT top, best_n = 0;
  R cap, a, best_cost = K(0.0), best_e = K(0.0);
  int j, best_m = 0;

  if (n == 0 || m == 0 || N < (NFFT_INT)1 || M < (NFFT_INT)1
      || !(goal > K(0.0)))
    return -1;

  /* The top rung always clears both the n > N and the 5/4 tests, so the
   * ladder is never empty and there is always a floor to cap against. */
  top = tune_dyadic_size(N, 2);
  if (top == (NFFT_INT)0)
    return -1;

  cap = tune_floor(&band[2], N, top, M, eps);
  a = goal > cap ? goal : cap;

  for (j = 0; j < 3; j++) {
    NFFT_INT nj = 0;
    int cand_m = 0;
    R cand_e = K(0.0), cost;

    if (Y(tune_dyadic_at)(N, M, adjoint, a, j, &nj, &cand_m, &cand_e) != 1)
      continue;

    cost = tune_cost(nj, cand_m, M);
    if (best_n == (NFFT_INT)0 || cost < best_cost) {
      best_cost = cost;
      best_n = nj;
      best_m = cand_m;
      best_e = cand_e;
    }
  }

  /* The cap is the top rung's own floor, so that rung reaches it and the loop
   * normally returns. It can still come back empty when the two are separated
   * by rounding, since -ffast-math lets one expression differ in the last bit
   * between call sites. Fall back to the top rung's best. */
  if (best_n == (NFFT_INT)0) {
    tune_scan(&band[2], (R)top / (R)N, eps, tune_m_max(top), K(0.0), N, M,
              &best_m, &best_e);
    best_n = top;
  }

  *n = best_n;
  *m = best_m;
  if (attained)
    *attained = best_e;
  return best_e <= goal ? 1 : 0;
}
