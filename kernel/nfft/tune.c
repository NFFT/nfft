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
 * Only the prefactors are fitted, over sigma in [5/4, 4], N in [32, 4096],
 * M/N in [1/4, 8], m in [1, 32] and all three precisions, then raised until
 * the formula dominates every measurement. It is an upper bound, so it errs
 * towards a larger m and towards calling a just-reachable goal unreachable.
 * See .scratch/sigma-m-study/ for the drivers, the data and the fit.
 */

#include <string.h>

#include "nfft3.h"
#include "infft.h"

/* Prefactors of the two branches, from the fit in gfit.py. */
typedef struct {
  /* truncation: a * u^p * m^r * exp(-D*m) * N^tn * M^tm, u = 1 - 1/sigma */
  R a, p, r, tn, tm;
  /* roundoff: c * eps * u^q * exp(alpha*A*m) * N^un * M^um */
  R c, alpha, q, un, um;
} tune_coeff;

static const tune_coeff tune_forward = {
     K(47.1246), K(0.93096),  K(-0.00821491), K(-0.577425),  K(0.0480177),
     K(216.943), K(0.951476), K(2.25144),     K(-0.0252356), K(0.11049)};

static const tune_coeff tune_adjoint = {
     K(24.7348), K(0.312453), K(0.434428), K(0.0374448), K(-0.506233),
     K(321.955), K(0.989474), K(1.48935),  K(0.415192),  K(-0.427671)};

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

int Y(tune)(NFFT_INT N, NFFT_INT n, NFFT_INT M, int adjoint, R goal, int *m,
            R *attained)
{
  const tune_coeff *k = adjoint ? &tune_adjoint : &tune_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  int m_max, best_m, hit;
  R best_e;

  if (m == 0 || N < (NFFT_INT)1 || M < (NFFT_INT)1 || n <= N
      || !(goal > K(0.0)))
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

  hit = tune_scan(k, (R)n / (R)N, eps, m_max, goal, N, M, &best_m, &best_e);

  if (hit) {
    *m = hit;
    if (attained)
      *attained = tune_error(k, (R)n / (R)N, eps, hit, N, M);
    return 1;
  }

  *m = best_m;
  if (attained)
    *attained = best_e;
  return 0;
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
static NFFT_INT tune_smallest_n(const tune_coeff *k, NFFT_INT N, NFFT_INT M,
                                R goal, R eps)
{
  NFFT_INT lo = tune_band_lo(N);
  NFFT_INT hi = (NFFT_INT)4 * N;

  if (hi < lo)
    hi = lo;
  if (tune_floor(k, N, hi, M, eps) > goal)
    return 0;

  while (hi - lo > (NFFT_INT)2) {
    const NFFT_INT mid =
         (lo + (hi - lo) / (NFFT_INT)2) / (NFFT_INT)2 * (NFFT_INT)2;
    if (mid <= lo || mid >= hi)
      break;
    if (tune_floor(k, N, mid, M, eps) <= goal)
      hi = mid;
    else
      lo = mid;
  }
  if (tune_floor(k, N, lo, M, eps) <= goal)
    hi = lo;
  return hi;
}

int Y(tune_sigma)(NFFT_INT N, NFFT_INT M, int adjoint, R goal, NFFT_INT *n,
                  R *attained)
{
  const tune_coeff *k = adjoint ? &tune_adjoint : &tune_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  NFFT_INT hit;

  if (n == 0 || N < (NFFT_INT)1 || M < (NFFT_INT)1 || !(goal > K(0.0)))
    return -1;

  hit = tune_smallest_n(k, N, M, goal, eps);
  if (hit == 0) {
    const NFFT_INT top = (NFFT_INT)4 * N;
    *n = top;
    if (attained)
      *attained = tune_floor(k, N, top, M, eps);
    return 0;
  }

  *n = hit;
  if (attained)
    *attained = tune_floor(k, N, hit, M, eps);
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

/* Number of factors of two.
 *
 * Two 5-smooth grids of nearly the same size do not cost the same to
 * transform: every FFT implementation optimises its radix-2 codelets hardest,
 * and measured per-point cost puts other 5-smooth sizes 5 to 40 % above a
 * power of two -- n = 480 loses to n = 512 despite being the smaller grid.
 * No operation count reproduces that ordering, since a radix-3 or radix-5
 * stage does more work per stage and counts out cheaper per point than
 * radix-2. So this is a tie-break and not a term in the cost: it applies only
 * between grids that need the same cut-off, where the node work is identical,
 * and it never overrides a grid more than TUNE_FFT_TIE cheaper.
 *
 * The spread is far wider in float, where n = 432 measures 1.7x the time of
 * n = 512 though it is the smaller grid. Widening the window to catch those
 * does not pay: at 1.25 every shape's median speed falls, because the grids
 * it then buys are dearer more often than they are faster. Choosing a grid on
 * measured time is what measured planning is for.
 */
static int tune_twos(NFFT_INT n)
{
  int k = 0;

  while (n % (NFFT_INT)2 == (NFFT_INT)0) {
    n /= (NFFT_INT)2;
    k++;
  }
  return k;
}

/* How much more FFT work the tie-break may accept for a richer power of two.
 * Set from the measured spread above, on the safe side of it. */
#define TUNE_FFT_TIE K(1.1)

int Y(tune_plan)(NFFT_INT N, NFFT_INT M, int adjoint, R goal, NFFT_INT *n,
                 int *m, R *attained)
{
  const tune_coeff *k = adjoint ? &tune_adjoint : &tune_forward;
  const R eps = Y(float_property)(NFFT_EPSILON);
  const NFFT_INT hi = (NFFT_INT)4 * N;
  R cap, a, best_cost = K(0.0), best_e = K(0.0);
  NFFT_INT n0, nn, best_n = 0;
  int best_m = 0, tries;

  if (n == 0 || m == 0 || N < (NFFT_INT)1 || M < (NFFT_INT)1
      || !(goal > K(0.0)))
    return -1;

  /* Cap the goal at what the band can actually deliver: the floor at the
   * widest oversampling on offer. Asking for less than this is asking for
   * accuracy no cut-off and no grid in range can reach. */
  cap = tune_floor(k, N, hi, M, eps);
  a = goal > cap ? goal : cap;

  /* Walk every even 5-smooth size in the band, only tens of them, take the
   * smallest sufficient cut-off at each and keep the cheapest pair. The
   * legacy size 2*next_power_of_2(N) is a power of two in [2N, 4N), so it is
   * always among the candidates and the answer never rates worse than it. */
  for (nn = tune_next_smooth(tune_band_lo(N)); nn != 0 && nn <= hi;
       nn = tune_next_smooth(nn + (NFFT_INT)1)) {
    const int m_max = tune_m_max(nn);
    int cand_m, hit, take;
    R cand_e, cost;

    if (m_max < 1)
      continue;
    hit = tune_scan(k, (R)nn / (R)N, eps, m_max, a, N, M, &cand_m, &cand_e);
    if (!hit)
      continue;

    cost = tune_cost(nn, hit, M);
    take = best_n == 0 || cost < best_cost;
    /* Same cut-off, near-equal FFT work, more factors of two: faster in
     * practice though the model rates it dearer. M cancels here, so the
     * choice stays monotone in the node count. */
    if (!take && hit == best_m && tune_twos(nn) > tune_twos(best_n)
        && tune_fft_cost(nn) <= TUNE_FFT_TIE * tune_fft_cost(best_n))
      take = 1;
    if (take) {
      best_cost = cost;
      best_n = nn;
      best_m = hit;
      best_e = tune_error(k, (R)nn / (R)N, eps, hit, N, M);
    }
  }

  if (best_n != 0) {
    *n = best_n;
    *m = best_m;
    if (attained)
      *attained = best_e;
    return best_e <= goal ? 1 : 0;
  }

  /* Nothing in the band cleared the capped goal, so the floor is only reached
   * past 4*N. Keep walking: more oversampling is never worse, and the first
   * size that clears it is the cheapest of those that do. */
  n0 = tune_smallest_n(k, N, M, a, eps);
  if (n0 == 0)
    n0 = hi;
  nn = tune_next_smooth(n0 > hi ? n0 : hi + (NFFT_INT)1);
  for (tries = 0; tries < 8 && nn != 0; tries++) {
    const int m_max = tune_m_max(nn);
    int cand_m, hit;
    R cand_e;

    if (m_max >= 1) {
      hit = tune_scan(k, (R)nn / (R)N, eps, m_max, a, N, M, &cand_m, &cand_e);
      if (hit) {
        const R got = tune_error(k, (R)nn / (R)N, eps, hit, N, M);
        *n = nn;
        *m = hit;
        if (attained)
          *attained = got;
        return got <= goal ? 1 : 0;
      }
    }
    nn = tune_next_smooth(nn + (NFFT_INT)1);
  }

  /* No smooth size cleared it: fall back to the unrounded n, whose floor was
   * established above. */
  {
    int fb_m;
    R fb_e;
    const int m_max = tune_m_max(n0);
    if (m_max < 1)
      return -1;
    tune_scan(k, (R)n0 / (R)N, eps, m_max, K(0.0), N, M, &fb_m, &fb_e);
    *n = n0;
    *m = fb_m;
    if (attained)
      *attained = fb_e;
    return fb_e <= goal ? 1 : 0;
  }
}

/* Measured refinement of the cut-off.
 *
 * The model is an upper envelope over every geometry it was fitted on, while
 * the accuracy that matters is that of the caller's own nodes, so the model
 * hands out about one cut-off more than the problem in front of it needs.
 * That gap cannot be closed by a formula -- at a fixed geometry the error
 * still varies by a factor of a few across node draws, which no a priori
 * expression can see. It can be closed by measuring.
 *
 * This runs the transform against the planner's own direct NDFT on the
 * caller's nodes and steps the cut-off down while the goal still holds. The
 * NDFT does not depend on n or m, so it is computed once: the refinement
 * costs one O(N*M) reference plus one transform per cut-off tried.
 *
 * Opt-in for that reason. It pays when the plan is reused, and pays most
 * where the node convolution dominates the cost -- which is exactly where an
 * extra cut-off is expensive. It measures arithmetic, not time, so it depends
 * on the machine no more than the transform itself does.
 *
 * One random probe vector is used, the same measure the error model was
 * fitted to: max output error over the l_1 norm of the input.
 */
/* Probes per candidate, and the headroom the worst of them must show.
 *
 * The nodes are the caller's, so the only thing left varying is the input
 * data -- and it varies a lot: at a fixed geometry the measured error spans a
 * median 1.55x and up to 6x over 24 random inputs. A single probe shaved to
 * the goal therefore misses it on the next vector, measurably so. The worst
 * of eight probes with a factor of two in hand does not: over 3840 fresh
 * draws on 96 refined geometries, none exceeded the goal.
 *
 * The margin is cheaper than it looks because m is quantised. Dropping one
 * cut-off multiplies the error by exp(D), which is 30 to 90 across the band,
 * so a cut-off is only ever wasted when the headroom exceeds that -- and the
 * margin merely raises the bar from exp(D) to 2*exp(D).
 */
#define TUNE_REFINE_PROBES 8
#define TUNE_REFINE_MARGIN K(2.0)

static R tune_probe(NFFT_INT N, NFFT_INT n, NFFT_INT M, int m, int adjoint,
                    const R *x, const C *in, const C *ref, C *out, C *fh,
                    C *ff, NFFT_INT in_len, NFFT_INT out_len)
{
  Y(plan_ng) * fast;
  R worst = K(0.0);
  int i;
  NFFT_INT j;

  fast = Y(plan_ng_guru)(1, &N, 0, &n, M, m, NFFT_WINDOW_KAISER_BESSEL, (R *)x,
                         (void *)fh, (void *)ff, 0u,
                         NFFT_ESTIMATE | NFFT_NO_DIRECT);
  if (!fast)
    return K(-1.0);
  Y(precompute)(fast);

  for (i = 0; i < TUNE_REFINE_PROBES; i++) {
    const C *pin = in + (size_t)i * (size_t)in_len;
    const C *pref = ref + (size_t)i * (size_t)out_len;
    R num = K(0.0), den = K(0.0);

    if (adjoint)
      Y(execute_adjoint_on)(fast, (void *)out, (void *)pin);
    else
      Y(execute_on)(fast, (void *)pin, (void *)out);

    for (j = 0; j < out_len; j++) {
      const R e = CABS(out[j] - pref[j]);
      if (e > num)
        num = e;
    }
    for (j = 0; j < in_len; j++)
      den += CABS(pin[j]);
    if (den > K(0.0) && num / den > worst)
      worst = num / den;
  }

  Y(plan_ng_destroy)(fast);
  return worst;
}

int Y(tune_refine)(NFFT_INT N, NFFT_INT M, int adjoint, R goal, const R *x,
                   NFFT_INT n, int *m, R *attained)
{
  const NFFT_INT in_len = adjoint ? M : N;
  const NFFT_INT out_len = adjoint ? N : M;
  const int m_max = tune_m_max(n);
  C *in, *ref, *out, *fh, *ff;
  Y(plan_ng) * direct;
  int cur, best, i;
  R best_e = K(-1.0), e0;

  if (x == 0 || m == 0 || N < (NFFT_INT)1 || M < (NFFT_INT)1 || n <= N
      || !(goal > K(0.0)) || *m < 1 || *m > m_max)
    return -1;
  if ((NFFT_INT)4 * n < (NFFT_INT)5 * N)
    return -1;

  in = (C *)Y(malloc)((size_t)TUNE_REFINE_PROBES * (size_t)in_len * sizeof(C));
  ref = (C *)Y(malloc)((size_t)TUNE_REFINE_PROBES * (size_t)out_len
                       * sizeof(C));
  out = (C *)Y(malloc)((size_t)out_len * sizeof(C));
  fh = (C *)Y(malloc)((size_t)N * sizeof(C));
  ff = (C *)Y(malloc)((size_t)M * sizeof(C));
  if (!in || !ref || !out || !fh || !ff) {
    Y(free)(in);
    Y(free)(ref);
    Y(free)(out);
    Y(free)(fh);
    Y(free)(ff);
    return -1;
  }
  memset(fh, 0, (size_t)N * sizeof(C));
  memset(ff, 0, (size_t)M * sizeof(C));

  /* The exact transform of each probe, once: it depends on the nodes and the
   * data, not on the cut-off. */
  direct = Y(plan_ng_guru)(1, &N, 0, &n, M, *m, NFFT_WINDOW_KAISER_BESSEL,
                           (R *)x, (void *)fh, (void *)ff, 0u,
                           NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
  if (!direct) {
    Y(free)(in);
    Y(free)(ref);
    Y(free)(out);
    Y(free)(fh);
    Y(free)(ff);
    return -1;
  }
  Y(precompute)(direct);
  Y(srand48)(1723);
  for (i = 0; i < TUNE_REFINE_PROBES; i++) {
    C *pin = in + (size_t)i * (size_t)in_len;
    C *pref = ref + (size_t)i * (size_t)out_len;
    Y(vrand_unit_complex)((void *)pin, in_len);
    if (adjoint)
      Y(execute_adjoint_on)(direct, (void *)pref, (void *)pin);
    else
      Y(execute_on)(direct, (void *)pin, (void *)pref);
  }
  Y(plan_ng_destroy)(direct);

  /* The model is an upper bound, so its cut-off should already measure well
   * inside the goal. If it does not, say so and change nothing. */
  e0 = tune_probe(N, n, M, *m, adjoint, x, in, ref, out, fh, ff, in_len,
                  out_len);
  best = *m;
  best_e = e0;
  if (e0 < K(0.0) || e0 > goal) {
    if (attained && e0 >= K(0.0))
      *attained = e0;
    Y(free)(in);
    Y(free)(ref);
    Y(free)(out);
    Y(free)(fh);
    Y(free)(ff);
    return 0;
  }

  for (cur = *m - 1; cur >= 1; cur--) {
    const R e = tune_probe(N, n, M, cur, adjoint, x, in, ref, out, fh, ff,
                           in_len, out_len);
    if (e < K(0.0) || e * TUNE_REFINE_MARGIN > goal)
      break;
    best = cur;
    best_e = e;
  }

  *m = best;
  if (attained)
    *attained = best_e;

  Y(free)(in);
  Y(free)(ref);
  Y(free)(out);
  Y(free)(fh);
  Y(free)(ff);
  return 1;
}
