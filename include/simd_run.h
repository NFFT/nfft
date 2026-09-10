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

/*! \file simd_run.h
 *  \brief The two vector primitives every window run reduces to, one
 *  instantiation per instruction set.
 *
 *  This header has no include guard on purpose: include it once per variant,
 *  with NFFT_SIMD_VARIANT set to one of the NFFT_SIMD_* ids from simd.h. Each
 *  inclusion defines
 *
 *      C    <name>_cdot (const R *psi, const C *g, INT len);
 *      void <name>_caxpy(C *g, const R *psi, C f, INT len);
 *
 *  where <name> is NFFT_SIMD_ID(nfft_simd) -- i.e. nfft_simd_cdot_avx2 and so
 *  on. cdot returns sum_l psi[l] * g[l] and caxpy adds psi[l] * f to g[l]; the
 *  two directions of the NFFT convolution step, both with a *real* weight
 *  vector against a complex data vector.
 *
 *  That is what makes the vector form simple: because psi is real, g is just
 *  2*len reals and each psi value applies to an adjacent pair of them. So an
 *  instruction set needs to supply little more than load, store, add and
 *  multiply-add; the three operations that are not that -- duplicating the
 *  weights into lane pairs, folding an accumulator down to one complex number,
 *  and spreading one complex number across the lanes -- are nfft_simd_dup2(),
 *  nfft_simd_reduce() and nfft_simd_set2() below.
 *
 *  NFFT_SIMD_ID, NFFT_SIMD_SFX and NFFT_SIMD_TGT stay defined afterwards so
 *  that the kernel bodies built on these primitives can be named and
 *  attributed the same way; the vector ops themselves are undefined again at
 *  the end, and the next inclusion clears the rest. Callers undefine
 *  NFFT_SIMD_VARIANT before that next inclusion.
 *
 *  Reading a C object through R lvalues assumes the layout C99 6.2.5 gives
 *  complex types, two adjacent reals; the OpenMP atomic paths in kernel/nfft
 *  already rely on it.
 */

#ifndef NFFT_SIMD_VARIANT
#error "define NFFT_SIMD_VARIANT before including simd_run.h"
#endif

#undef NFFT_SIMD_CAT_
#undef NFFT_SIMD_CAT
#undef NFFT_SIMD_ID
#undef NFFT_SIMD_SFX
#undef NFFT_SIMD_TGT
#undef NFFT_V
#undef NFFT_VLEN
#undef NFFT_VZERO
#undef NFFT_VLOAD
#undef NFFT_VSTORE
#undef NFFT_VADD
#undef NFFT_VFMA

#define NFFT_SIMD_CAT_(a, b) a ## b
#define NFFT_SIMD_CAT(a, b) NFFT_SIMD_CAT_(a, b)
#define NFFT_SIMD_ID(name) NFFT_SIMD_CAT(name, NFFT_SIMD_SFX)

#if NFFT_SIMD_VARIANT == NFFT_SIMD_SCALAR

#define NFFT_SIMD_SFX _scalar
#define NFFT_SIMD_TGT

#elif NFFT_SIMD_VARIANT == NFFT_SIMD_SSE2

#include <emmintrin.h>
#define NFFT_SIMD_SFX _sse2
/* Harmless where SSE2 is already the baseline (x86-64, or -msse2 on i386):
 * target options add to the command line, they do not replace it. */
#define NFFT_SIMD_TGT __attribute__((target("sse2")))
#if defined(NFFT_SINGLE)
#define NFFT_V __m128
#define NFFT_VLEN 4
#define NFFT_VZERO() _mm_setzero_ps()
#define NFFT_VLOAD(p) _mm_loadu_ps(p)
#define NFFT_VSTORE(p, v) _mm_storeu_ps((p), (v))
#define NFFT_VADD(a, b) _mm_add_ps((a), (b))
#define NFFT_VFMA(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#else
#define NFFT_V __m128d
#define NFFT_VLEN 2
#define NFFT_VZERO() _mm_setzero_pd()
#define NFFT_VLOAD(p) _mm_loadu_pd(p)
#define NFFT_VSTORE(p, v) _mm_storeu_pd((p), (v))
#define NFFT_VADD(a, b) _mm_add_pd((a), (b))
#define NFFT_VFMA(a, b, c) _mm_add_pd(_mm_mul_pd((a), (b)), (c))
#endif

#elif NFFT_SIMD_VARIANT == NFFT_SIMD_AVX || NFFT_SIMD_VARIANT == NFFT_SIMD_AVX2

#include <immintrin.h>
#if NFFT_SIMD_VARIANT == NFFT_SIMD_AVX
#define NFFT_SIMD_SFX _avx
#define NFFT_SIMD_TGT __attribute__((target("avx")))
#else
#define NFFT_SIMD_SFX _avx2
#define NFFT_SIMD_TGT __attribute__((target("avx2,fma")))
#endif
#if defined(NFFT_SINGLE)
#define NFFT_V __m256
#define NFFT_VLEN 8
#define NFFT_VZERO() _mm256_setzero_ps()
#define NFFT_VLOAD(p) _mm256_loadu_ps(p)
#define NFFT_VSTORE(p, v) _mm256_storeu_ps((p), (v))
#define NFFT_VADD(a, b) _mm256_add_ps((a), (b))
#if NFFT_SIMD_VARIANT == NFFT_SIMD_AVX
#define NFFT_VFMA(a, b, c) _mm256_add_ps(_mm256_mul_ps((a), (b)), (c))
#else
#define NFFT_VFMA(a, b, c) _mm256_fmadd_ps((a), (b), (c))
#endif
#else
#define NFFT_V __m256d
#define NFFT_VLEN 4
#define NFFT_VZERO() _mm256_setzero_pd()
#define NFFT_VLOAD(p) _mm256_loadu_pd(p)
#define NFFT_VSTORE(p, v) _mm256_storeu_pd((p), (v))
#define NFFT_VADD(a, b) _mm256_add_pd((a), (b))
#if NFFT_SIMD_VARIANT == NFFT_SIMD_AVX
#define NFFT_VFMA(a, b, c) _mm256_add_pd(_mm256_mul_pd((a), (b)), (c))
#else
#define NFFT_VFMA(a, b, c) _mm256_fmadd_pd((a), (b), (c))
#endif
#endif

#elif NFFT_SIMD_VARIANT == NFFT_SIMD_NEON

#include <arm_neon.h>
#define NFFT_SIMD_SFX _neon
#define NFFT_SIMD_TGT
#if defined(NFFT_SINGLE)
#define NFFT_V float32x4_t
#define NFFT_VLEN 4
#define NFFT_VZERO() vdupq_n_f32(0.0f)
#define NFFT_VLOAD(p) vld1q_f32(p)
#define NFFT_VSTORE(p, v) vst1q_f32((p), (v))
#define NFFT_VADD(a, b) vaddq_f32((a), (b))
#define NFFT_VFMA(a, b, c) vfmaq_f32((c), (a), (b))
#else
#define NFFT_V float64x2_t
#define NFFT_VLEN 2
#define NFFT_VZERO() vdupq_n_f64(0.0)
#define NFFT_VLOAD(p) vld1q_f64(p)
#define NFFT_VSTORE(p, v) vst1q_f64((p), (v))
#define NFFT_VADD(a, b) vaddq_f64((a), (b))
#define NFFT_VFMA(a, b, c) vfmaq_f64((c), (a), (b))
#endif

#else
#error "unknown NFFT_SIMD_VARIANT"
#endif

#if NFFT_SIMD_VARIANT == NFFT_SIMD_SCALAR

/* The fallback: no intrinsics, no assumption about the layout of C, and the
 * same expression the transforms used before there was a vector path. */
static inline C NFFT_SIMD_ID(nfft_simd_cdot)(const R *psi, const C *g,
    const INT len)
{
  C s = K(0.0);
  INT l;

  for (l = 0; l < len; l++)
    s += psi[l] * g[l];

  return s;
}

static inline void NFFT_SIMD_ID(nfft_simd_caxpy)(C *g, const R *psi, const C f,
    const INT len)
{
  INT l;

  for (l = 0; l < len; l++)
    g[l] += psi[l] * f;
}

#else

/* Loads NFFT_VLEN/2 weights and duplicates each into an adjacent lane pair, so
 * that lane 2i and 2i+1 both hold psi[i] and line up with the real and the
 * imaginary part of g[i]. */
NFFT_SIMD_TGT static inline NFFT_V NFFT_SIMD_ID(nfft_simd_dup2)(const R *psi)
{
#if NFFT_SIMD_VARIANT == NFFT_SIMD_SSE2
#if defined(NFFT_SINGLE)
  const __m128 v = _mm_castsi128_ps(
      _mm_loadl_epi64((const __m128i *)(const void *)psi));
  return _mm_unpacklo_ps(v, v);
#else
  return _mm_load1_pd(psi);
#endif
#elif NFFT_SIMD_VARIANT == NFFT_SIMD_AVX || NFFT_SIMD_VARIANT == NFFT_SIMD_AVX2
#if defined(NFFT_SINGLE)
  const __m128 v = _mm_loadu_ps(psi);
  const __m256 t = _mm256_insertf128_ps(_mm256_castps128_ps256(v), v, 1);
  /* Per 128-bit lane: the low lane picks psi[0],psi[0],psi[1],psi[1] and the
   * high lane psi[2],psi[2],psi[3],psi[3], which one immediate cannot express
   * -- hence the index vector rather than _mm256_permute_ps. */
  return _mm256_permutevar_ps(t, _mm256_setr_epi32(0, 0, 1, 1, 2, 2, 3, 3));
#else
  const __m128d v = _mm_loadu_pd(psi);
  const __m256d t = _mm256_insertf128_pd(_mm256_castpd128_pd256(v), v, 1);
  return _mm256_permute_pd(t, 0xC);
#endif
#else /* NEON */
#if defined(NFFT_SINGLE)
  const float32x2_t v = vld1_f32(psi);
  return vcombine_f32(vdup_lane_f32(v, 0), vdup_lane_f32(v, 1));
#else
  return vld1q_dup_f64(psi);
#endif
#endif
}

/* Horizontal tail of a dot product: the accumulator holds partial real parts
 * in the even lanes and partial imaginary parts in the odd ones. Kept in
 * registers -- a run is only 2m+2 points long, so a round trip through memory
 * here would cost as much as the multiplications it closes. */
NFFT_SIMD_TGT static inline C NFFT_SIMD_ID(nfft_simd_reduce)(NFFT_V a)
{
#if NFFT_SIMD_VARIANT == NFFT_SIMD_AVX || NFFT_SIMD_VARIANT == NFFT_SIMD_AVX2
#if defined(NFFT_SINGLE)
  const __m128 h = _mm_add_ps(_mm256_castps256_ps128(a),
      _mm256_extractf128_ps(a, 1));
  const __m128 t = _mm_add_ps(h, _mm_movehl_ps(h, h));
  return _mm_cvtss_f32(t)
      + II * _mm_cvtss_f32(_mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 1, 1, 1)));
#else
  const __m128d h = _mm_add_pd(_mm256_castpd256_pd128(a),
      _mm256_extractf128_pd(a, 1));
  return _mm_cvtsd_f64(h) + II * _mm_cvtsd_f64(_mm_unpackhi_pd(h, h));
#endif
#elif NFFT_SIMD_VARIANT == NFFT_SIMD_SSE2
#if defined(NFFT_SINGLE)
  const __m128 t = _mm_add_ps(a, _mm_movehl_ps(a, a));
  return _mm_cvtss_f32(t)
      + II * _mm_cvtss_f32(_mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 1, 1, 1)));
#else
  return _mm_cvtsd_f64(a) + II * _mm_cvtsd_f64(_mm_unpackhi_pd(a, a));
#endif
#else /* NEON */
#if defined(NFFT_SINGLE)
  const float32x2_t t = vadd_f32(vget_low_f32(a), vget_high_f32(a));
  return vget_lane_f32(t, 0) + II * vget_lane_f32(t, 1);
#else
  return vgetq_lane_f64(a, 0) + II * vgetq_lane_f64(a, 1);
#endif
#endif
}

/** The scalar f spread over the vector as re, im, re, im, ... -- the shape a
 *  duplicated weight vector multiplies into. */
NFFT_SIMD_TGT static inline NFFT_V NFFT_SIMD_ID(nfft_simd_set2)(const R re,
    const R im)
{
#if NFFT_SIMD_VARIANT == NFFT_SIMD_AVX || NFFT_SIMD_VARIANT == NFFT_SIMD_AVX2
#if defined(NFFT_SINGLE)
  const __m128 v = _mm_set_ps(im, re, im, re);
  return _mm256_insertf128_ps(_mm256_castps128_ps256(v), v, 1);
#else
  const __m128d v = _mm_set_pd(im, re);
  return _mm256_insertf128_pd(_mm256_castpd128_pd256(v), v, 1);
#endif
#elif NFFT_SIMD_VARIANT == NFFT_SIMD_SSE2
#if defined(NFFT_SINGLE)
  return _mm_set_ps(im, re, im, re);
#else
  return _mm_set_pd(im, re);
#endif
#else /* NEON */
#if defined(NFFT_SINGLE)
  const float32x2_t v = vset_lane_f32(im, vdup_n_f32(re), 1);
  return vcombine_f32(v, v);
#else
  return vcombine_f64(vdup_n_f64(re), vdup_n_f64(im));
#endif
#endif
}

NFFT_SIMD_TGT static inline C NFFT_SIMD_ID(nfft_simd_cdot)(const R *psi,
    const C *g, const INT len)
{
  const R *gr = (const R *)g;
  const INT step = NFFT_VLEN / 2; /* complex points per vector */
  NFFT_V a0 = NFFT_VZERO(), a1 = NFFT_VZERO();
  C s;
  INT l = 0;

  /* Two accumulators: a run is short enough that the latency of one chain of
   * dependent multiply-adds, not the throughput, is what it costs. */
  for (; l + 2 * step <= len; l += 2 * step)
  {
    a0 = NFFT_VFMA(NFFT_SIMD_ID(nfft_simd_dup2)(psi + l),
        NFFT_VLOAD(gr + 2 * l), a0);
    a1 = NFFT_VFMA(NFFT_SIMD_ID(nfft_simd_dup2)(psi + l + step),
        NFFT_VLOAD(gr + 2 * (l + step)), a1);
  }

  if (l + step <= len)
  {
    a0 = NFFT_VFMA(NFFT_SIMD_ID(nfft_simd_dup2)(psi + l),
        NFFT_VLOAD(gr + 2 * l), a0);
    l += step;
  }

  s = NFFT_SIMD_ID(nfft_simd_reduce)(NFFT_VADD(a0, a1));

  for (; l < len; l++)
    s += psi[l] * g[l];

  return s;
}

NFFT_SIMD_TGT static inline void NFFT_SIMD_ID(nfft_simd_caxpy)(C *g,
    const R *psi, const C f, const INT len)
{
  R *gr = (R *)g;
  const INT step = NFFT_VLEN / 2;
  const NFFT_V fv = NFFT_SIMD_ID(nfft_simd_set2)(CREAL(f), CIMAG(f));
  INT l = 0;

  for (; l + step <= len; l += step)
    NFFT_VSTORE(gr + 2 * l, NFFT_VFMA(NFFT_SIMD_ID(nfft_simd_dup2)(psi + l),
        fv, NFFT_VLOAD(gr + 2 * l)));

  for (; l < len; l++)
    g[l] += psi[l] * f;
}

#endif /* NFFT_SIMD_VARIANT != NFFT_SIMD_SCALAR */

#undef NFFT_V
#undef NFFT_VLEN
#undef NFFT_VZERO
#undef NFFT_VLOAD
#undef NFFT_VSTORE
#undef NFFT_VADD
#undef NFFT_VFMA
