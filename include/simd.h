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

/*! \file simd.h
 *  \brief Which SIMD instruction sets this build can use, and which one the
 *  host actually runs.
 *
 *  Include this for the *detection* side only; simd_run.h holds the vector
 *  primitives themselves. It is included from infft.h, after the precision and
 *  the name mangling are fixed, and needs both.
 *
 *  The scheme has two independent halves:
 *
 *  - **Compile time.** NFFT_SIMD_HAVE_<isa> says a translation unit may hold
 *    code for that instruction set. Every such variant is compiled *in
 *    addition to* the scalar one, never instead of it, so a build always
 *    carries the plain C fallback.
 *  - **Run time.** Y(simd_isa)() asks the CPU what it supports and returns the
 *    widest variant that is both compiled in and available on this host. A
 *    binary built on an AVX2 machine therefore still runs on a host without
 *    AVX2 -- it takes the SSE2 or the scalar path instead.
 *
 *  Widening beyond the ISA the build targets rests on GCC's and clang's
 *  function target attribute, which lets one object file hold functions
 *  compiled for different instruction sets. Without it (other compilers, or
 *  --disable-simd) nothing but the scalar variant is built and Y(simd_isa)()
 *  answers NFFT_SIMD_SCALAR.
 */
#ifndef __NFFT_SIMD_H__
#define __NFFT_SIMD_H__

/* Instruction set identifiers. Ordered by increasing capability *within one
 * architecture*: a larger id on the same machine subsumes the smaller ones, so
 * the detection code can compare them. Ids from different architectures never
 * meet on one host. */
#define NFFT_SIMD_SCALAR 0
#define NFFT_SIMD_SSE2 1
#define NFFT_SIMD_AVX 2
#define NFFT_SIMD_AVX2 3
#define NFFT_SIMD_NEON 4

/* The vector primitives work on R as a plain array of reals and rely on the
 * complex type having the layout of two adjacent reals. That holds for float
 * and double; long double is padded (10 bytes of value in 12 or 16), so it
 * stays scalar, as does any build configured with --disable-simd. */
#if defined(NFFT_ENABLE_SIMD) && !defined(NFFT_LDOUBLE)

/* Compilers whose function target attribute we rely on. The version gate is
 * GCC's: 4.9 is the first release that accepts intrinsics for an instruction
 * set the command line did not enable. Compilers that merely define __GNUC__
 * for compatibility are excluded rather than assumed to follow. */
#if defined(__clang__)
#define NFFT_SIMD_TARGET_ATTR 1
#elif defined(__INTEL_COMPILER) || defined(__NVCOMPILER) || defined(__PGI) \
    || defined(__IBMC__) || defined(__SUNPRO_C)
/* No target attribute assumed. */
#elif defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 9))
#define NFFT_SIMD_TARGET_ATTR 1
#endif

#if defined(NFFT_SIMD_TARGET_ATTR) && (defined(__i386__) || defined(__x86_64__))
#define NFFT_SIMD_HAVE_SSE2 1
#define NFFT_SIMD_HAVE_AVX 1
#define NFFT_SIMD_HAVE_AVX2 1
#endif

/* Advanced SIMD is part of the AArch64 baseline, in both single and double
 * precision, so it needs no target attribute and no run-time check. 32-bit ARM
 * is left out: NEON is optional there and single precision only. */
#if defined(__aarch64__) || defined(_M_ARM64)
#define NFFT_SIMD_HAVE_NEON 1
#endif

#endif /* NFFT_ENABLE_SIMD && !NFFT_LDOUBLE */

#ifndef NFFT_SIMD_HAVE_SSE2
#define NFFT_SIMD_HAVE_SSE2 0
#endif
#ifndef NFFT_SIMD_HAVE_AVX
#define NFFT_SIMD_HAVE_AVX 0
#endif
#ifndef NFFT_SIMD_HAVE_AVX2
#define NFFT_SIMD_HAVE_AVX2 0
#endif
#ifndef NFFT_SIMD_HAVE_NEON
#define NFFT_SIMD_HAVE_NEON 0
#endif

/* The widest variant this build carries. Equal to NFFT_SIMD_SCALAR when the
 * build has no SIMD at all, which lets callers drop the dispatch entirely. */
#if NFFT_SIMD_HAVE_NEON
#define NFFT_SIMD_MAX NFFT_SIMD_NEON
#elif NFFT_SIMD_HAVE_AVX2
#define NFFT_SIMD_MAX NFFT_SIMD_AVX2
#elif NFFT_SIMD_HAVE_AVX
#define NFFT_SIMD_MAX NFFT_SIMD_AVX
#elif NFFT_SIMD_HAVE_SSE2
#define NFFT_SIMD_MAX NFFT_SIMD_SSE2
#else
#define NFFT_SIMD_MAX NFFT_SIMD_SCALAR
#endif

/** Instruction set the transforms use, resolved against the host on first call
 *  and cached. Returns NFFT_SIMD_SCALAR when nothing wider is both compiled in
 *  and supported here. The environment variable NFFT_SIMD overrides the
 *  detection: "scalar" (equivalently "0", "off", "none"), "sse2", "avx",
 *  "avx2" or "neon", each clamped down to what this build actually has. */
int Y(simd_isa)(void);

/** Whether \p isa is compiled into this build *and* supported by this host. */
int Y(simd_isa_available)(int isa);

/** Pins the instruction set Y(simd_isa)() reports, clamped to what is
 *  available. Plans pick the kernels up when they are initialised, so call
 *  this before X(init...). Meant for tests and for reproducing a run on a
 *  narrower host; there is no need to call it in normal use. */
void Y(simd_force_isa)(int isa);

/** Lower-case name of \p isa ("scalar", "sse2", "avx", "avx2", "neon"), or
 *  "unknown" for an id this build does not know. Never NULL. */
const char *Y(simd_isa_name)(int isa);

#endif /* __NFFT_SIMD_H__ */
