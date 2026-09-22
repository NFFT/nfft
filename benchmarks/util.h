/*
 * Copyright (c) 2025 Jens Keiner, Stefan Kunis, Daniel Potts
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

#ifndef NFFT_BENCHMARKS_UTIL_H
#define NFFT_BENCHMARKS_UTIL_H

#include "config.h"

#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

// Macro to register benchmark with optional prefix.
#define BENCH(function, suffix) BENCHMARK(function)->Name(BENCHMARKS_PREFIX #function suffix)

// Caps the OpenMP team unless the caller asked for a size. A small fixed team 
// is steadier and independent of the host.
#define NFFT_BENCH_MAX_THREADS 2

static inline void nfft_bench_cap_threads(void)
{
#ifdef _OPENMP
    if (std::getenv("OMP_NUM_THREADS") != nullptr)
        return;

    {
        const int max = omp_get_max_threads();
        omp_set_num_threads(max < NFFT_BENCH_MAX_THREADS ? max
                                                         : NFFT_BENCH_MAX_THREADS);
    }
#endif
}

// Aligns every plan allocation to (hopefully) a cache line and rounds its size up to one.
//
// nfft_malloc is fftw_malloc, which promises 16 bytes: enough for the SIMD
// codelets, but it leaves the plan's small heap arrays wherever the allocator
// puts them. Change the size of any early allocation and x, f and c_phi_inv all
// shift, which moves a cache-resident benchmark by several percent for reasons
// that have nothing to do with the code under test. Rounding the size up is
// what absorbs the change; aligning alone does not.
//
// Costs a little memory and measures a slightly better aligned configuration
// than a caller who uses nfft_malloc directly gets.
//
// A cache line is not enough for the larger arrays. An array of a megabyte or
// more comes from mmap, whose placement moves with the state of the allocator,
// so the cache sets the array lands on change when an unrelated allocation
// changes size. That decides how much of the array stays resident and moves a
// benchmark at the cache boundary by more than ten percent. Each array is
// therefore aligned to its own size, rounded down to a power of two and capped
// at the large page size, which fixes where in the cache it sits.
#define NFFT_BENCH_ALIGNMENT 64
#define NFFT_BENCH_MAX_ALIGNMENT ((size_t)2 * 1024 * 1024)

static inline void *nfft_bench_aligned_malloc(size_t n)
{
    size_t a = NFFT_BENCH_ALIGNMENT;
    void *p = NULL;

    if (n == 0)
        n = 1;
    while (a < NFFT_BENCH_MAX_ALIGNMENT && a * 2 <= n)
        a *= 2;
    n = (n + (a - 1)) & ~(a - 1);

#ifdef _WIN32
    p = _aligned_malloc(n, a);
#else
    if (posix_memalign(&p, a, n) != 0)
        p = NULL;
#endif
    if (p == NULL)
        std::abort();

    return p;
}

static inline void nfft_bench_aligned_free(void *p)
{
#ifdef _WIN32
    _aligned_free(p);
#else
    std::free(p);
#endif
}

// Both hooks go in together: nfft_free falls back to fftw_free when only one is
// set, which would not match the allocation.
static inline void nfft_bench_align_allocations(void)
{
    NFFT(malloc_hook) = nfft_bench_aligned_malloc;
    NFFT(free_hook) = nfft_bench_aligned_free;
}

#endif // NFFT_BENCHMARKS_UTIL_H
