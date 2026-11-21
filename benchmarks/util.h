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

#ifdef __linux__
#include <sys/mman.h>
#endif

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

// Gives every plan array its own 2 MiB huge pages.
//
// The cache sets an array maps to depend on its physical address. With 4 KiB
// pages that address is random per process, and a plan rebuilt in the same
// process gets the same pages back, so all rounds of one run share a single
// placement. Runs then differ by up to 18% in discrete steps while the spread
// within a run stays near 0.1%. Inside a huge page the physical offset equals
// the virtual one, so the placement is the same in every run.
//
// Costs up to 2 MiB per array and fewer TLB misses than a caller with 4 KiB
// pages gets. Without transparent huge pages the alignment alone still fixes
// the virtual placement.
#define NFFT_BENCH_ALIGNMENT ((size_t)2 * 1024 * 1024)

static inline void *nfft_bench_aligned_malloc(size_t n)
{
    const size_t a = NFFT_BENCH_ALIGNMENT;
    void *p = NULL;

    n = n == 0 ? a : (n + (a - 1)) & ~(a - 1);

#ifdef _WIN32
    p = _aligned_malloc(n, a);
#else
    if (posix_memalign(&p, a, n) != 0)
        p = NULL;
#endif
    if (p == NULL)
        std::abort();
#ifdef __linux__
    madvise(p, n, MADV_HUGEPAGE);
#endif

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
