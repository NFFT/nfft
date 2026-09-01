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

#endif // NFFT_BENCHMARKS_UTIL_H
