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
#include <benchmark/benchmark.h>
#include "config.h"

#include <complex.h>
#include <stdlib.h>

#include "nfft3.h"
#include "infft.h"

#include "util.h"

#ifdef _OPENMP
  #define SUFFIX "_omp"
#else
  #define SUFFIX ""
#endif

static void DoSetup(const benchmark::State& state) {
    nfft_bench_cap_threads();
    nfft_bench_align_allocations();
    #ifdef _OPENMP
    #ifdef HAVE_FFTW_THREADS
    FFTW(init_threads)();
    #endif  
    #endif
}

static void DoTeardown(const benchmark::State& state) {
    #ifdef _OPENMP
    #ifdef HAVE_FFTW_THREADS
    FFTW(cleanup_threads)();
    #endif
    #endif
}

/* Direct NNFFT on the geometries of bench_nfft_direct.cpp, NV = N_total. */
static void init_plan(nnfft_plan* plan, int d, int* N, int M) {
    int NV = 1;
    for (int t = 0; t < d; t++)
        NV *= N[t];
    nnfft_init(plan, d, NV, M, N);
    NFFT(vrand_shifted_unit_double)(plan->x, d * plan->M_total);
    NFFT(vrand_shifted_unit_double)(plan->v, d * plan->N_total);
    NFFT(vrand_unit_complex)(plan->f_hat, plan->N_total);
    NFFT(vrand_unit_complex)(plan->f, plan->M_total);
}

static void run(benchmark::State& state, int d, bool adjoint) {
    int N[3];
    for (int t = 0; t < d; t++)
        N[t] = (int)state.range(t);

    nnfft_plan plan;
    init_plan(&plan, d, N, (int)state.range(d));

    for (auto _ : state) {
        if (adjoint)
            nnfft_adjoint_direct(&plan);
        else
            nnfft_trafo_direct(&plan);
    }

    nnfft_finalize(&plan);
}

static void nnfft_forward_direct_1d(benchmark::State& state) { run(state, 1, false); }
static void nnfft_adjoint_direct_1d(benchmark::State& state) { run(state, 1, true); }
static void nnfft_forward_direct_2d(benchmark::State& state) { run(state, 2, false); }
static void nnfft_adjoint_direct_2d(benchmark::State& state) { run(state, 2, true); }
static void nnfft_forward_direct_3d(benchmark::State& state) { run(state, 3, false); }
static void nnfft_adjoint_direct_3d(benchmark::State& state) { run(state, 3, true); }

BENCH(nnfft_forward_direct_1d, SUFFIX)
    ->Args({128, 400})
    ->Args({512, 1600})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nnfft_adjoint_direct_1d, SUFFIX)
    ->Args({128, 400})
    ->Args({512, 1600})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nnfft_forward_direct_2d, SUFFIX)
    ->Args({32, 32, 1000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nnfft_adjoint_direct_2d, SUFFIX)
    ->Args({32, 32, 1000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nnfft_forward_direct_3d, SUFFIX)
    ->Args({8, 8, 8, 500})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nnfft_adjoint_direct_3d, SUFFIX)
    ->Args({8, 8, 8, 500})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

// Main function.
BENCHMARK_MAIN();
