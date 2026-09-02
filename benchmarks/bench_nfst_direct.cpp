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
#include <time.h>

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

// Helper function to initialize random data
static void init_random_data(NFST(plan)* plan) {
    NFFT(vrand_real)(plan->x, plan->d * plan->M_total, K(0.0), K(0.5));
    NFFT(vrand_real)(plan->f_hat, plan->N_total, K(0.0), K(1.0));
    NFFT(vrand_real)(plan->f, plan->M_total, K(0.0), K(1.0));
}

// Benchmark for NFST direct transform (trafo_direct)
static void nfst_forward_direct_1d(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);

    NFST(plan) plan;
    NFST(init_1d)(&plan, N, M);
    init_random_data(&plan);

    for (auto _ : state) {
        NFST(trafo_direct)(&plan);
    }

    NFST(finalize)(&plan);
}

// Benchmark for NFST adjoint direct transform (adjoint_direct)
static void nfst_adjoint_direct_1d(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);

    NFST(plan) plan;
    NFST(init_1d)(&plan, N, M);
    init_random_data(&plan);

    for (auto _ : state) {
        NFST(adjoint_direct)(&plan);
    }

    NFST(finalize)(&plan);
}

// Benchmark for 2D NFST direct transform
static void nfst_forward_direct_2d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int M = state.range(2);

    NFST(plan) plan;
    NFST(init_2d)(&plan, N1, N2, M);
    init_random_data(&plan);

    for (auto _ : state) {
        NFST(trafo_direct)(&plan);
    }

    NFST(finalize)(&plan);
}

// Benchmark for 2D NFST adjoint direct transform
static void nfst_adjoint_direct_2d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int M = state.range(2);

    NFST(plan) plan;
    NFST(init_2d)(&plan, N1, N2, M);
    init_random_data(&plan);

    for (auto _ : state) {
        NFST(adjoint_direct)(&plan);
    }

    NFST(finalize)(&plan);
}

// Benchmark for 3D NFST direct transform
static void nfst_forward_direct_3d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int N3 = state.range(2);
    int M = state.range(3);

    NFST(plan) plan;
    NFST(init_3d)(&plan, N1, N2, N3, M);
    init_random_data(&plan);

    for (auto _ : state) {
        NFST(trafo_direct)(&plan);
    }

    NFST(finalize)(&plan);
}

// Benchmark for 3D NFST adjoint direct transform
static void nfst_adjoint_direct_3d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int N3 = state.range(2);
    int M = state.range(3);

    NFST(plan) plan;
    NFST(init_3d)(&plan, N1, N2, N3, M);
    init_random_data(&plan);

    for (auto _ : state) {
        NFST(adjoint_direct)(&plan);
    }

    NFST(finalize)(&plan);
}

// Register benchmarks for direct transforms
BENCH(nfst_forward_direct_1d, SUFFIX)
    ->Args({128, 400})
    ->Args({512, 1600})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nfst_adjoint_direct_1d, SUFFIX)
    ->Args({128, 400})
    ->Args({512, 1600})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nfst_forward_direct_2d, SUFFIX)
    ->Args({32, 32, 1000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nfst_adjoint_direct_2d, SUFFIX)
    ->Args({32, 32, 1000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nfst_forward_direct_3d, SUFFIX)
    ->Args({8, 8, 8, 500})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCH(nfst_adjoint_direct_3d, SUFFIX)
    ->Args({8, 8, 8, 500})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

// Main function.
BENCHMARK_MAIN();
