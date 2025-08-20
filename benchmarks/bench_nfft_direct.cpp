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

// Helper function to initialize random data
static void init_random_data(nfft_plan* plan) {
    NFFT(vrand_shifted_unit_double)(plan->x, plan->d * plan->M_total);
    NFFT(vrand_unit_complex)(plan->f_hat, plan->N_total);
    NFFT(vrand_unit_complex)(plan->f, plan->M_total);
}

// Benchmark for NFFT direct transform (trafo_direct)
static void BM_NFFT_TrafoDirect(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);
    
    nfft_plan plan;
    nfft_init_1d(&plan, N, M);
    init_random_data(&plan);
    
    for (auto _ : state) {
        nfft_trafo_direct(&plan);
    }
    
    nfft_finalize(&plan);
    state.SetComplexityN(N * M);
}

// Benchmark for NFFT adjoint direct transform (adjoint_direct)
static void BM_NFFT_AdjointDirect(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);
    
    nfft_plan plan;
    nfft_init_1d(&plan, N, M);
    init_random_data(&plan);
    
    for (auto _ : state) {
        nfft_adjoint_direct(&plan);
    }
    
    nfft_finalize(&plan);
    state.SetComplexityN(N * M);
}

// Benchmark for 2D NFFT direct transform
static void BM_NFFT_2D_TrafoDirect(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int M = state.range(2);
    
    nfft_plan plan;
    nfft_init_2d(&plan, N1, N2, M);
    init_random_data(&plan);
    
    for (auto _ : state) {
        nfft_trafo_direct(&plan);
    }
    
    nfft_finalize(&plan);
    state.SetComplexityN(N1 * N2 * M);
}

// Benchmark for 2D NFFT adjoint direct transform
static void BM_NFFT_2D_AdjointDirect(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int M = state.range(2);
    
    nfft_plan plan;
    nfft_init_2d(&plan, N1, N2, M);
    init_random_data(&plan);
    
    for (auto _ : state) {
        nfft_adjoint_direct(&plan);
    }
    
    nfft_finalize(&plan);
    state.SetComplexityN(N1 * N2 * M);
}

// Benchmark for 3D NFFT direct transform
static void BM_NFFT_3D_TrafoDirect(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int N3 = state.range(2);
    int M = state.range(3);
    
    nfft_plan plan;
    nfft_init_3d(&plan, N1, N2, N3, M);
    init_random_data(&plan);
    
    for (auto _ : state) {
        nfft_trafo_direct(&plan);
    }
    
    nfft_finalize(&plan);
    state.SetComplexityN(N1 * N2 * N3 * M);
}

// Benchmark for 3D NFFT adjoint direct transform
static void BM_NFFT_3D_AdjointDirect(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int N3 = state.range(2);
    int M = state.range(3);
    
    nfft_plan plan;
    nfft_init_3d(&plan, N1, N2, N3, M);
    init_random_data(&plan);
    
    for (auto _ : state) {
        nfft_adjoint_direct(&plan);
    }
    
    nfft_finalize(&plan);
    state.SetComplexityN(N1 * N2 * N3 * M);
}

// Register benchmarks for direct transforms
BENCHMARK(BM_NFFT_TrafoDirect)
    ->Args({32, 100})
    ->Args({64, 200})
    ->Args({128, 400})
    ->Args({256, 800})
    ->Args({512, 1600})
    ->Complexity();

BENCHMARK(BM_NFFT_AdjointDirect)
    ->Args({32, 100})
    ->Args({64, 200})
    ->Args({128, 400})
    ->Args({256, 800})
    ->Args({512, 1600})
    ->Complexity();

BENCHMARK(BM_NFFT_2D_TrafoDirect)
    ->Args({16, 16, 500})
    ->Args({32, 32, 1000})
    ->Args({64, 64, 2000})
    ->Complexity();

BENCHMARK(BM_NFFT_2D_AdjointDirect)
    ->Args({16, 16, 500})
    ->Args({32, 32, 1000})
    ->Args({64, 64, 2000})
    ->Complexity();

BENCHMARK(BM_NFFT_3D_TrafoDirect)
    ->Args({4, 4, 4, 250})
    ->Args({8, 8, 8, 500})
    ->Args({16, 16, 16, 1000})
    ->Complexity();

BENCHMARK(BM_NFFT_3D_AdjointDirect)
    ->Args({4, 4, 4, 250})
    ->Args({8, 8, 8, 500})
    ->Args({16, 16, 16, 1000})
    ->Complexity();

// Main function.
BENCHMARK_MAIN();