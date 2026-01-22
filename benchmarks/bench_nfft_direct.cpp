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

// Helper function to initialize random data
static void NFFT(init_random_data)(NFFT(plan)* plan) {
    NFFT(vrand_shifted_unit_double)(plan->x, plan->d * plan->M_total);
    NFFT(vrand_unit_complex)(plan->f_hat, plan->N_total);
    NFFT(vrand_unit_complex)(plan->f, plan->M_total);
}

// Benchmark for NFFT direct transform (trafo_direct)
static void nfft_forward_direct_1d(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);
    
    NFFT(plan) plan;
    NFFT(init_1d)(&plan, N, M);
    NFFT(init_random_data)(&plan);
    
    for (auto _ : state) {
        NFFT(trafo_direct)(&plan);
    }
    
    NFFT(finalize)(&plan);
    state.SetComplexityN(N * M);
}

// Benchmark for NFFT adjoint direct transform (adjoint_direct)
static void nfft_adjoint_direct_1d(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);
    
    NFFT(plan) plan;
    NFFT(init_1d)(&plan, N, M);
    NFFT(init_random_data)(&plan);
    
    for (auto _ : state) {
        NFFT(adjoint_direct)(&plan);
    }
    
    NFFT(finalize)(&plan);
    state.SetComplexityN(N * M);
}

// Benchmark for 2D NFFT direct transform
static void nfft_forward_direct_2d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int M = state.range(2);
    
    NFFT(plan) plan;
    NFFT(init_2d)(&plan, N1, N2, M);
    NFFT(init_random_data)(&plan);
    
    for (auto _ : state) {
        NFFT(trafo_direct)(&plan);
    }
    
    NFFT(finalize)(&plan);
    state.SetComplexityN(N1 * N2 * M);
}

// Benchmark for 2D NFFT adjoint direct transform
static void nfft_adjoint_direct_2d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int M = state.range(2);
    
    NFFT(plan) plan;
    NFFT(init_2d)(&plan, N1, N2, M);
    NFFT(init_random_data)(&plan);
    
    for (auto _ : state) {
        NFFT(adjoint_direct)(&plan);
    }
    
    NFFT(finalize)(&plan);
    state.SetComplexityN(N1 * N2 * M);
}

// Benchmark for 3D NFFT direct transform
static void nfft_forward_direct_3d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int N3 = state.range(2);
    int M = state.range(3);
    
    NFFT(plan) plan;
    NFFT(init_3d)(&plan, N1, N2, N3, M);
    NFFT(init_random_data)(&plan);
    
    for (auto _ : state) {
        NFFT(trafo_direct)(&plan);
    }
    
    NFFT(finalize)(&plan);
    state.SetComplexityN(N1 * N2 * N3 * M);
}

// Benchmark for 3D NFFT adjoint direct transform
static void nfft_adjoint_direct_3d(benchmark::State& state) {
    int N1 = state.range(0);
    int N2 = state.range(1);
    int N3 = state.range(2);
    int M = state.range(3);
    
    NFFT(plan) plan;
    NFFT(init_3d)(&plan, N1, N2, N3, M);
    NFFT(init_random_data)(&plan);
    
    for (auto _ : state) {
        NFFT(adjoint_direct)(&plan);
    }
    
    NFFT(finalize)(&plan);
    state.SetComplexityN(N1 * N2 * N3 * M);
}

// Register benchmarks for direct transforms
BENCH(nfft_forward_direct_1d)
    ->Args({32, 100})
    ->Args({64, 200})
    ->Args({128, 400})
    ->Args({256, 800})
    //->Args({512, 1600})
    ->Complexity();

BENCH(nfft_adjoint_direct_1d)
    ->Args({32, 100})
    ->Args({64, 200})
    ->Args({128, 400})
    ->Args({256, 800})
    //->Args({512, 1600})
    ->Complexity();

BENCH(nfft_forward_direct_2d)
    ->Args({16, 16, 500})
    ->Args({32, 32, 1000})
    //->Args({64, 64, 2000})
    ->Complexity();

BENCH(nfft_adjoint_direct_2d)
    ->Args({16, 16, 500})
    ->Args({32, 32, 1000})
    //->Args({64, 64, 2000})
    ->Complexity();

BENCH(nfft_forward_direct_3d)
    ->Args({4, 4, 4, 250})
    ->Args({8, 8, 8, 500})
    //->Args({16, 16, 16, 1000})
    ->Complexity();

BENCH(nfft_adjoint_direct_3d)
    ->Args({4, 4, 4, 250})
    ->Args({8, 8, 8, 500})
    //->Args({16, 16, 16, 1000})
    ->Complexity();

// Main function.
BENCHMARK_MAIN();