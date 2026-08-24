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

/* Planner-native (plan_ng) counterpart of bench_nfft_direct.cpp: the same six
 * transforms at the same sizes, so legacy and planner direct NDFT are directly
 * comparable in CodSpeed. NFFT_NO_FAST_NATIVE leaves only the direct solvers,
 * and NFFT_ESTIMATE skips the measured race, which would otherwise both cost
 * setup time and write into f_hat/f. */

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

// plan_ng allocates none of the user arrays, so the benchmark owns x, f_hat
// and f. Oversampling and window cut-off repeat the NFFT(init) defaults that
// bench_nfft_direct.cpp gets implicitly through init_1d/2d/3d.
struct ng_fixture {
    R *x;
    FC *f_hat;
    FC *f;
    NFFT(plan_ng) *p;
};

static void ng_setup(ng_fixture* fx, int d, const NFFT_INT* N, NFFT_INT M) {
    NFFT_INT n[3];
    NFFT_INT N_total = 1;

    for (int t = 0; t < d; t++) {
        N_total *= N[t];
        n[t] = 2 * NFFT(next_power_of_2)(N[t]);
    }

    fx->x = (R*) NFFT(malloc)((size_t)(d * M) * sizeof(R));
    fx->f_hat = (FC*) NFFT(malloc)((size_t)N_total * sizeof(FC));
    fx->f = (FC*) NFFT(malloc)((size_t)M * sizeof(FC));

    // x is copied into the plan at construction, so fill it first.
    NFFT(vrand_shifted_unit_double)(fx->x, d * M);

    fx->p = NFFT(plan_ng_guru)(d, N, NULL, n, M, WINDOW_HELP_ESTIMATE_m,
                               NFFT(get_window_id)(), fx->x, fx->f_hat, fx->f,
                               0u, NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);

    NFFT(vrand_unit_complex)(fx->f_hat, N_total);
    NFFT(vrand_unit_complex)(fx->f, M);

    NFFT(precompute)(fx->p);
}

static void ng_teardown(ng_fixture* fx) {
    NFFT(plan_ng_destroy)(fx->p);
    NFFT(free)(fx->f);
    NFFT(free)(fx->f_hat);
    NFFT(free)(fx->x);
}

static void nfft_ng_forward_direct_1d(benchmark::State& state) {
    NFFT_INT N = (NFFT_INT) state.range(0);
    NFFT_INT M = (NFFT_INT) state.range(1);

    ng_fixture fx;
    ng_setup(&fx, 1, &N, M);

    for (auto _ : state) {
        NFFT(execute)(fx.p);
    }

    ng_teardown(&fx);
    state.SetComplexityN(N * M);
}

static void nfft_ng_adjoint_direct_1d(benchmark::State& state) {
    NFFT_INT N = (NFFT_INT) state.range(0);
    NFFT_INT M = (NFFT_INT) state.range(1);

    ng_fixture fx;
    ng_setup(&fx, 1, &N, M);

    for (auto _ : state) {
        NFFT(execute_adjoint)(fx.p);
    }

    ng_teardown(&fx);
    state.SetComplexityN(N * M);
}

static void nfft_ng_forward_direct_2d(benchmark::State& state) {
    NFFT_INT N[2] = {(NFFT_INT) state.range(0), (NFFT_INT) state.range(1)};
    NFFT_INT M = (NFFT_INT) state.range(2);

    ng_fixture fx;
    ng_setup(&fx, 2, N, M);

    for (auto _ : state) {
        NFFT(execute)(fx.p);
    }

    ng_teardown(&fx);
    state.SetComplexityN(N[0] * N[1] * M);
}

static void nfft_ng_adjoint_direct_2d(benchmark::State& state) {
    NFFT_INT N[2] = {(NFFT_INT) state.range(0), (NFFT_INT) state.range(1)};
    NFFT_INT M = (NFFT_INT) state.range(2);

    ng_fixture fx;
    ng_setup(&fx, 2, N, M);

    for (auto _ : state) {
        NFFT(execute_adjoint)(fx.p);
    }

    ng_teardown(&fx);
    state.SetComplexityN(N[0] * N[1] * M);
}

static void nfft_ng_forward_direct_3d(benchmark::State& state) {
    NFFT_INT N[3] = {(NFFT_INT) state.range(0), (NFFT_INT) state.range(1),
                     (NFFT_INT) state.range(2)};
    NFFT_INT M = (NFFT_INT) state.range(3);

    ng_fixture fx;
    ng_setup(&fx, 3, N, M);

    for (auto _ : state) {
        NFFT(execute)(fx.p);
    }

    ng_teardown(&fx);
    state.SetComplexityN(N[0] * N[1] * N[2] * M);
}

static void nfft_ng_adjoint_direct_3d(benchmark::State& state) {
    NFFT_INT N[3] = {(NFFT_INT) state.range(0), (NFFT_INT) state.range(1),
                     (NFFT_INT) state.range(2)};
    NFFT_INT M = (NFFT_INT) state.range(3);

    ng_fixture fx;
    ng_setup(&fx, 3, N, M);

    for (auto _ : state) {
        NFFT(execute_adjoint)(fx.p);
    }

    ng_teardown(&fx);
    state.SetComplexityN(N[0] * N[1] * N[2] * M);
}

BENCH(nfft_ng_forward_direct_1d, SUFFIX)
    ->Args({32, 100})
    ->Args({64, 200})
    ->Args({128, 400})
    ->Args({256, 800})
    ->Args({512, 1600})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown)
    ->Complexity();

BENCH(nfft_ng_adjoint_direct_1d, SUFFIX)
    ->Args({32, 100})
    ->Args({64, 200})
    ->Args({128, 400})
    ->Args({256, 800})
    ->Args({512, 1600})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown)
    ->Complexity();

BENCH(nfft_ng_forward_direct_2d, SUFFIX)
    ->Args({16, 16, 500})
    ->Args({32, 32, 1000})
    ->Args({64, 64, 2000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown)
    ->Complexity();

BENCH(nfft_ng_adjoint_direct_2d, SUFFIX)
    ->Args({16, 16, 500})
    ->Args({32, 32, 1000})
    ->Args({64, 64, 2000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown)
    ->Complexity();

BENCH(nfft_ng_forward_direct_3d, SUFFIX)
    ->Args({4, 4, 4, 250})
    ->Args({8, 8, 8, 500})
    ->Args({16, 16, 16, 1000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown)
    ->Complexity();

BENCH(nfft_ng_adjoint_direct_3d, SUFFIX)
    ->Args({4, 4, 4, 250})
    ->Args({8, 8, 8, 500})
    ->Args({16, 16, 16, 1000})
    ->Setup(DoSetup)
    ->Teardown(DoTeardown)
    ->Complexity();

BENCHMARK_MAIN();
