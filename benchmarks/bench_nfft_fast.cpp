/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* Benchmarks for the fast NFFT with the standard PRE_PSI precomputation
 * scheme. Precomputation and transform are measured by separate benchmarks;
 * the transform benchmarks reuse the psi table left behind by the
 * precomputation benchmark for the same geometry. */

#include <benchmark/benchmark.h>
#include "config.h"

#include <complex.h>
#include <string.h>
#include <stdlib.h>

#include "nfft3.h"
#include "infft.h"

#include "util.h"

#ifdef _OPENMP
  #define SUFFIX "_omp"
#else
  #define SUFFIX ""
#endif

#define BENCH_MAX_D 4
#define BENCH_SEED 4711L

/* The window cutoff nfft_init_* picks for the configured window and precision. */
#define DEFAULT_M WINDOW_HELP_ESTIMATE_m
/* Sweep cutoffs, relative to the default so no sweep case collides with the
 * default-m case in any window/precision configuration. */
#define SMALL_M (DEFAULT_M / 2)
#define LARGE_M (DEFAULT_M * 2)

struct Geometry {
    int d;
    int N[BENCH_MAX_D];
    int M;
    int m;
};

static bool same_geometry(const Geometry& a, const Geometry& b) {
    if (a.d != b.d || a.M != b.M || a.m != b.m)
        return false;
    for (int t = 0; t < a.d; t++)
        if (a.N[t] != b.N[t])
            return false;
    return true;
}

/* Plan cache. At most one plan is alive at any time: requesting a
 * different geometry finalizes the previous one. Benchmarks are registered
 * grouped by geometry, so precomputation and both transforms for one geometry
 * share a single init. */
struct PlanSlot {
    bool valid = false;
    bool psi_ready = false;
    Geometry geom = {};
    NFFT(plan) plan = {};
};

static PlanSlot slot;
static bool fftw_threads_started = false;

static void release_plan(void) {
    if (slot.valid) {
        NFFT(finalize)(&slot.plan);
        slot.valid = false;
        slot.psi_ready = false;
    }
}

/* fftw_cleanup_threads() invalidates every plan created so far, so it must not
 * run between benchmarks while a cached plan is still alive. Both it and the
 * final release happen once, at exit. */
struct SlotCleanup {
    ~SlotCleanup() {
        release_plan();
        #ifdef _OPENMP
        #ifdef HAVE_FFTW_THREADS
        if (fftw_threads_started)
            FFTW(cleanup_threads)();
        #endif
        #endif
    }
};
static SlotCleanup slot_cleanup;

static void DoSetup(const benchmark::State& state) {
    #ifdef _OPENMP
    #ifdef HAVE_FFTW_THREADS
    if (!fftw_threads_started) {
        FFTW(init_threads)();
        fftw_threads_started = true;
    }
    #endif
    #endif
}

/* Never timed. With fresh, the plan is rebuilt even when the cache already
 * holds this geometry, so its psi table has never been filled. */
static NFFT(plan)* acquire_plan(const Geometry& geom, bool fresh = false) {
    if (!fresh && slot.valid && same_geometry(slot.geom, geom))
        return &slot.plan;

    release_plan();

    NFFT(srand48)(BENCH_SEED);

    int N[BENCH_MAX_D];
    int n[BENCH_MAX_D];
    for (int t = 0; t < geom.d; t++) {
        N[t] = geom.N[t];
        n[t] = 2 * (int)NFFT(next_power_of_2)(geom.N[t]);  /* sigma = 2 */
    }

    NFFT(init_guru)(&slot.plan, geom.d, N, geom.M, n, geom.m,
        PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
        | FFTW_INIT | FFT_OUT_OF_PLACE,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    NFFT(vrand_shifted_unit_double)(slot.plan.x, slot.plan.d * slot.plan.M_total);

    slot.valid = true;
    slot.psi_ready = false;
    slot.geom = geom;
    return &slot.plan;
}

/* Leaves x untouched: psi depends only on x, so a cached psi table stays valid. */
static void reseed_data(NFFT(plan)* plan) {
    NFFT(srand48)(BENCH_SEED);
    NFFT(vrand_unit_complex)(plan->f_hat, plan->N_total);
    NFFT(vrand_unit_complex)(plan->f, plan->M_total);
}

static void ensure_psi(NFFT(plan)* plan) {
    if (!slot.psi_ready) {
        NFFT(precompute_one_psi)(plan);
        slot.psi_ready = true;
    }
}

static Geometry geometry_from(const benchmark::State& state, int d) {
    Geometry geom = {};
    geom.d = d;
    for (int t = 0; t < d; t++)
        geom.N[t] = (int)state.range(t);
    geom.M = (int)state.range(d);
    geom.m = (int)state.range(d + 1);
    return geom;
}

/* Zero the psi table to incur initial page faults outside the timed loop.
 * len must stay equal to the PRE_PSI allocation in init_help (nfft.c),
 * M_total * d * (2m+2) reals. */
static void prefault_psi(NFFT(plan)* plan) {
    const size_t len = (size_t)plan->M_total * (size_t)plan->d
        * (size_t)(2 * plan->m + 2);
    memset(plan->psi, 0, len * sizeof(R));
}

/* Times one nfft_precompute_one_psi on a plan whose psi table has never been
 * filled. Nothing in the API promises that a second call on the same plan
 * redoes the work, so every measured call gets its own plan; Iterations(1) at
 * the registration keeps that true, and repetitions, which re-enter this
 * function, supply the statistics. Building the plan here rather than inside
 * the loop keeps init and allocation out of the measurement in both walltime
 * and simulation mode. */
static void run_precompute(benchmark::State& state, int d) {
    NFFT(plan)* plan = acquire_plan(geometry_from(state, d), /*fresh=*/true);
    prefault_psi(plan);

    for (auto _ : state) {
        NFFT(precompute_one_psi)(plan);
        benchmark::ClobberMemory();
    }

    slot.psi_ready = true;
}

static void run_trafo(benchmark::State& state, int d) {
    NFFT(plan)* plan = acquire_plan(geometry_from(state, d));
    ensure_psi(plan);
    reseed_data(plan);

    for (auto _ : state) {
        NFFT(trafo)(plan);
        benchmark::ClobberMemory();
    }
}

static void run_adjoint(benchmark::State& state, int d) {
    NFFT(plan)* plan = acquire_plan(geometry_from(state, d));
    ensure_psi(plan);
    reseed_data(plan);

    for (auto _ : state) {
        NFFT(adjoint)(plan);
        benchmark::ClobberMemory();
    }
}

#define DEFINE_DIM(tag, dim) \
    static void nfft_fast_precompute_psi_##tag(benchmark::State& state) { \
        run_precompute(state, dim); \
    } \
    static void nfft_fast_trafo_##tag(benchmark::State& state) { \
        run_trafo(state, dim); \
    } \
    static void nfft_fast_adjoint_##tag(benchmark::State& state) { \
        run_adjoint(state, dim); \
    }

DEFINE_DIM(1d, 1)
DEFINE_DIM(2d, 2)
DEFINE_DIM(3d, 3)
DEFINE_DIM(4d, 4)

/* Args are N[0..d-1], M, m. */
#define REGISTER_CASE(tag, ...) \
    BENCH(nfft_fast_precompute_psi_##tag, SUFFIX) \
        ->Args({__VA_ARGS__})->Iterations(1)->Setup(DoSetup); \
    BENCH(nfft_fast_trafo_##tag, SUFFIX) \
        ->Args({__VA_ARGS__})->Setup(DoSetup); \
    BENCH(nfft_fast_adjoint_##tag, SUFFIX) \
        ->Args({__VA_ARGS__})->Setup(DoSetup);

/* 1d size sweep. */
REGISTER_CASE(1d, 1024, 1024, DEFAULT_M)
REGISTER_CASE(1d, 8192, 8192, DEFAULT_M)
REGISTER_CASE(1d, 65536, 65536, DEFAULT_M)

/* 1d, off the M = N_total diagonal: separates FFT-phase from B-phase cost. */
REGISTER_CASE(1d, 8192, 1024, DEFAULT_M)
REGISTER_CASE(1d, 8192, 65536, DEFAULT_M)

/* 1d cutoff sweep: B-phase cost scales as (2m+2)^d. */
REGISTER_CASE(1d, 8192, 8192, SMALL_M)
REGISTER_CASE(1d, 8192, 8192, LARGE_M)

/* 2d size sweep. */
REGISTER_CASE(2d, 32, 32, 1024, DEFAULT_M)
REGISTER_CASE(2d, 128, 128, 16384, DEFAULT_M)
REGISTER_CASE(2d, 256, 256, 65536, DEFAULT_M)

/* 3d size sweep. */
REGISTER_CASE(3d, 8, 8, 8, 512, DEFAULT_M)
REGISTER_CASE(3d, 16, 16, 16, 2048, DEFAULT_M)
REGISTER_CASE(3d, 32, 32, 32, 4096, DEFAULT_M)

/* d = 4 uses the generic path instead of the specialized 1d/2d/3d kernels. */
REGISTER_CASE(4d, 8, 8, 8, 8, 256, DEFAULT_M)

BENCHMARK_MAIN();
