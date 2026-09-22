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
 * scheme. Precomputation and transform are measured by separate benchmarks.
 *
 * Every round builds its own plan behind a spacer of varying size, so buffer
 * placement varies across the rounds of one run and lands inside the reported
 * spread. Placement moves these transforms by up to 17%, far more than the
 * noise within one placement. */

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

/* At most one plan is alive at any time: each acquire finalizes the previous. */

/* Which benchmark last seeded f_hat and f. Neither transform modifies its own
 * input so repeated entries into one benchmark need no reseeding. Moving between the
 * two does. */
enum DataRole { DATA_NONE, DATA_TRAFO, DATA_ADJOINT };

struct PlanSlot {
    bool valid = false;
    bool psi_ready = false;
    DataRole data_role = DATA_NONE;
    Geometry geom = {};
    NFFT(plan) plan = {};
};

static PlanSlot slot;
static bool fftw_threads_started = false;

/* Rebuilding a plan of the same geometry hands back the same addresses, so the
 * spacer must stay alive across the round. The size sequence is fixed, so both
 * sides of a comparison walk the same placements. */
static void *layout_spacer = NULL;
static unsigned long spacer_state = 1UL;

static void shift_layout(void) {
    free(layout_spacer);
    spacer_state = spacer_state * 1103515245UL + 12345UL;
    layout_spacer = malloc(64 + (size_t)((spacer_state >> 13) & 0x3fff) * 64);
}

static void release_plan(void) {
    if (slot.valid) {
        NFFT(finalize)(&slot.plan);
        slot.valid = false;
        slot.psi_ready = false;
        slot.data_role = DATA_NONE;
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

/* Flag to indicate if OpenMP thread team has been created. It's a one-off cost 
 * of the first parallel region in the process so it should be done once outside
 * a timing loop. */
static bool omp_team_started = false;

static void warm_omp_team(void) {
    #ifdef _OPENMP
    int sink = 0;
    #pragma omp parallel reduction(+:sink)
    {
        sink += 1;
    }
    benchmark::DoNotOptimize(sink);
    #endif
}

static void DoSetup(const benchmark::State& state) {
    nfft_bench_cap_threads();
    nfft_bench_align_allocations();
    #ifdef _OPENMP
    #ifdef HAVE_FFTW_THREADS
    if (!fftw_threads_started) {
        FFTW(init_threads)();
        fftw_threads_started = true;
    }
    #endif
    if (!omp_team_started) {
        warm_omp_team();
        omp_team_started = true;
    }
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
    slot.data_role = DATA_NONE;
    slot.geom = geom;
    return &slot.plan;
}

/* Leaves x untouched: psi depends only on x, so a cached psi table stays valid. */
static void reseed_data(NFFT(plan)* plan) {
    NFFT(srand48)(BENCH_SEED);
    NFFT(vrand_unit_complex)(plan->f_hat, plan->N_total);
    NFFT(vrand_unit_complex)(plan->f, plan->M_total);
}

/* Seeds only when the data on the plan belongs to another benchmark. */
static void ensure_data(NFFT(plan)* plan, DataRole role) {
    if (slot.data_role == role)
        return;
    reseed_data(plan);
    slot.data_role = role;
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

/* Manual warm-up, untimed: Iterations() excludes MinTime(), which
 * MinWarmUpTime() requires. Two faults in the fresh buffers. */
#define BENCH_WARMUP_ITERS 2

/* Times nfft_precompute_one_psi on a plan built here, so init and allocation
 * stay out of the measurement. */
static void run_precompute(benchmark::State& state, int d) {
    shift_layout();
    NFFT(plan)* plan = acquire_plan(geometry_from(state, d), /*fresh=*/true);
    prefault_psi(plan);

    for (auto _ : state) {
        NFFT(precompute_one_psi)(plan);
        benchmark::ClobberMemory();
    }

    slot.psi_ready = true;
}

static void run_trafo(benchmark::State& state, int d) {
    shift_layout();
    NFFT(plan)* plan = acquire_plan(geometry_from(state, d), /*fresh=*/true);
    ensure_psi(plan);
    ensure_data(plan, DATA_TRAFO);

    for (int w = 0; w < BENCH_WARMUP_ITERS; w++)
        NFFT(trafo)(plan);

    for (auto _ : state) {
        NFFT(trafo)(plan);
        benchmark::ClobberMemory();
    }
}

static void run_adjoint(benchmark::State& state, int d) {
    shift_layout();
    NFFT(plan)* plan = acquire_plan(geometry_from(state, d), /*fresh=*/true);
    ensure_psi(plan);
    ensure_data(plan, DATA_ADJOINT);

    for (int w = 0; w < BENCH_WARMUP_ITERS; w++)
        NFFT(adjoint)(plan);

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

#define BENCH_BUDGET(name, iters) \
    BENCH(name, SUFFIX)->Iterations(iters)->Setup(DoSetup)

/* Iterations per round, stated per case as precompute, trafo, adjoint. Fixed
 * rather than derived from a time budget, so both sides of a comparison run
 * the same shape. Each entry is worth about 20 ms per round before
 * BENCH_ITER_DIV.
 *
 * Per precision, because the default cutoff is not. Long double takes the
 * double counts.
 *
 * Re-derive after changing a geometry or the window: configure a walltime tree
 * of that precision, run
 *   bench_nfft_fast --benchmark_min_time=0.02s --benchmark_repetitions=3
 * and read the iteration column, with a floor of 4. */
#if defined(NFFT_SINGLE)
#define IT_1D_1024        400, 3200, 3400
#define IT_1D_8192         48,  230,  230
#define IT_1D_65536         6,   12,   13
#define IT_1D_8192_M1024  400,  330,  300
#define IT_1D_8192_M65536   6,   62,   67
#define IT_1D_SMALL_M      64,  260,  240
#define IT_1D_LARGE_M      35,  170,  190
#define IT_2D_32          190,  300,  480
#define IT_2D_128          12,    8,   10
#define IT_2D_256          12,    4,    4
#define IT_3D_8           260,   38,  100
#define IT_3D_16           64,   10,   20
#define IT_3D_32          124,    6,    7
#define IT_4D_8           420,    7,   13
#else
#define IT_1D_1024        320, 1600, 1450
#define IT_1D_8192         43,   51,   50
#define IT_1D_65536         5,    4,    4
#define IT_1D_8192_M1024  350,   66,   58
#define IT_1D_8192_M65536   5,   23,   20
#define IT_1D_SMALL_M      57,   57,   51
#define IT_1D_LARGE_M      29,   45,   42
#define IT_2D_32          170,   69,  120
#define IT_2D_128          10,    4,    4
#define IT_2D_256          10,    4,    4
#define IT_3D_8           230,   21,   13
#define IT_3D_16           57,    4,    4
#define IT_3D_32          115,    4,    4
#define IT_4D_8           350,    5,    4
#endif

/* A round only has to be long enough to time cleanly. Raise the divisor to buy
 * CI time, raise --benchmark_repetitions in .github/workflows/bench-linux.yml
 * to buy accuracy. */
#define BENCH_ITER_DIV 5
#define BENCH_MIN_ITERS 2

#define BENCH_ITERS(n) \
    ((n) / BENCH_ITER_DIV < BENCH_MIN_ITERS ? BENCH_MIN_ITERS \
                                            : (n) / BENCH_ITER_DIV)

#define REGISTER_CASE(tag, iters, ...) REGISTER_CASE_(tag, iters, __VA_ARGS__)

/* Trailing args are N[0..d-1], M, m. */
#define REGISTER_CASE_(tag, ipre, itrafo, iadj, ...) \
    BENCH_BUDGET(nfft_fast_precompute_psi_##tag, BENCH_ITERS(ipre))->Args({__VA_ARGS__}); \
    BENCH_BUDGET(nfft_fast_trafo_##tag, BENCH_ITERS(itrafo))->Args({__VA_ARGS__}); \
    BENCH_BUDGET(nfft_fast_adjoint_##tag, BENCH_ITERS(iadj))->Args({__VA_ARGS__});

/* 1d size sweep. */
REGISTER_CASE(1d, IT_1D_1024,       1024,         1024, DEFAULT_M)
REGISTER_CASE(1d, IT_1D_8192,       8192,         8192, DEFAULT_M)
REGISTER_CASE(1d, IT_1D_65536,      65536,       65536, DEFAULT_M)

/* 1d, off the M = N_total diagonal: separates FFT-phase from B-phase cost. */
REGISTER_CASE(1d, IT_1D_8192_M1024, 8192,         1024, DEFAULT_M)
REGISTER_CASE(1d, IT_1D_8192_M65536, 8192,       65536, DEFAULT_M)

/* 1d cutoff sweep: B-phase cost scales as (2m+2)^d. */
REGISTER_CASE(1d, IT_1D_SMALL_M,    8192,         8192, SMALL_M)
REGISTER_CASE(1d, IT_1D_LARGE_M,    8192,         8192, LARGE_M)

/* 2d size sweep. */
REGISTER_CASE(2d, IT_2D_32,         32, 32,       1024, DEFAULT_M)
REGISTER_CASE(2d, IT_2D_128,        128, 128,    16384, DEFAULT_M)
/* M is a quarter of the grid: at M = N_total one transform is too slow to time
 * in a short round. */
REGISTER_CASE(2d, IT_2D_256,        256, 256,    16384, DEFAULT_M)

/* 3d size sweep. */
REGISTER_CASE(3d, IT_3D_8,          8, 8, 8,       512, DEFAULT_M)
REGISTER_CASE(3d, IT_3D_16,         16, 16, 16,   2048, DEFAULT_M)
REGISTER_CASE(3d, IT_3D_32,         32, 32, 32,   1024, DEFAULT_M)

/* d = 4 uses the generic path instead of the specialized 1d/2d/3d kernels. */
REGISTER_CASE(4d, IT_4D_8,          8, 8, 8, 8,    256, DEFAULT_M)

BENCHMARK_MAIN();
