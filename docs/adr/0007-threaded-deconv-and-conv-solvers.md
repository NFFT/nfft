# ADR-0007: Threaded DECONV and CONV solvers in `plan_ng`

## Status
Accepted (2026-09-13).

## Context
ADR-0005 gave `plan_ng` FFTW's patience lattice, `NO_NONTHREADEDP`, and an
add-on threading library — but its own roster (`kernel/threads/conf.c`) was
empty: no threaded solver had been written, so `X(init_threads)` always
refused and `pl->nthr` could never rise above 1. ADR-0006 gave every solver a
composable operation count, but `Y(planner_mkplan)` still only ever measured
at the top level; a recursive call picked its winner by estimate alone.

This work fills both gaps: a threaded DECONV solver (Step A, every rank up to
`THREADS_MAX_RANK`) and a threaded CONV solver (Step C, every rank, forward
and adjoint), plus a real
recursive race in `Y(planner_mkplan)` so that a raised thread count now
affects every node of the plan tree, not only the one the guru measures
directly. `nfft_solver_fast_native` is reclassified from a serial solver to a
composer as part of the same change, since its own decline under
`NO_NONTHREADEDP` was blocking exactly the children this work adds.

## Decision

### DECONV (Step A)
- **Rows, not a flattened frequency index.** Legacy recovers each element's
  multi-index independently, paying `d` divisions per element. The threaded
  solver's `omp for` runs over *rows* (all axes but the last) and recovers a
  row's multi-index once, paying `d-1` divisions per row instead of `d` per
  element (`kernel/threads/deconv.c:run`, the `for (t = d - 2; t >= 0; t--)`
  loop). Rank 1 has one row and distributes its frequencies directly instead.
- **The clear is inside the parallel region.** Forward zeroes the untouched
  grid before writing (the oversampling gap DECONV never fills), as its own
  `omp for` with an implicit barrier ahead of the write loop, never with
  `nowait`. In rank 1 this clears only the untouched run
  (`[Nposl, nl - Nnegl)`); at rank $\geq 2$ it is cheaper to clear the whole
  grid in one pass than to compute the gap's shape per row. Legacy clears the
  whole grid single-threaded before any DECONV work starts.
- **No clear in the adjoint**, matching `deconv-nd.c`'s serial solver: the
  adjoint only ever reads grid cells DECONV's own forward would have written,
  so there is nothing to clear on that side either.

### CONV (Step C)
- **Bands, not atomics.** Legacy offers both accumulation strategies for its
  adjoint; the atomic variant makes the grid's accumulation order depend on
  the OpenMP schedule, which is the source of the run-to-run ULP drift
  documented for the legacy accuracy report. The threaded solver instead
  partitions the grid's axis 0 into disjoint bands, each written by exactly
  one band-owning iteration, giving bit-exact reproducibility at any thread
  count — the plan's own tests (`compare_conv` in `tests/ng_omp.c`) assert
  `memcmp`-equality against the serial solver at 2, 4 and 7 threads.
- **A counting sort by band, not a radix sort of the linearised window
  start.** The adjoint needs, per band, the list of nodes whose window
  touches it — a locality order is not what it needs. Legacy's radix sort
  over the flattened window start (plus a binary search) reaches a stronger
  result than the adjoint requires; `fill_bands`'s counting sort by
  `u[j*d] / band` is the part of that machinery this problem actually needs,
  and it is stable, so each band's list stays ascending in node index `j` —
  which is what lets the merge below reproduce the serial accumulation order.
  Sorting nodes for cache locality is ADR-0008.
- **Band ownership follows the `omp for`'s iteration variable, not the
  thread id.** `apply_adjoint`'s band loop indexes bands `0..nb-1` directly
  (`schedule(dynamic, 1)`); a team smaller than `nb` still visits every band
  exactly once, just with some threads taking more than one.
- **Bands are at least `2m+2` rows wide.** A node's window support spans at
  most `len = 2m+2` grid rows, so from a start in band `b` it can only reach
  band `b` itself or its cyclic predecessor — never further — provided every
  band is at least `len` wide. `nb = min(n0/len, 4*nthr)` guarantees
  `band >= len` by construction; merging a band's own node list with its
  predecessor's, both ascending in `j`, is then exhaustive. Below two bands
  (`nb == 1`, i.e. the axis is too narrow to give every thread a band this
  wide) the adjoint's band loop degenerates to one band and runs on a single
  thread; the forward direction still uses the whole team regardless, since
  it has no such constraint.
- **The shared per-node kernels.** `kernel/conv/conv.h`'s
  `Y(conv_node_1d/2d/3d/nd)` are the same tap-accumulation code the serial
  CONV solvers already used, now factored out so the threaded solver calls
  them too — this is what makes bit-exact reproduction possible in the first
  place; a second, independently-written kernel could not have promised it.
  Lifting axis 0 out of the nd solver's flat multi-index increment loop (needed so the
  shared kernel could take a `[lo, hi)` clip on axis 0 alone) moved rank
  $\geq 4$ results by ULPs relative to the previous nd loop structure — an
  accepted, already-recorded rounding change, not a new one introduced here.

### Both
- **They decline below a two-thread budget** (`pl->nthr < 2`), exactly as
  ADR-0005's serial solvers decline above it: linking the add-on library and
  raising the count to 1 changes no plan until the caller actually asks for
  more than one thread.
- **No tuning constant anywhere.** An operation count cannot express a
  thread count, so both threaded solvers report the *serial* op counts
  (`kernel/threads/deconv.c`, `kernel/threads/conv.c` copy the serial
  formulas verbatim), plus, for CONV, an honest `other` term for the
  straddle cost a node crossing a band boundary pays (its per-node run setup
  runs twice; its tap work is split, not repeated). Below `NFFT_PATIENT` the
  roster alone decides which solver runs, because `NO_NONTHREADEDP` makes
  every serial DECONV/CONV solver decline once `pl->nthr > 1` — no cost
  comparison happens, there is nothing to compare against once the serial
  candidate is gone. At `NFFT_PATIENT` and above, both compete and
  `Y(planner_mkplan)`'s race measures real elapsed time, exactly as the
  top-level guru already did for the whole NFFT problem (ADR-0005) — this is
  the recursive race this ADR's title work refers to. **Divergence
  from FFTW, recorded here as it was found and decided during this work's
  Task 7:** a lone surviving candidate — whether because it was the only
  applicable solver, or because every other candidate failed the
  operation-count prune gate — is adopted *untimed*; its `pcost` stays at
  its estimate (rank 0 in `Y(iestimate_cost)`, effectively). FFTW measures a
  lone candidate anyway, because its `pcost` field starts at exactly 0.0 and
  its `evaluate_plan` always measures a plan whose cost reads as zero. Here
  the guru's own top-level lone-survivor path already did this before this
  work started; Task 7 made the same choice for the *recursive* case rather
  than introduce an inconsistency between the two call sites.
- **`nfft_solver_fast_native` is a composer.** It runs no loop of its own —
  `apply`/`apply_adjoint` delegate entirely to its DECONV child, an internal
  FFTW plan, and its CONV child. Before this work it still tested
  `NO_NONTHREADEDP` itself and declined below `NFFT_PATIENT` whenever
  `pl->nthr > 1`, which was simply wrong: it is not a serial solver with a
  loop to parallelise, so the condition the guard exists to express (this
  solver does not know how to use more than one thread) was never true of
  it. The guard's effect was to make the *entire* fast NFFT unavailable
  below `NFFT_PATIENT` the moment a caller raised the thread count, which
  defeated the point of adding threaded DECONV/CONV children at all: with
  `nfft_solver_fast_native` declining, its children were never even reached.
  Removing the decline lets the composer answer at every patience level; its
  children's own `NO_NONTHREADEDP` guards (now genuinely applicable, since
  they are genuinely serial) are what continue to gate the serial-vs-threaded
  choice one level down, exactly where that decision belongs.
- **Cost is measured forward only**, an existing, explicitly accepted
  limitation this work inherits and does not revisit: one problem, one plan,
  one race. `prob[ADJ]` is never constructed for the NFFT kind, and
  `execute_adjoint` runs the forward plan's own `apply_adjoint`.
  `Y(plan_measure_cost)` times `apply` (the forward direction) exclusively,
  so the whole tree — including which of a serial or threaded DECONV/CONV
  child wins — is chosen on forward cost and then serves the adjoint
  unconditionally. This is unchanged by adding threaded children: they are
  raced by the same forward-only measurement as everything else. Separate
  per-direction solvers are the way out, and the `dir[ADJ]` scaffolding
  already present in the guru (unused today) is where they would land.
- **The measured scaling tables.** See below.

### The FP-contraction fix in `kernel/conv/conv.h` (required addition)
Task 6 (the threaded CONV adjoint) found that its own given implementation,
built unmodified with this project's default optimized flags
(`AX_CC_MAXOPT`: `-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math`,
GCC 13.3.0, aarch64), did not reproduce the serial CONV solver's output
bit-exactly. Every mismatch was 1-2 ULP, on cells the merge logic provably
computes from the same nodes in the same order as the serial loop — ruling
out a logic error — and disappeared entirely when the affected translation
units were rebuilt with `-ffp-contract=off`, appearing again when that flag
was removed, reproduced twice. The mechanism: under `-ffast-math` GCC ignores
the standard `#pragma STDC FP_CONTRACT OFF` and fuses a tap's
`g[idx] += fj * psi[idx]` into a single rounded FMA or not depending on
whether the call's clip bound (`lo`) is a compile-time literal `0` (every
existing serial call site) or a value computed at runtime (the threaded
adjoint's per-band `lo = b * band`) — the same source line, compiled two
different ways depending on what its caller passes.

The fix, added to `kernel/conv/conv.h` immediately after its include guard
and still in the tree:

```c
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC optimize("fp-contract=off")
#endif
```

Because it is a `#pragma` rather than a build flag, it applies to *every*
translation unit that includes the header — every serial CONV solver
(`conv-1d/2d/3d/nd.c`) as well as the threaded one — confirmed during Task 6's
code review by objdump: without the pragma those four files' compiled tap
call sites contain 8, 11, 11 and 1 FMA instructions respectively; with it,
zero, across all of them. `#if defined(__GNUC__) && !defined(__clang__)`
scopes it to GCC specifically. Clang is also detected as GNU-vendor by this
project's `m4/ax_compiler_vendor.m4` and so also receives `-ffast-math` from
`AX_CC_MAXOPT`; Task 6's code review built and tested clang 18.1.3 (aarch64)
without the pragma and found it fuses consistently across call shapes
regardless of whether `lo` is a literal or a runtime value — no divergence,
and so no gap left by scoping the pragma to GCC. This is one compiler
version on one architecture, empirically checked, not a guarantee that every
clang release behaves the same.

**The measured cost of the fix (this task).** A throwaway micro-benchmark
(not committed) called the shared rank-3 tap kernel directly
(`Y(conv_node_3d)`, forward and adjoint, `m=6` so `len=14`, a `128^3`
oversampled grid, `M=200000` synthetic nodes) against two copies of
`conv.h` — the shipped one, and a copy with the pragma's body removed —
compiled with the project's exact `AX_CC_MAXOPT` flags. objdump confirmed
the expected mechanism (0 vs 29 FMA instructions in the resulting binaries).
Six interleaved runs (alternating pragma / no-pragma, each an internal
median of 3) gave, taking the median across the six run-medians:

| direction | pragma (`fp-contract=off`) | no pragma (FMA allowed) | cost of the fix |
|---|---|---|---|
| forward | 0.4845 s | 0.4470 s | +8.4% |
| adjoint | 0.6185 s | 0.6090 s | +1.6% |

Run-to-run variance across the six interleaved samples (roughly 0.37-0.62 s
per side) was comparable in magnitude to the effect itself, so these
percentages should be read as "a small, single-digit-percent cost on this
host," not as a precise constant. What the measurement rules out is the
naive expectation of up to a 2x cost from banning FMA in a
multiply-accumulate-shaped loop: the tap kernel's actual bottleneck on this
host is the strided gather/scatter against the oversampled grid, not the
multiply-add itself, consistent with ADR-0006's observation that these
kernels are more bandwidth-bound than compute-bound. Only rank 3, double
precision, one host, one problem size was checked — the mechanism (a
`#pragma` covering every TU that includes the header) is rank- and
precision-independent by construction, so the *presence* of some cost at
every rank and precision is expected, but this measurement does not license
extrapolating the specific percentage to rank 1, rank $\geq 4$, or to float
or long double, where the tap loop's arithmetic-to-memory-traffic ratio
differs.

**The tradeoff, stated plainly.** This is a deliberate correctness-for-
performance decision: every serial GCC-built CONV solver, at every
precision, gives up FMA fusion in its hottest loop so that the threaded
adjoint's band-clipped calls reproduce the serial accumulation bit-exactly.
GCC's own documentation for the `optimize` pragma and its function-attribute
equivalent does not treat it as production-safe in general — it is not
guaranteed to produce correct or even consistently-defined code for
arbitrary optimization strings, and its exact effect has changed across GCC
releases, unlike a stable command-line flag. Nothing in this project's build
currently threatens that in practice (the pragma turns one specific,
GCC-versioned string off, not on, and both GCC 13's frontend and the two
reviewers who independently re-derived this fix on this dev container agree
on what `fp-contract=off` does here), but a future GCC could in principle
change that pragma's behaviour in a way `-ffp-contract=off` as a `CFLAGS`
entry would not be exposed to. This is the price of putting the fix in the
header rather than in the build system: it travels with the source (correct
wherever `conv.h` is compiled, including out-of-tree consumers who build
with their own flags) but rests on a mechanism GCC documents as less stable
than a compiler flag.

**Decision (final whole-branch review).** Keeping the pragma, and keeping
`compare_conv`'s `memcmp`-equality assertion exactly as written, is a
deliberate choice, not an incidental side effect of Task 6's debugging. R11
is literal and binding: "The adjoint result is exactly equal to the serial
result, at every thread count and every band count." Bit-equality between
the serial and threaded CONV solvers on GCC is preserved because R11
requires it in those words, at the measured cost above — roughly 1-8% on
the forward and adjoint tap kernels on this host, GCC 13.3.0, aarch64 (see
the table above). The rejected alternative, raised during review, was to
drop the pragma and relax `compare_conv`'s comparison to a small
relative-ULP tolerance against the serial solver, asserting only
thread-count invariance instead of serial-vs-threaded bit-equality. That was
rejected for two reasons: it would directly contradict R11's own "exactly
equal" wording, and it would make this the only approximate-exactness test
in a suite that is otherwise exact or `memcmp`-based throughout — R5's
thread-count-invariance test and R11's own CONV test both compare by
`memcmp` today. Relaxing one of the two to a tolerance would introduce an
inconsistent precedent into the same test file rather than follow the one
the PRD already set. This decision is not architecturally locked in: if the
project later decides the performance cost outweighs literal bit-equality,
dropping the pragma and relaxing the test to a tolerance is a small,
localised change, confined to `kernel/conv/conv.h` and `compare_conv`.

### ThreadSanitizer run (final whole-branch review)

PRD.md's acceptance criterion asks for "one ThreadSanitizer run against
clang's `libomp` during the adjoint task, recorded in the ADR rather than
wired into CI — `libgomp` is not instrumented and reports false positives on
every OpenMP barrier." This was attempted here, on the dev container
(clang 18.1.3, aarch64, Ubuntu 24.04). The short version: clang does default
to `libomp`, not `libgomp`, as the PRD assumed it would — but the underlying
worry the PRD raised for `libgomp` turns out to apply to `libomp` too, for
the same reason (barrier synchronisation invisible to TSan's happens-before
model), and no `libomp` build available in this environment — the distro
package or one built from clang 18's own matching release source with
`-fsanitize=thread` — comes back clean for even a trivial OpenMP program.

**Confirming clang's default runtime.** `clang -fopenmp -fsanitize=thread -###
-o /tmp/x /dev/null` shows the link line pulling
`libclang_rt.tsan-aarch64.a` and `-lomp` — clang's own driver links `libomp`
for `-fopenmp` on this distribution, never `libgomp`, so the PRD's premise
("clang's libomp") holds without any special flag.

**Attempt 1: the distribution's `libomp.so.5` (`libomp5-18`, not built with
TSan instrumentation).** A trivial `#pragma omp parallel for
reduction(+:s)` program, compiled `clang -fopenmp -fsanitize=thread -O1 -g`,
reports a data race before a single line of user code runs:
```
WARNING: ThreadSanitizer: data race (pid=...)
  Atomic read of size 1 ... by main thread:
    #0 pthread_mutex_lock ...
    #1 <null> <null> (libomp.so.5+0xb7064)
  Previous write of size 1 ... by thread T8:
    #0 pthread_mutex_init ...
    #1 <null> <null> (libomp.so.5+0xb3aa8)
```
This is a mutex-reuse false positive inside `libomp`'s own thread-pool
bookkeeping (`pthread_mutex_init` on memory a pooled worker's prior mutex
occupied, raced against a later `pthread_mutex_lock` from the same pool) —
`libomp.so.5` itself carries no TSan instrumentation, so its internal
synchronisation is invisible to the tool from the outside. A suppression
targeting `libomp.so.5` by name (`called_from_lib:libomp.so.5`,
`race:libomp.so.5`) silences it, but a second experiment — the same
technique applied to a program with a *genuine*, deliberately planted user
race (`g[i % 4] += i` across an unsynchronised `omp for`) — also went
silent under that suppression. A module-name suppression broad enough to
hide `libomp`'s internal noise is also broad enough to hide a real user race
whose call stack passes through `libomp` (which every OpenMP-parallel write
does, via the fork/dispatch frames), so this route was abandoned as unsound
rather than merely noisy.

**Attempt 2: a `libomp` built from source with `-fsanitize=thread`.** To
rule out "the distro package just isn't instrumented" as the whole story,
`openmp-18.1.3.src.tar.xz` (matching the installed clang 18.1.3 exactly) was
fetched and built standalone (`cmake ... -DCMAKE_C_COMPILER=clang
-DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_FLAGS=-fsanitize=thread
-DCMAKE_CXX_FLAGS=-fsanitize=thread -DLIBOMP_ENABLE_SHARED=ON
-DOPENMP_ENABLE_LIBOMPTARGET=OFF`), producing a genuinely TSan-instrumented
`libomp.so`, confirmed by TSan's own reports resolving into the fetched
`openmp-src/runtime/src/*.cpp` sources by file and line once loaded via
`LD_LIBRARY_PATH` (with a `libomp.so.5` symlink so the dynamic loader
actually picks it over the system copy). Even with a real,
source-instrumented runtime, the trivial reduction program above still
races — this time inside the runtime's own barrier and thread-launch code:
```
WARNING: ThreadSanitizer: data race (pid=...)
  Write of size 8 ... by main thread:
    #0 __kmp_allocate_thread kmp_runtime.cpp:4559
    #1 __kmp_allocate_team kmp_runtime.cpp:5390
    ...
  Previous read of size 8 ... by thread T1:
    #0 __kmp_check_stack_overlap kmp_runtime.cpp:321
    #1 __kmp_launch_worker(void*) z_Linux_util.cpp:562
```

**Running the actual target: `checkall_ng_omp` under the from-source
instrumented `libomp`.** The project was rebuilt with `CC=clang
CFLAGS="-g -O1 -fsanitize=thread -fno-omit-frame-pointer"
LDFLAGS="-fsanitize=thread"` (overridden on the `make` command line, since
this project's `configure.ac` runs `AX_CC_MAXOPT` regardless of an
inherited `CFLAGS` — an unrelated, pre-existing autoconf-generated-script
quirk, `ac_test_CFLAGS` compares against a literal that the current
autoconf version no longer produces), then `tests/checkall_ng_omp` was run
directly (CUnit drives every registered case, including
`check_ng_omp_conv_adjoint` and its narrow-grid/clustered/wrap siblings, at
2, 4 and 7 threads internally) with the instrumented `libomp.so`
preloaded via `LD_LIBRARY_PATH`. The run reports on the order of a thousand
warnings (983 in one representative run), but inspecting where the two
racing accesses actually happen shows the same pattern throughout: every
block's top frame is either (a) literally inside `libomp`'s own runtime —
`kmp_barrier.cpp`, `kmp_runtime.cpp`, `kmp_wait_release.h`,
`z_Linux_util.cpp` — including several that are TSan flagging the barrier's
*own release operation* (`__kmp_hyper_barrier_release`,
`flag_traits<...>::test_then_add4`) as racing against a worker spinning on
the same flag, which is the barrier mechanism working as designed, not a
bug; or (b) a legitimate serial read in the test harness (e.g.
`compare_rank`'s post-`run_once` comparison loop, `tests/ng_omp.c:141`) of
memory a worker thread wrote *inside* the now-finished parallel region
(`kernel/threads/deconv.c`'s or `kernel/threads/conv.c`'s outlined
lambda) — i.e. TSan reporting a race between a parallel write and the very
serial read the `#pragma omp parallel` region's implicit join barrier exists
to make safe. Category (b) is the "false positives on every OpenMP barrier"
failure mode the PRD named for `libgomp`, reproduced here against `libomp`,
including one built from source with instrumentation: `libomp`'s spin-wait
barrier releases a worker via a plain atomic flag flip
(`kmp_wait_release.h`) without the explicit TSan happens-before annotation
(`__tsan_acquire`/`__tsan_release`) that would be needed for the tool to
recognise it as a synchronisation point, so every read that follows a
parallel write — which is every test in this suite — is flagged, drowning
out anything a real race in `kernel/threads/conv.c` or `deconv.c` would look
like.

**What this does and does not show.** It does not deliver a clean run — no
`libomp` available in this environment, distribution-packaged or built from
the matching release source with `-fsanitize=thread`, currently comes back
clean for *any* OpenMP program under TSan on this host, before any NFFT code
runs at all. It does confirm the PRD's premise that clang links `libomp` by
default, and it does show, by direct inspection of every reported race's
call stacks across a real run of `checkall_ng_omp`, that none of the
roughly-a-thousand reports point at a plausible genuine cross-thread
conflict over shared NFFT data (`g`, `f_hat`, `f`, the band/`u`/`psi`
tables) — every one resolves to either `libomp`'s own internals or the
barrier-invisibility pattern above. TSan detects races on instrumented code
its happens-before model can see; it cannot currently see past `libomp`'s
own barrier implementation on this toolchain, instrumented or not, which
means a TSan run in this environment cannot today distinguish "no race" from
"a race hidden behind the same barrier that hides all this noise." The
corroborating evidence for actual thread-safety is the non-TSan evidence
this branch already has in quantity: R5's and R11's `memcmp`-exact
comparisons at 2, 4 and 7 threads, run for real (not simulated) in
`make check` across double, float and long-double builds (see the test
results recorded elsewhere in this ADR and in the branch's own task
reports), which is precisely the kind of check a genuine race would be
likely to eventually break and has not.

## The measured scaling tables
Throwaway harness (not committed; scratchpad only), private planners built
directly against `Y(mkproblem_deconv)`/`Y(mkproblem_conv)`, one at `nthr = 1`
and one at `nthr = T` with `PLNR_NO_NONTHREADED` forced so the threaded
solver is the one actually exercised at every `T`. 200 `apply` /
`apply_adjoint` calls per cell, median of 3 runs, idle dev container
(aarch64, 10 cores). One representative geometry per rank (not tuned for a
best case): DECONV `d=1`: `N=2^18, n=2^19`; `d=2`: `N=512^2, n=1024^2`;
`d=3`: `N=96^3, n=192^3`. CONV `d=1`: `n=2^14, M=100000`; `d=2`:
`n=192^2, M=30000`; `d=3`: `n=64^3, M=8000`. All times in seconds for 200
calls.

**DECONV, forward (`apply`):**

| d | T=1 | T=2 | T=4 | T=8 |
|---|---|---|---|---|
| 1 | 0.0987 | 0.0954 | 0.0575 | 0.0386 |
| 2 | 0.0395 | 0.0528 | 0.0374 | 0.0639 |
| 3 | 0.4858 | 0.6286 | 0.5641 | 0.5546 |

**DECONV, adjoint:**

| d | T=1 | T=2 | T=4 | T=8 |
|---|---|---|---|---|
| 1 | 0.0212 | 0.0110 | 0.0074 | 0.0130 |
| 2 | 0.0223 | 0.0108 | 0.0061 | 0.0096 |
| 3 | 0.2501 | 0.0771 | 0.0702 | 0.0885 |

**CONV, forward, uniform nodes:**

| d | T=1 | T=2 | T=4 | T=8 |
|---|---|---|---|---|
| 1 | 0.1067 | 0.0765 | 0.0466 | 0.0634 |
| 2 | 0.3020 | 0.1577 | 0.0887 | 0.0922 |
| 3 | 1.6321 | 0.9022 | 0.5375 | 0.3911 |

**CONV, adjoint, uniform nodes:**

| d | T=1 | T=2 | T=4 | T=8 |
|---|---|---|---|---|
| 1 | 0.1626 | 0.1756 | 0.1101 | 0.1298 |
| 2 | 0.4701 | 0.3094 | 0.2121 | 0.1698 |
| 3 | 2.6798 | 1.1045 | 0.6820 | 0.7167 |

**CONV, forward, nodes clustered into a 2% slice per axis:**

| d | T=1 | T=2 | T=4 | T=8 |
|---|---|---|---|---|
| 1 | 0.1007 | 0.0602 | 0.0337 | 0.0446 |
| 2 | 0.2998 | 0.1619 | 0.0948 | 0.0949 |
| 3 | 1.5918 | 0.9660 | 0.5980 | 0.4174 |

**CONV, adjoint, nodes clustered into a 2% slice per axis:**

| d | T=1 | T=2 | T=4 | T=8 |
|---|---|---|---|---|
| 1 | 0.1455 | 0.1711 | 0.1557 | 0.1764 |
| 2 | 0.5217 | 0.4244 | 0.3383 | 0.3557 |
| 3 | 2.3059 | 1.1730 | 1.1579 | 1.1799 |

**What the tables show.** CONV's forward direction — a plain `omp for` over
independent nodes — scales the most cleanly and identically whether nodes
are uniform or clustered (as expected: the forward gather does not care
where a node's support lands relative to other nodes), reaching roughly
4-4.5x at `T=8` for `d=3`. The adjoint tells the more interesting story,
which is what "what band oversubscription bought" (the brief's phrase) means
concretely: on **uniform** nodes the `d=3` adjoint scales similarly to the
forward (2.4x at `T=4`) because nodes spread themselves across most of the
`nb = min(n0/len, 4*nthr)` oversubscribed bands. On **clustered** nodes —
all `M` nodes packed into a 2% slice of one axis — the `d=3` adjoint's
speedup collapses after `T=2` (1.97x, then flat: 1.99x and 1.95x at `T=4`
and `T=8`) because virtually every node now falls into the same one or two
physical bands regardless of how many bands the team is nominally divided
into; the busiest band's serial work sets the floor no team size below "one
band per node" can beat. Oversubscription (four bands per thread, per
`kernel/threads/conv.c`'s comment) buys resilience against moderate,
realistic skew — it is why the uniform case does not need every node placed
by hand to scale — but it cannot rescue an adversarial worst case where the
geometry itself concentrates all the work in one place; that ceiling is
inherent to spatial partitioning, not a flaw in the band mechanism.

DECONV is smaller relative to its own per-call overhead at these sizes and
shows the opposite lesson: `d=1` and `d=2` barely benefit or even regress at
low thread counts (own-plan/OpenMP-team startup cost dominating an
already-cheap diagonal scale), and `d=3` forward *regresses* at `T=2`
(0.486s -> 0.629s) before flattening — consistent with DECONV forward being
almost pure memory traffic (a zero-fill pass over the whole oversampled grid
plus one more pass writing it) so that adding threads adds contention for
the same memory channels without adding usable parallelism once bandwidth,
not the loop itself, is the bottleneck. This is not a regression introduced
by this work in the sense of a defect — the serial solver is simply
narrower than DECONV's threaded solver is built to help with at these
sizes; a caller planning at `NFFT_PATIENT` gets the real measured cost and
the race picks the serial solver in that regime automatically.

## Deliberate divergences from FFTW
1. **A lone surviving candidate is adopted untimed**, at every level the
   search recurses into now, not only at the guru's own top level. FFTW's
   `pcost` starts at exactly `0.0` and its `evaluate_plan` always measures
   a plan whose cost reads as zero, so FFTW effectively measures every fresh
   candidate including a lone one; here `Y(iestimate_cost)` already gives a
   lone candidate a real, nonzero, comparable cost, so nothing is gained by
   spending a measurement on it. See the "Both" section above for why this
   consistency (not adding a special case for the recursive call this work
   introduces) was the deliberate choice made in Task 7.
2. **`nfft_solver_fast_native`'s reclassification as a composer** has no
   FFTW analogue to diverge from as such, but is recorded here as a
   correction to this project's own prior modelling (ADR-0005 treated every
   NFFT-kind solver, including this one, as subject to the same
   `NO_NONTHREADEDP` test): FFTW has no single "solver" that is purely a
   dispatcher over the same kind's own child search the way this composer
   recurses into DECONV/CONV, so there was no precedent to follow either
   way; the fix was to recognise the category error rather than to import a
   rule from FFTW.

## Consequences
- **Raising the thread count now does something below `NFFT_PATIENT`.**
  Before this work, `pl->nthr > 1` below `NFFT_PATIENT` made every
  NFFT-kind solver decline (including, incorrectly, the composer itself),
  so a raised count was inert everywhere except at `NFFT_PATIENT` and above,
  and even there it had nothing threaded to race against. Now it is not: the
  serial DECONV/CONV solvers step aside for their threaded counterparts at
  every level below `PATIENT`, and race them at `PATIENT` and above.
- **CONV's atomic-free design removes a source of nondeterminism that
  legacy has**: the run-to-run ULP drift in the legacy accuracy report
  traced to atomic accumulation order is structurally impossible for the
  threaded CONV solver, at any thread count, by construction (bands, not
  atomics) — verified by the plan's own `memcmp` tests, not merely argued.
- **Every serial GCC-built CONV solver loses FMA fusion in its hottest
  loop**, at every precision, as the price of that same bit-exactness
  guarantee reaching the threaded solver's band-clipped calls. Measured
  cost on this host: single-digit percent (see above), not the naive 2x a
  compute-bound reading of "no more FMA" would suggest, because the kernel
  is bandwidth-bound in practice. Clang needs no such guard on the one
  version tested.
- **Wisdom entries for DECONV/CONV problems now sometimes point at a
  threaded solver.** A wisdom file exported by a program that used more
  than one thread stores that fact implicitly, through which solver won;
  importing it into a program that has not linked the threading add-on (and
  so has no threaded solver registered) makes that specific entry a clean
  cache miss on replay, per the existing configuration-signature mechanism
  (ADR-0005) — the roster itself is part of that signature, so linking or
  not linking the add-on already partitions wisdom correctly without any
  new mechanism.

## Known gaps
- `NFFT_NO_FAST_NATIVE` below `NFFT_PATIENT` with a raised thread count
  still returns `NULL` — forbidding the only fast solver leaves the guru
  with no plan to serve a request that also insists on more than one thread
  being used, since the direct NDFT solvers are genuinely serial and still
  decline.
- The internal FFTW child plan inside `nfft_solver_fast_native` still runs
  at one thread regardless of `pl->nthr` — nothing in this work touches
  FFTW's own thread count.
- **Task 7's timelimit-division-among-survivors plumbing in
  `Y(planner_mkplan)` (the `budget / nsurvive` divide around the race loop)
  is currently a no-op for DECONV and CONV**, since both are leaf problem
  kinds with no further recursive `Y(planner_mkplan)` call inside their own
  `mkplan` to inherit the divided budget, and the race loop itself enforces
  no timeout of its own — not a defect, just dormant until a threaded
  solver kind can itself recurse into a further child race.
- **`check_ng_omp_nfft_all_levels` (Task 8, `tests/ng_omp.c:522`) fails in a
  float `--enable-openmp` build**, found while running this task's own
  full-matrix sanity check. Its per-node correctness assertion
  (`tests/ng_omp.c:573`) uses a fixed absolute tolerance,
  `K(1.0e-12) * (CABS(f_ref[j]) + K(1.0))` — a bound written with only
  double precision in mind. `1.0e-12` is roughly five orders of magnitude
  below float's own machine epsilon (~1.19e-7), so the assertion cannot
  pass in a float build even for a correct result: any single-ULP rounding
  difference between the reference call (built at `nthr = 1`, before
  `X(plan_with_nthreads)(4)` is called) and the loop's `nthr = 4` calls —
  expected and harmless, since a different DECONV/CONV solver answers below
  `NFFT_PATIENT` — already exceeds it at float precision, while comfortably
  clearing it at double or long-double precision (confirmed: double and
  long-double both pass this test cleanly in this task's matrix run; only
  float fails, and only this one assertion). This is a test-tolerance
  defect in Task 8's test, not a solver correctness defect, and parallel in
  kind to the codegen-sensitivity that got `check_nplan_ndft_accuracy`
  removed (ADR-0006): that test compared two implementations of the same
  recurrence at round-off level and codegen flipped the comparison — but
  this is a *different* test, not previously documented, and this task's
  brief scopes only `docs/adr/` and `.claude/skills/` changes, so it is
  recorded here rather than fixed.

## Not decided here
- Racing several thread counts against each other, or searching for a good
  one, rather than taking `pl->nthr` as an explicit, caller-given parameter
  (unchanged from ADR-0005).
- Threading FFTW's own internal child plan inside `nfft_solver_fast_native`.
- Per-direction solvers that would let the adjoint race on its own measured
  cost rather than inherit the forward winner (unchanged from ADR-0005 and
  ADR-0006).
- Making the timelimit-division-among-survivors plumbing actually bind for
  DECONV/CONV — it needs a threaded solver kind that itself recurses into a
  further child race, which none of today's do.
