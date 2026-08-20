# Review: `feature/planner` vs `develop`

Date: 2026-08-21. Scope: the full changeset (`git diff develop...HEAD`,
126 files, +39,528 / -273 lines, one squashed commit `a03a26b8`).

Method. Read of every new or changed C source and header under `include/`,
`kernel/planner/`, `kernel/nfft/`, `kernel/deconv/`, `kernel/conv/`,
`kernel/util/`, the tests, examples, build files and docs. Verification on
this host: clean double build (`--enable-all --enable-openmp --enable-tests`),
float build, `--enable-debug` (ASan/UBSan) build, `tests/checkall_ng` and
`checkall_ng_threads` runs in all three, a CMake configure and build,
`make distdir`, valgrind on `examples/nfft/nfast_native`, `nm` on the built
libraries, clang-format dry runs, an md5 test-vector program, and small
programs against the public API (wisdom export, plan print, planning time).
Every finding below cites a file and line and was checked against the code.

Verdict in one line: the architecture is the right one and the numerics are
correct, but the branch is not at merge quality. Section 1 lists what must
change before merge. Sections 2-8 list what I would have done differently,
ranked. Section 9 has the build-matrix evidence. Section 10 is a suggested
order of work.

Ranking: **blocker** = must fix before merge. **major** = fix before merge
or file an issue with an owner. **minor** = cleanup pass. **nit** = taste.

---

## 0. What is good

- The problem / solver / plan trinity with const vtables, `mkplan` returning
  `NULL` for "not applicable", `awake()` for node-dependent tables, and a
  data-blind 128-bit size-class key is the right shape for a self-tuning
  NFFT. It composes: `fast_native` recurses into the planner for DECONV and
  CONV children.
- `x` copied at construction with unit-axis elision, `f_hat`/`f` borrowed,
  `_on` new-array variants: a clean ownership story.
- Wisdom import is atomic (snapshot and rollback), rejects on configuration
  signature mismatch, and cannot be made to yield a wrong plan from a
  crafted file. md5 verified against all RFC 1321 vectors.
- Numerics verified by read: the blocked 1D recurrence reseeds per block
  with an FMA-reduced phase (error does not grow with N); index ranges are
  right for even, odd, type-II and unit axes; adjoints mirror forwards; the
  Kaiser-Bessel self-normalisation is exact in theory and cancellation-free
  in code; the fast path is correct on the domain `guards_ok` admits.
- 96 new CUnit tests, 0 failures in double, float, OpenMP and ASan/UBSan
  builds, 1.3 s wall. Valgrind: 0 errors, 0 bytes definitely lost.
- New code compiles without `-Wconversion`/`-Wsign-conversion` noise under
  `--enable-debug`; the only new warning there is one unused parameter.
- `.gitignore` `nfft-*` -> `/nfft-*` and the `ticks.h.in` stray-semicolon fix
  are correct and necessary.

---

## 1. Blockers

### 1.1 CMake build does not link

**Status 2026-08-22: fixed.** `conv-1d.c`/`deconv-1d.c` added to `kernel/CMakeLists.txt`; `planner.*`/`nplan.*` removed from the legacy `checkall` list; `example_nfft_nfast_native` target added. CMake build links; `ctest -R checkall_ng` passes 2/2.

`kernel/CMakeLists.txt:15-16` list `conv-2d/3d/nd.c` and
`deconv-2d/3d/nd.c` but not `conv/conv-1d.c` and `deconv/deconv-1d.c`,
which Autotools lists (`kernel/conv/Makefile.am:4`,
`kernel/deconv/Makefile.am:4`) and the rosters reference
(`kernel/conv/solver.c:30`). A CMake build yields
`U nfft_conv_solver_1d_register` in `libnfft3.so` and every executable fails
to link. The CodSpeed benchmark CI (`.github/workflows/bench-linux.yml`)
uses this build.

`tests/CMakeLists.txt:11-12` adds `planner.c nplan.c` to the legacy
`checkall` sources; `nplan.c` needs `nplan_perm.c`, which is only in the
`checkall_ng` list, so `checkall` fails to link too. Autotools does not make
this mistake (`tests/Makefile.am:38`).

Fix: add the two `-1d.c` files to `NFFT_KERNEL_SOURCES`; remove
`planner.*`/`nplan.*` from `NFFT_TEST_SOURCES`; add a CMake target for
`examples/nfft/nfast_native.c` (parity with Autotools); run the CMake lane
in CI on every PR, not only the bench lane.

### 1.2 `make dist` tarball cannot compile

**Status 2026-08-22: fixed.** `ndft.h` added to `libnfft_la_SOURCES`; empty `xcheck.h` deleted with its two includes. `make distdir` verified.

`kernel/nfft/Makefile.am:10` lists only `.c` files; `ndft.h` and `xcheck.h`
are in no `_SOURCES` or `EXTRA_DIST`. `make distdir` confirms both missing.
Fix: append them to `libnfft_la_SOURCES` as `conv/Makefile.am` does with
`conv.h` (or delete `xcheck.h`, see 5.2).

### 1.3 Measured mode (the default) stalls for a minute on ordinary sizes

**Status 2026-08-22: fixed (option 1).** `PLNR_PRUNE_RATIO = 8.0` (`iplanner.h`); the race in `kernel/nfft/plan.c` skips candidates whose `pcost` exceeds 8x the cheapest estimate and adopts a lone survivor untimed. Same 2D case: guru 49.1 s -> 0.001 s. Regression test `check_nplan_measured_prunes_by_estimate` (slow test solver with `pcost = 1e18` must never be applied). Two tests that relied on the blocked NDFT being timed at N=64, M=8192 now use M=64 so two candidates survive the gate.

`NFFT_MEASURE == 0` is the default. The race times every applicable
candidate, including the O(N*M) direct NDFT, with no pruning by analytic
cost. Measured on this host, 2D, N = 256x256, n = 512x512, M = 65536, m = 6,
`-O3`, double:

| planning flag | guru | precompute | execute |
|---|---|---|---|
| `NFFT_ESTIMATE` | 0.001 s | 0.035 s | 0.049 s |
| `NFFT_MEASURE` (default) | **49.1 s** | 0.015 s | 0.022 s |
| `NFFT_NO_DIRECT` | 0.001 s | 0.022 s | 0.026 s |

Analytic `pcost` of the direct candidate is ~4.3e9 vs ~5e7 for the fast one,
a factor ~100, and the planner still runs one full direct apply. The
per-candidate budget `PLNR_TIME_LIMIT_SECONDS = 2.0`
(`kernel/planner/timer.c:87,142`) is checked only after a batch completes,
so one slow apply is never cut short. FFTW does not have this problem only
because all of its candidates are O(n log n); here the candidate set mixes
complexity classes, so the estimate must gate the race.

Fix (1 is the minimum):
1. Prune before timing: skip any candidate whose `pcost` exceeds
   `K * min_pcost` (K ~ 4-10) unless an exhaustive flag is set. This is the
   use the unused `PLNR_BELIEVE_PCOST` bit was meant for.
2. Time the cheapest-by-estimate candidate first and give each later
   candidate a budget of a few multiples of the current best.
3. Until fixed, make `NFFT_ESTIMATE` the default and document `MEASURE` as
   opt-in.

Add a regression test: planning time for a mid-size 2D problem must stay
below a small multiple of one fast execute.

### 1.4 Unmangled global symbol `hash` exported from the library

**Status 2026-08-22: fixed.** `static unsigned reg_nam_hash()`. `nm -D` no longer lists `hash`. The CI symbol check is not added yet.

`kernel/planner/planner.c:270` `unsigned hash(const char *s)`: external
linkage, no `Y()`. `nm -D .libs/libnfft3.so` lists `T hash`. Breaks the
mangling rule in CONVENTIONS.md, collides with any application symbol named
`hash`, and with itself when `libnfft3` and `libnfft3f` are both linked.
`kernel/nfft/problem.c:47` has a `static void hash(...)` of the same name.
Fix: `static unsigned reg_nam_hash(...)`. Add a CI check that
`nm -D --defined-only` on new objects yields only `nfft[fl]?_`-prefixed
symbols.

### 1.5 Incompatible-pointer warnings in new code

**Status 2026-08-22: fixed.** `(FC *)` casts in `nfft-nd.c`; zero `-Wincompatible-pointer-types` in Autotools and CMake builds.

`kernel/nfft/nfft-nd.c:220-221` passes `C *` to `FFTW(plan_dft)`, which
takes `fftw_complex *`: four `-Wincompatible-pointer-types` warnings in
double and float builds. gcc-14 makes this an error under `-Werror` and
several distributions make it an error by default. `infft.h` introduces `FC`
for this purpose and then does not use it here. Fix: cast via `(FC *)`.

### 1.6 Documentation cites records that are not in the repository

**Status 2026-08-22: fixed by stripping.** Removed the `docs/adr/0005` header line from 14 `kernel/planner/*.c` files; rewrote the affected sentences in `CONTEXT.md` (also corrected `plan_ng.c`/`api_ng.c` file names, "aliased not copied", and the two-plan claim in the same sentences), `AGENTS.md`, `tests/nplan.c`, and the planner skill (provenance section reworded, ADR quick-map removed). Only ADR-0002/0004 references remain; both exist. No ADRs were written: their content is not in the tree to recover.

All 14 `kernel/planner/*.c` headers cite
`docs/adr/0005-clean-room-planner-infrastructure.md`. `docs/adr/` contains
0001-0004. `CONTEXT.md`, `AGENTS.md:90`, and the skill docs cite
ADR-0005..0018, "spec 05/06/07/08", decisions "D8/D11/D12/D13",
"Phase 2/2.5/3/4", commit `408f1bee`, and
`.scratch/planner-v2/specs/00-clean-room-protocol.md`; `git ls-files
.scratch` is empty. A reader of `develop` after merge can resolve none of
them. Fix: commit the ADRs (at least the provenance one, see 2.5) or strip
the references and fold surviving content into `CONTEXT.md`.

---

## 2. Architecture and design

### 2.1 Dormant two-direction scaffolding in the bundle (major)

`kernel/nfft/plan.c` models the bundle as `problem *prob[2]; plan *dir[2]`
with `FWD`/`ADJ`, loops `for (dirn = FWD; dirn <= ADJ; ...)`, five negative
sentinels in `ncands[]`, a `DIR_IS_ABSENT` macro, and `select_estimate()`
that plans both directions, while `prob[ADJ]` is always `NULL` by
construction ("dormant for a potential option with two plans"). The plan
tree printed to users ends in `(adj (null))`. The guru is ~350 lines; the
reachable logic is about a third.

I would delete the second direction: one problem, one plan,
`apply`/`apply_adjoint`. An adjoint-specialised race later is a small
additive change. Speculative generality has already produced wrong
documentation (CONTEXT.md: "one bundle owns two per-direction problems and
two independently planner-selected plans").

### 2.2 Guru error handling by repetition (major)

`Y(plan_ng_guru)` restores `pl->flags` and destroys `prob[]`/`p` in four
separately written paths; the `need_restart` block re-implements the
estimate setup inline. One `fail:` label, or split into `build_problem`,
`select_estimate`, `select_measured`, each with one exit.

### 2.3 Global mutable planner state as a calling convention (major)

Impatience bounds reach solvers by mutating `the_planner()->flags` before a
call and restoring after. `mkplan_native_fast` (`nfft-nd.c:232,237`)
ignores its `pl` parameter and calls `Y(the_planner)()` again. Not
re-entrant, not thread-safe, and save/restore order is part of correctness.
Pass the bounds explicitly (`planner_mkplan(pl, p, &flags)`), keep
`planner_s` free of per-call state, use the `pl` that is already threaded
through. Then either add `nfft_make_planner_thread_safe()` (FFTW parity) or
state in `nfft3.h` that the guru must be serialised; today nothing public
says so (`theplanner.c` says it in a comment).

### 2.4 `mkplan` violates its own contract (major)

`iplanner.h:170` says `mkplan` "must be cheap and must not build
node-dependent state". `mkplan_native_fast` allocates two `ntot`-sized
scratch buffers and creates two FFTW plans with the user's `fftw_flags`
(seconds with `FFTW_MEASURE`); every CONV solver allocates `psi`
(`M*d*(2m+2)` reals) at mkplan (`conv-1d.c:220`, `conv-2d.c:287`,
`conv-3d.c:584`, `conv-nd.c:197`). Estimate mode builds every candidate
and destroys all but one; measured mode builds all before timing any; a
wisdom hit re-runs the winner's mkplan. Move buffers, FFTW plans and `psi`
into `awake()`; `mkplan` becomes an applicability check plus cost model.

### 2.5 Provenance: the "clean-room" claim does not match the artifact (major)

Compared against the FFTW 3 source in `.fftw3-reference/` (git-excluded,
present on this host). Facts:
- `struct printer_s` and `struct scanner_s`: same members in the same
  order as FFTW `ifftw.h` with two fields renamed.
- `md5uint`/`md5sig`/`md5` including the `UINT_MAX >= 0xffffffffUL` selector
  and the 32-bit-only length word in `md5_end`.
- `flags_t`: same five bitfields, same order, same names modulo
  `hash_info`->`info`, `timelimit_impatience`->`timelimit_imp`.
- `PLNR_BELIEVE_PCOST = 1, PLNR_ESTIMATE = 2`, `BLESSING/H_VALID/H_LIVE =
  1/2/4` with FFTW's comments; `LEQ`, `INFEASIBLE_SLVNDX`, `slvdesc`,
  `solvtab`/`SOLVTAB`/`SOLVTAB_END`, `FORALL_SOLVERS`, `htab_blessed/
  unblessed`, `cur_reg_nam/cur_reg_id`, `nslvdesc/slvdescsiz`, `nthr`.
- Wisdom grammar: export line `(%s %d #x%x #x%x #x%x #x%w #x%w #x%w #x%w)`
  equals FFTW's with `%M`->`%w`; import loop, 64-byte name buffer,
  snapshot-and-rollback; `config_signature` = FFTW's
  `signature_of_configuration`.
- Hash table: probe `s[0] % size`, step `1 + s[1] % (size-1)`, load bound
  `nelem + nelem/8`, growth to `next_prime` of the bound squared.
- Timer: doubling `n`, best-of-8 repeats, floors 100 ticks / 1e-3 s, 2.0 s
  limit, `-1.0` sentinel.
- Genuine differences: djb2 string hash; `tensor`/`mvdim` with rectangular
  factors is a real generalisation; `%w` not `%M`; `amnesia` enumerators.

Both projects are GPL-2.0-or-later, so distribution is lawful. But "no FFTW
source (verbatim or adapted)" is not a defensible description of this code,
and 14 file headers assert it while pointing at a document that is not in
the tree. I would: remove the word "clean-room" from every file header; add
one provenance statement in `CONTEXT.md` and an ADR ("architecture and
planner kernel derived from FFTW 3.x (Frigo/Johnson, MIT), GPL-2.0+;
`include/cycle.h` is FFTW's MIT file"); add FFTW's copyright notice to the
derived files as GPL section 2 requires; delete the flags that exist only
because FFTW has them (2.6).

### 2.6 Over-engineering for a five-solver roster (major)

NFFT has 5 solvers; DECONV and CONV one per rank. For that roster the
branch carries: a growable prime-sized double-hashed wisdom table with
kill/rehash and a subsumption lattice over 16-bit bounds; a per-kind linked
list through a growable descriptor array; a refcounted solver base; a
registrar-name hash to avoid `strcmp` on import; generation counters for
malloc ABA; a candidate cap of 8 guarded by `A(0)`. Defined and never used
anywhere in `kernel/`: `PLNR_BELIEVE_PCOST`, `PLNR_NO_UGLY`,
`PLNR_ALLOW_PRUNING`, `FORALL_SOLVERS`, `PLNR_FORGET_UNBLESSED`
(tests only), `PLNR_NO_SLOW` (tests only, as an arbitrary bit), the 12-bit
`timelimit_imp` (only ever written as 0), `plan.uses_x` (written twice,
read nowhere), six of seven `Y(tensor_*)` functions (`copy`, `equal`,
`squarep`, `kosherp`, `compress`, `compress_contiguous`: tests only; keys
hash the raw tensor, so `compress` never reaches wisdom).

None of this is wrong; all of it is cost with no present benefit, and the
lattice (estimate in `u` only, `LEQ` subsumption, mutual-subsumption
healing) has semantics that only the tests explain. I would keep the store
and registry but implement them as the simplest correct thing (exact-match
lookup on md5, flags as one `unsigned`) and grow into a lattice when a
second impatience level exists. If the lattice stays, write its rules once
in `iplanner.h` and delete every unused name.

### 2.7 OpenMP: the new fast path is serial in `libnfft3_omp` (major)

The legacy path calls `FFTW(plan_with_nthreads)` (`nfft.c:6026`); the
native fast solver never does (`nfft-nd.c:220-221`), and no new file has a
`#pragma omp`. Yet `plan.c:98` keys `nthr = get_num_threads()` into wisdom,
and `kernel/Makefile.am:99` links serial `libconv.la`/`libdeconv.la` into
`libkernel_threads.la` while `planner/` gets a `_threads` twin compiled with
`OPENMP_CFLAGS`. Either thread the FFT child and the CONV/DECONV loops
(then the key is honest) or remove `nthr` from the key and the `_threads`
twins, and say in the docs that the new API is serial.

### 2.8 Timer clock (major)

`timer.c:48-50` routes the budget and the slow-path fine clock through
`Y(clock_gettime_seconds)`, which is `CLOCK_REALTIME`
(`kernel/util/time.c:36`). An NTP step or suspend moves the 2 s budget, and
on the `HAVE_CLOCK_GETTIME` arm a negative `batch_sec` becomes `batch_min`
(`timer.c:139-140`). The tick arm does not reject a negative `elapsed()`
either (`timer.c:84-85`). Use a planner-private `CLOCK_MONOTONIC` reader
and treat non-positive readings as "no reading". `timer.c` is also two
45-line copies differing only in clock type; one `fine_now()` and one floor
constant remove the duplication, and the comment justifying the split
("getticks() to 0U compiles only as a whole statement") is stale once
`ticks.h.in` drops the semicolon, which this branch does.

### 2.9 Cost units mixed in one tree (minor)

After a measured race `pln->pcost` is overwritten with ticks; children
keep analytic costs. `nfft_fprint_plan` prints
`(nfft_solver_fast_native pcost=7 (deconv ... pcost=64) ... (conv ... pcost=800))`
(verified). `iplanner.h:140` says "measured seconds"; the tick arm returns
raw ticks and `TICKS_PER_SECOND` is unused. Keep `pcost` analytic, store the
measurement in a separate field, print both with units.

### 2.10 Runtime `window` ordinal vs compile-time window (minor)

The library is compiled for one window (`--with-window`), yet the guru
takes `int window`, solvers dispatch through `Y(window_phi)` /
`Y(window_phi_hut)` switch tables that duplicate the formulas of
`infft.h:175-238`, and `nfft_get_window_id()` exists only to pass the
compile-time choice back in. The runtime design is sound (only `phi`,
`phi_hut` and `b = pi(2 - N/n)` depend on the window on the planner path),
but it is half done: `Y(get_default_window_cut_off)` is compile-time, so a
non-default runtime window has no default `m`; unknown or Dirac ordinals
return 0 silently from `window.c:219,247` and enter the wisdom key
unguarded (only `nfft-nd.c:173` checks the range; the NDFT solvers accept
anything). Either finish it (guard in `mkproblem_nfft`, per-window default
`m`, re-express the legacy macros through the new functions) or drop the
parameter.

### 2.11 No `nfft_cleanup` (minor)

`the_planner_destroy()` is internal only. FFTW has `fftw_cleanup()`;
valgrind shows 54 KB "still reachable" at exit. Export `nfft_cleanup()`
and document that plans must be destroyed first.

---

## 3. Public API (`include/nfft3.h:897-962`)

- **(major)** `nfft_plan_ng` in a permanent public header. `_ng` is a
  branch label, not a name. CONTEXT.md says `_ng` marks "colliding nouns
  only", which concedes the point. Pick a noun that ages: `nfft_plan2`,
  `nfft_tplan`, or a namespaced prefix. Decide before merge; renaming after
  release is an ABI break.
- **(major)** `int X(get_window_id)();` (`:958`): K&R declarator and no
  `NFFT_EXTERN`. Under `NFFT_DLL` this is an unexported symbol on Windows.
  Same pre-existing defect at `:889-890`.
- **(major)** Doc comments say nothing about the contract the skill calls
  "the #1 thing to get right": mandatory `precompute`, `x` copied,
  `f_hat`/`f` borrowed and clobbered by a measured race, no thread safety
  of the guru. These belong on the prototypes, not in a skill file.
- **(minor)** `//` comments in a header that is `/* */` elsewhere
  (`-std=c89 -pedantic` consumers break).
- **(minor)** `set_timelimit` comment: "applies to all kinds
  (NFFT/NFCT/NFST)". NFCT/NFST have no planner API.
- **(minor)** Flag names leak implementation: `NFFT_NO_FAST_NATIVE`,
  `NFFT_NO_NDFT_BLOCKED`, `NFFT_NO_NDFT_PLAIN`. "native" is meaningful only
  relative to a removed design. FFTW names user intent. Suggest
  `NFFT_NO_FAST`, `NFFT_NO_DIRECT`; keep per-algorithm bits internal.
- **(minor)** `NFFT_WINDOW_*` ordinals duplicate the configure-time window
  selection with a second spelling; `NFFT_WINDOW_DIRAC_DELTA` is public
  and "unsupported here".
- **(minor)** `include/plan_ng_test.h` / `nfft_plan_ng_test_awake_state`
  compile a test hook into the shipped library (`nm`: exported). Nothing
  test-only belongs in `libnfft3.so`; use a test-build define or link the
  test against the object.
- **(nit)** `int d, int m, int window` next to `NFFT_INT M`; say why.
- **(nit)** Comment register varies: "Create transform plan." / "printplan
  tree" / "build the node-dependent tables, call after filling x" (x is
  read at construction).

---

## 4. Planner core (`kernel/planner/`, `include/iplanner.h`)

- **(major)** Registrar names persisted in wisdom are double-mangled.
  `SOLVTAB(Y(nfft_solver_fast_native_register))` stringizes after
  expansion, and the inner `nfft_` is part of the identifier, so the wisdom
  file contains `nfft_nfft_solver_fast_native_register` (verified by export;
  `nfftf_nfft_...` in float). The same doubling affects every exported
  symbol declared as `Y(nfft_...)` in `iplanner.h:428-450`
  (`nfft_nfft_ensure_registered`, `nfft_nfft_x_md5`, ...). Fix before any
  wisdom file exists: drop the inner `nfft_` (`Y(solver_fast_native_register)`
  after checking no clash with `kernel/solver`), and
  `#define SOLVTAB(s) { Y(s), #s }` with bare names.
- **(major)** `planner_candidates` drops candidates beyond `cap` in release
  builds (`A(0)` is a no-op). Return the needed count or allocate.
- **(minor)** Import masks `l`, `u` to 16 bits and `tli` to 12
  (`planner.c:389-391`) instead of rejecting out-of-range; `read_hex`
  (`scan.c:56-78`) accepts unbounded digits and masks. Reject.
- **(minor)** Wisdom key feeds raw `int`/`INT` bytes (`md5.c:195-203`) so
  keys differ by endianness and `sizeof(INT)`, but `config_signature`
  covers only `sizeof(R)`. Add `sizeof(INT)` and an endian word.
- **(minor)** A blessed entry whose solver now returns `NULL` is memoised
  unblessed on the miss path and keeps winning the blessed-first lookup
  (`planner.c:553-564`): one wasted `mkplan` per call forever. Kill it.
- **(minor)** `nthr` is in the core key but only the NFFT guru sets it.
  Own it in the planner or let the problem hash it.
- **(minor)** `printer_create_str` is unbounded by design; tests use fixed
  64-512-byte stack buffers with comments admitting the risk. Add
  `printer_create_strn(buf, cap)`.
- **(minor)** `theplanner.c` documents the thread-safety contract in a
  comment; it is not in the public header (see 2.3).
- **(nit)** `planner.c:235` example says `nfft_wisdom-3.5.4alpha`; the
  version is `3.6.0alpha`. `iplanner.h` has 14 `//` lines. `flags_t`
  "packs into 64 bits" is a comment, not a guarantee. Typedef names
  (`plan`, `problem`, `solver`, `planner`, `printer`, `scanner`, `md5`,
  `tensor`) are unmangled and identical to FFTW's; `iplanner.h` is never
  installed, so the only hazard is a TU that includes FFTW's `ifftw.h`
  alongside; note it in the header.

---

## 5. NFFT bundle and solvers (`kernel/nfft/`)

- **(major)** Docs and the code disagree on file names. The skill,
  `plan_ng_test.h:21`, `xcheck.h:21`, `AGENTS.md:90`, `CONTEXT.md:159,167`
  reference `kernel/nfft/plan_ng.c`, `api_ng.c`, `nfast.c`, `ndft.c`,
  `nconst.c`; the files are `plan.c`, `api.c`, `nfft-nd.c`, `ndft-1d.c`,
  `rnk0.c`. Code comments too: `ndft.h:20` ("ndft.c"), `conf.c:51`,
  `deconv-2d.c:195`, `deconv-3d.c:241`, `tests/nplan.c:300,2629`,
  `tests/nfast.c:1127` ("nfast.c"). A rename happened after the docs were
  written and nobody re-grepped.
- **(major)** `kernel/nfft/xcheck.h` is empty (guards and a comment); its
  declarations live in `iplanner.h`. The `NFFT_DEBUG` md5 guard over `x` in
  the race is, by its own comment, "a tautology — intentional substrate for
  a future permuting solver". Delete the header, the include, and the
  speculative guard with `xcheck.c`; `AGENTS.md:87-93` then needs rewriting.
- **(major)** Four names for one solver: file `rnk0.c`, registrar
  `nfft_solver_rnk0_register`, printed `nfft_solver_const_0d`, wisdom
  `nfft_nfft_solver_rnk0_register`; tests use the third. The NDFT plans
  print registrar-name-minus-`_register`; make rnk0 follow that rule.
- **(minor)** `nfft-nd.c:60` `// TODO: Defer to cost estimates for DECONV
  and CONV`; `problem.c:82` `TODO: hash the pointers' alignment class`.
  Resolve or file; no TODOs in merged code.
- **(minor)** `conf.c:44-66` explains a use-after-free that lazy
  registration inside the recursion would cause. Root cause is
  `FORALL_SOLVERS_OF_KIND` holding a raw pointer into a growable array;
  close registration before any planning and the comment goes away.
- **(minor)** `api.c`: `fprint_plan` calls `nfft_ensure_registered()` and
  only prints. `export_wisdom_to_filename` does count + string + `fwrite`
  although `printer_create_file` exists.
- **(minor)** `keyable_fftw_flags` strips `DESTROY/PRESERVE_INPUT` in
  `plan.c` while `nfft-nd.c:226` re-normalises them when planning the FFT.
  One place: the problem constructor.
- **(minor)** `rnk0.c:77` `pcost = 2.0*M` but pcost is per apply (one pass);
  `uses_x` left at 1 though rnk0 never reads `x` (dead field anyway).
- **(minor)** `ndft-1d.c:26-40` header comment is copied from `ndft-nd.c`
  and describes an odometer carry that does not exist in 1D; two duplicate
  "Block size" comments. The plain 1D variant uses an unreduced phase
  `2 pi k x` (error ~eps*pi*N at large N); say so, since the blocked
  variant reduces it.
- **(nit)** `ndft.h` exists to share two functions between two files
  under 260 lines; one `ndft.c` removes the header and two exports. The
  `#undef X / #define X(name) NFFT(name)` boilerplate is copied into
  `api.c` and `plan.c` while the solver files avoid `X()` entirely; pick
  one. `ndft-nd.c:apply` is a near-copy of legacy `X(trafo_direct)`
  (`nfft.c:145-210`); acceptable by design, but say that rather than "two
  changes compared to legacy".

---

## 6. DECONV and CONV (`kernel/deconv/`, `kernel/conv/`)

- **(major, latent wrong result) Fixed (`0f1bcb00`).** `deconv-1d.c:74`
  zero-padded from `N/2`, but the last non-negative frequency lands at
  `N-1-N/2`; for odd N the cell `n - N/2 - 1` was never written. Fixed to
  `memset(g + Npos, ...)` with `Npos` hoisted into the plan, plus a
  `return 0` decline for `n < N`.
- **(major, latent wrong result) Fixed (`eb9d0851`, `d3c679cd`,
  `024f4fa3`).** `deconv-2d.c`, `deconv-3d.c`, `deconv-nd.c` ignored
  `variant` (type-II) entirely and did not decline; only `deconv-1d.c`
  applied the shift. A type-II rank-2 DECONV problem hashed to its own key
  and got a type-I plan. All three now carry the per-axis `Nneg`/`Npos`
  slot split, same as `deconv-1d.c`.
- **(major)** `conv-2d.c:30-34` `Y(conv_b_pcost)` hard-codes rank 2
  (`2*M*(2m+2)^2`) and is called for rank 3 (`conv-3d.c:586`) and rank >= 4
  (`conv-nd.c:200`); `nfft-nd.c:66` models the same step as `(2m+2)^d`.
  Harmless today (one solver per rank) and wrong the day the `TODO` in
  `nfft-nd.c:60` is done. Use `pow(2m+2, rnk)`, define it once in `conv.c`;
  same for `deconv_d_pcost` (`deconv-2d.c:30` -> `deconv.c`).
- **(major)** Duplication with no measured justification. Forward and
  adjoint bodies are mirror copies differing in the innermost assignment:
  `conv-3d.c` 73-285 vs 288-499 (212 lines twice), `conv-2d.c` 77-143 vs
  146-210, `deconv-2d.c` 67-112 vs 115-156, `deconv-3d.c` 69-135 vs
  138-200; `deconv-nd.c:63` and `conv-nd.c:72` already use one
  `run(pln, p, forward)`. `uo2` is duplicated verbatim (`conv-2d.c:47`,
  `conv-3d.c:41`). `conv-3d.c` spends 425 of 593 lines on 16 hand-unrolled
  wrap branches that a per-axis run list replaces in ~40 lines with the
  same memory order. The window-range check is copy-pasted 8 times.
  Nothing in `docs/`, `benchmarks/`, `CONTEXT.md` compares 2d/3d against
  the nd path. I would ship 1d and nd only and reintroduce 2d/3d on a
  CodSpeed number.
- **(minor)** `conv-1d.c:25-32,137-180`: a dead reference implementation
  behind `#define CONV_OPTIMIZED`; `tests/nfast.c:633-650` already has an
  independent reference. Delete. The 1D optimisation (wrapped start `u[j]`
  at awake) is absent in 2d/3d/nd, which recompute `FLOOR+LRINT+2 %` per
  node per apply (`conv-2d.c:87`, `conv-3d.c:84`, `conv-nd.c:104-110`).
- **(minor)** `deconv-nd.c:86` `memset(f_hat)` in the adjoint although
  every cell is written once for even N; `deconv-1d.c:65-67,86-88`
  re-derive N, n, type_ii from the problem per apply though the plan holds
  them.
- **(minor) Partially fixed (`706c1ad7`).** Release builds had no guards:
  `A(N % 2 == 0)` and the VLA bound in `conv-nd.c:194` were no-ops; direct
  callers with odd N, `n <= 2m+1` (`uo2` wrap needs `m <= n/2`) or `n < N`
  got silent wrong results. The DECONV `mkplan`s now `return 0` for `n < N`
  and the composed fast-NFFT path additionally requires `sigma = n/N > 1`
  strictly (`guards_ok` in `kernel/nfft/nfft-nd.c`) — the oversampling guard
  named here is in place. The `conv-nd.c:194` VLA bound and the CONV-side
  `n <= 2m+1` case are unchanged.
- **(minor)** Naming: 1d files use `deconv_plan`/`conv_plan`, the others
  `deconv_2d_*`; `conv_adjoint_2d_compute_serial` has no omp twin;
  `deconv_d_pcost`/`conv_b_pcost` infixes are opaque. `conv.c:42` hashes an
  `INT` with `md5_put_int` (narrowing); `cv_print`/`dv_print` omit `window`
  (and `N`) though both are hashed, so wisdom dumps cannot tell windows
  apart.
- **(minor)** Stale and provenance comments: "guards_ok (nfast.c)"
  (`deconv-2d.c:195`, `deconv-3d.c:241`); legacy macro names as comments
  (`MACRO_D_init_result_A`, "not ported: the d==4 and d==5 hand-unrolled
  branches", `deconv-nd.c:82-114`, `conv-nd.c:25-129`); "(sparse PRE_PSI
  strategy in legacy code)" in all four `conv-*.c` headers;
  "defense-in-depth"; "nD, n>=4" meaning d.
- **(nit)** `conv.c:2` copyright misspells an author ("Jens Keiter").
  Mixed declarations and statements (`conv-2d.c:59-63`, `conv-3d.c:53-60`).
  Trailing whitespace in all eight headers.

Wisdom consequence of the three fixes above: none of them touch the wisdom
key (`variant[]` and `N` parity were already hashed, and the configuration
signature only covers `sizeof(R)` plus the solver roster, neither of which
changed), so no existing wisdom file is invalidated. What changes is that an
entry blessed before this work, for an odd-`N` or type-II problem, names a
direct NDFT solver the planner would no longer pick — it stays correct, just
slower. A user who wants the fast path for a shape they already measured
must call `nfft_forget_wisdom()` or delete their wisdom file.

Verified OK: every register uses `REGISTER_SOLVER`; all non-static symbols
mangled; no `#pragma omp`; adjoint symmetry and quadrant/octant index maps
for even N; ownership in `destroy` including partial construction.

---

## 7. Window and Bessel (`kernel/util/window.c`, `bessel_i0.c`)

- **(major, by omission)** The Kaiser-Bessel self-normalisation
  (`phi` and `phi_hut` both scaled by `1/I0(mb)`) is justified: the factor
  cancels exactly and `phi_hut^d` at d = 3 leaves float range without it.
  But it changes the numeric contract of `psi` relative to the legacy API
  and is documented only in `window.c:67-75` and `tests/nfast.h:94-106`.
  It needs an ADR and a line in the API docs.
- **(minor)** `kb_phi_eval` (`window.c:92-117`): the two `a > 0` branches
  are algebraically identical; keep the stable one. The `xpk <= SPLIT`
  branch `0.5*(EXP(u-lg) - EXP(-u-lg))` cancels as `u -> 0`; use
  `SINH(u)*EXP(-lg_peak)`.
- **(minor)** `bessel_i0.c:351-460`: four copies of the same range split;
  factor one `i0_parts()`. `:444` "cancellation-free" is false on the
  `x <= SPLIT` branch. `infft.h:1491-1493` "single source of truth (kept in
  sync with bessel_i0.c)" is self-contradictory; the real constraint is
  that P2/Q2 are fitted on (15, inf), so SPLIT is not tunable. No test
  checks continuity of `log`/`scaled`/`logtail` across SPLIT; the tests sit
  in `tests/nfast.c`, not beside `check_bessel_i0`.
- **(minor)** `window_phi_precompute` divides `idx/n` per tap and
  `kb_phi_eval` multiplies by `n` again; pass `n*x - idx`.
- **(nit)** `kb_consts(int window, ...)`: `window` unused (the one new
  `-Wextra` warning). `window.c:84` `&& xpk > SPLIT` is implied.
  `window.c:82` "SQRT is safe" holds only for `|k| <= N/2`, `N <= n`; add
  `A()`. `kernel/nfct/Makefile.am` and `kernel/nfst/Makefile.am` diffs are
  whitespace-only; drop from the branch.

---

## 8. Style, formatting, comments (systemic)

- **(major)** ~1,100 occurrences of
  ```c
  Y(free)
  (ptr);
  ```
  across the branch (386 in `tests/nplan.c`, 238 in `tests/nfast.c`, 187 in
  `tests/planner.c`, 52 in `kernel/planner/`, 58 in conv/deconv, 27 in
  `plan.c`, ...), including under unbraced `for` bodies (`conv-3d.c:61-65`).
  Cause: clang-format treats `Y(free)` as a statement macro and the repo
  `.clang-format` does not declare the mangling macros. Verified with
  clang-format 18.1.3 in the container: adding
  ```yaml
  Macros:
    - 'Y(x)=x'
    - 'X(x)=x'
    - 'NFFT(x)=x'
    - 'FFTW(x)=x'
  ```
  repairs the broken form automatically. One mechanical commit before
  content review.
- **(major)** The new code does not follow the repo `.clang-format` apart
  from the macro issue: function braces are cuddled
  (`static int subsumes(...) {`) where the style says `AfterFunction: true`
  and CONVENTIONS.md:56 says BSD. Dry-run drift: `planner.c` 109 lines,
  `nfft-nd.c` 92, `ndft-1d.c` 67, `conv-3d.c` 62, `tensor.c` 63,
  `window.c` 177, `bessel_i0.c` 648. Either the code or the style file is
  wrong; the repo already said which.
- **(minor)** Comment register. Many comments narrate process rather than
  describe code: "the accepted one-extra build", "Phase 2.5", "for
  FFTW-trained readers", sentinel tables, "not ported". The repo ships a
  `deslop` skill for exactly this; run it on the branch and then decide
  whether the skill itself belongs in the repo (its samples cite
  `nproblem.c`, "Item 15", "Task 9": names from this branch's history).
- **(minor)** `//` comments: 6 in `nfft3.h`, 14 in `iplanner.h`, 6 in
  `bessel_i0.c`; the rest of the tree is `/* */`.
- **(nit)** Trailing whitespace in most new headers; `tests/Makefile.am:67`
  recipe indentation changed from one tab to two.

---

## 9. Build matrix and evidence

| configuration | library | `checkall_ng` | notes |
|---|---|---|---|
| double, `--enable-all --enable-openmp --enable-tests`, gcc-14 `-O3` | ok | 96/96, 1.3 s; `_threads` 96/96 | 4 new `-Wincompatible-pointer-types` in `nfft-nd.c` |
| float, `--enable-float` | ok | 96/96 | same 4 warnings (`fftwf_plan_dft`) |
| double, `--enable-debug` (ASan/UBSan, `-Wall -Wextra -Wconversion`) | ok | 96/96 | no sanitizer reports; new code adds 1 warning (`window.c` unused parameter) |
| CMake, `-DNFFT_ENABLE_TESTS=ON -DNFFT_ENABLE_OPENMP=ON` | **fails to link** (fixed 2026-08-22: builds, ctest 2/2) | n/a | 1.1 |
| `make distdir` | **incomplete** (fixed 2026-08-22) | n/a | 1.2 |
| valgrind `examples/nfft/nfast_native` | 0 errors, 0 definitely lost, 54 KB reachable | | 2.11 |
| `nm -D libnfft3.so` | exports `hash`, `nfft_plan_ng_test_awake_state`, 91 `nfft_`-prefixed internals | | 1.4, 3 |
| long double | not built (too slow on this host; CI covers it) | | |

Pre-existing and out of scope: `make -j` race in `examples/nfsft`
("file too short" on a shared object file), `fscanf` unused-result
warnings, `libnfft_threads_la_LDLAGS` typo (`kernel/nfft/Makefile.am:14`;
the branch edits the adjacent lines, fix it while there), `configure.ac:387`
replacing user `CFLAGS` and adding `-v` under `--enable-debug`.

---

## 10. Docs, examples, repository hygiene

- **(major)** CONTEXT.md contradicts code: `:169` "made `x` ... aliased
  not copied" (code: `plan.c:131` `copy_x=1`); `:72-75` "ctest suite covers
  util / planner / nplan / nfft / nfct / nfst" (planner/nplan are in
  `checkall_ng`, never mentioned, the likely origin of 1.1's second half);
  `:167` "six wisdom/print wrappers in `api_ng.c`" (seven, in `api.c`);
  "two independently planner-selected plans" (one). AGENTS.md:75,119-134,193
  omit `checkall_ng`, `checkall_ng_threads` and the `_ng` XML names.
  README.md says nothing about the new API.
- **(major)** `examples/nfft/nfast_native.c` is a benchmark harness, not a
  showcase: 578 lines, Welford statistics, a 464 KB data file baked in by
  absolute path (`Makefile.am:44`), includes `ticks.h`/`config.h`, and has a
  type bug at `:396,:473`: `(int *)n` where `n` is `NFFT_INT n[1]`
  (`ptrdiff_t`), passed to legacy `init_guru(int *)`; works on little-endian
  only. `ndft_fast.c` is unchanged from `develop` and does not use the new
  API, yet the skill presents it as a planner example. Move `nfast_native.c`
  to `benchmarks/` or `tests/`, fix the cast, and add a ~60-line
  `simple_test_ng.c` (guru -> precompute -> execute -> destroy on random
  data, no file).
- **(major)** `examples/nfft/data/nfft_1d_8192_128.txt` is byte-for-byte
  `tests/data/nfft_1d_8192_128.txt` (`cmp` identical, 8,457 lines each).
  One copy.
- **(minor)** `.gitignore:184-186` ignores two named files under
  `.scratch/planner-v2/` while `.scratch/` itself is neither tracked nor
  ignored; `:181-182` ignore tool-private dirs (`.worktrees`,
  `.superpowers`) that belong in `.git/info/exclude`. Ignore `.scratch/`
  whole or nothing.
- **(minor)** Committed skills: `.claude/skills/understanding-the-planner-api/`
  (1,100 lines) and `deslop/`. Defensible (the repo commits `docs/agents/`),
  but the planner skill cites the renamed files and missing ADRs throughout
  (`history-and-drift.md:3,20,32-43,71-78` is an "ADR status quick-map" of
  documents that do not exist), its "minimal example" calls the internal
  `Y(nfft_ensure_registered)()`, and it says "five-way" where the example
  header says "Seven-way". Skill prose rots faster than code; keep it short
  and generate the file list from the tree.
- **(minor)** `tests/check_ng.c:45-50` re-registers the util suite that
  `checkall` already runs.
- **(nit)** `docs/agents/domain.md:18-19` ADR names do not match
  `docs/adr/` (pre-existing).

---

## 11. Tests (`tests/`)

| file | test functions | lines | API used | under `make check` |
|---|---|---|---|---|
| `tests/planner.c` | 20 | 1,610 | `iplanner.h` only | yes (`checkall_ng`, `_threads`) |
| `tests/nplan.c` | 44 | 2,795 | ~20 `nfft3.h`-only, 24 `iplanner.h` | yes |
| `tests/nplan_data.c` | 1 (100 cases) | 254 | `nfft3.h` (+`the_planner_destroy`) | yes |
| `tests/nplan_perm.c` | test solver | 107 | `iplanner.h` | n/a |
| `tests/nfast.c` | 25 | 2,135 | ~13 `nfft3.h`-only, 12 `iplanner.h` | yes |
| `tests/check_ng.c` | driver | 157 | | yes; CMake: link failure (1.1) |

Of 90 test functions, 46 touch only internals. `checkall_ng`: 96 tests,
2,726 assertions, 0 failures, 1.2 s. Where the tests compare a transform
against an oracle (reference files, direct transform, classic API) they are
good. The rest:

- **(major)** refgen is no longer the single source of truth.
  `tests/refgen/registration.py:70-77 render_extra_dist()` iterates
  `_legacy_grid(module)` only; a regeneration drops 34 native-only files
  (`*_t21*`, odd-N) plus `nfft_1d_8192_128.txt` from
  `tests/data/Makefile.am`, so the committed Makefile.am is hand-edited and
  a regen breaks `make dist` for `checkall_ng`. The generated
  `tests/data/generated/nfft_native_testcases.h` says "GENERATED — do not
  edit" and differs from the generator output on 120 lines (it was
  clang-formatted). No CI step regenerates and diffs. Fix: iterate
  `G.GRIDS[module]` with `variant`; make the generator emit the final form;
  add `python -m tests.refgen.generate && git diff --exit-code tests/data`
  to CI. Update `docs/agents/test-methodology.md` and `tests/refgen/README.md`,
  which mention none of the native header, the `_t21` tag, or `checkall_ng`.
- **(major)** `tests/data/nfft_1d_8192_128.txt` is referenced by no test;
  it exists only as the duplicate of the example's data file (10). Not in
  `grids.py`, so not generatable.
- **(major)** A test that cannot fail: `nfast.c:1051-1053`
  `check_nfast_native_declines_window` plans with
  `NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE` and asserts `NULL`; the window is
  irrelevant. Same at `nfast.c:1110-1113`. Its premise (native fast declines
  a non-compile-time window) is false: `check_nfast_native_window_select`
  shows all four windows plan. Delete.
- **(major)** Vestigial tests for a removed design: `nplan.c:1771-1794`
  `check_nplan_per_plan_core` asserts nothing beyond non-NULL;
  `nplan.c:1796-1818` `check_nplan_core_owns_no_data_arrays` has
  `TODO: assert the core flags` and names `kernel/nfft/nsolver.c`, which does
  not exist; `nplan.c:1905-2002` repeats `check_case_against_direct`.
  `nplan.h:80-88` describes "wrapper plan owns its legacy core"; no such
  thing exists in `kernel/`.
- **(major)** Three divergent tolerance tables with one name: `err_bound`
  in `nplan.c:700-702` (KB `b = 2100`), `nfast.c:94-96` (`b = 3000`),
  `nfft.c:315` (`2100`); Gaussian `a = 0.41` vs `1.0`. One `tests/bounds.c`.
- **(major)** Copy-paste volume: four case readers (`nplan_data.c:22`,
  `nfast.c:749`, `nfast.c:806`, `nfft.c:431`), `rel_max_err` three times,
  wisdom string round-trip helpers twice, the fprint-to-buffer block five
  times in `nfast.c`, the `NFFT(init_guru)` oracle block 13 times, the
  `tol = 1e-10; if (1e4*EPSILON > tol)` idiom 11 times; `nfast.c:1128-1147`
  hand-lists 2D files the generated roster already holds. Helpers in
  `util.c`; iterate `native_testcases[]` filtered by `d`/`variant`.
- **(major)** Implementation-detail assertions that freeze what sections
  2.1 and 2.6 say should change: exact plan print strings (`nplan.c:52,76`,
  `nfast.c:183,261`, `planner.c:1054`); roster size
  `CU_ASSERT_EQUAL(the_planner()->nslvdesc, 13u)` (`nplan.c:496,501`);
  hashtable `nelem`/`nrehash` counters (`planner.c:409-420`,
  `nplan.c:1003-1048,1163-1234,1460`); enum literal values
  (`nplan.c:2233-2238`); reverse registration order (`planner.c:1188,1216,1475`).
  Express wisdom effects through behaviour (re-plan does not re-run
  `mkplan`), not counters.
- **(minor)** Float coverage holes: under `NFFT_SINGLE` the 3D file and
  classic-API comparisons are `(void)err;` (`nfast.c:1379-1383, 1438-1443,
  1500-1504, 1532-1537`), so 3D tests assert only the printed plan. Eight
  absolute `K(1e-12)` tolerances in the DECONV/CONV solver tests sit below
  float eps and pass only because solver and test share one evaluator.
  Scale by `EPSILON`.
- **(minor)** `CK()` used as a test assertion (`planner.c:1596-1602`,
  `nplan.c:1515,1545,1596`) aborts the binary instead of recording a CUnit
  failure. 46 `the_planner_destroy()` calls in `nplan.c`, 16 in `nfast.c`;
  a FATAL exit leaves `set_timelimit(0.0)` (`nplan.c:1486`) or the perm
  solver (`:2373`) in place for later tests; use suite setup/cleanup.
  `nplan_data.c:224-237` runs 100 cases in one CU test, so a failure does
  not name the case. `nfast.c:749-864` readers leak on partial read.
  `nplan.c:1899` asserts only `err_blocked < err_plain`, no absolute bound;
  a tie is possible in long double.
- **(minor)** CI never shows `checkall_ng` output:
  `.github/workflows/build-linux.yml:175-211` prints and uploads
  `checkall.log`/`checkall_threads.log` only; `.github` is untouched.
- **(minor)** Stale and task-numbered comments: `nplan.h:163` "completing
  the P1.3 NULL/d checks"; `nplan.h:41,47,51` "legacy direct wrapper
  retired"; `nfast.h:40` (false premise, see above); `nplan.c:927-931,
  1026, 1213, 1798` (union core, `build_core`, `nsolver.c`). `nplan.h:106-110`
  hard-codes "66 reused + 34 new". "legacy" appears 68 times for the
  still-shipping API; "classic" reads better.
- **(minor)** Formatting: 781 `Y(x)\n(args)` artifacts in the tests; 674
  changed lines in `nplan.c`, 527 in `nfast.c` under the repo
  `.clang-format`; 131/50/28/42 lines over 79 columns. See 8.
- **(nit)** `tests/check.c` diff is a whitespace-only comment edit;
  `tests/data/Makefile.am` is a full rewrite from single- to double-space
  continuations; `planner.c:121` determinism check compares md5 words 0
  and 3 only. `check_ng.c:45-50` re-runs the util suite.

What I would do: keep `nplan_data.c` (file oracle) and the oracle-based
tests in `nplan.c`/`nfast.c` as the public acceptance suite; move every
`iplanner.h` test into a separate binary that is allowed to change with the
internals; delete the tautological and vestigial tests; one bounds table,
one reader, one plan-to-string helper.

---

## 12. Suggested order of work

1. Mechanical, one commit each: `.clang-format` `Macros:` + reformat;
   `hash` -> static; `(FC *)` casts in `nfft-nd.c`; CMake source lists;
   `ndft.h`/`xcheck.h` in `_SOURCES` or deleted; `SOLVTAB` stringizing and
   the inner `nfft_` prefix; delete unused `PLNR_*` flags, `FORALL_SOLVERS`,
   `uses_x`, unused `tensor_*`, `CONV_OPTIMIZED` branch, empty `xcheck.h`,
   whitespace-only Makefile diffs, duplicate data file.
2. Default-mode stall (1.3): estimate-gated race + regression test.
   Monotonic clock and non-positive-reading rejection in `timer.c`.
3. Simplify the bundle (2.1, 2.2); explicit bounds instead of global
   mutation (2.3); move buffers, FFTW plans and `psi` from `mkplan` to
   `awake` (2.4). Expect test deletions, not additions.
4. DECONV/CONV: decline odd N and type-II where unimplemented; rank-correct
   `conv_pcost`; fold forward/adjoint mirrors; decide 2d/3d on a benchmark.
5. Public names (3), window parameter (2.10), `nfft_cleanup` (2.11),
   OpenMP story (2.7).
6. Tests: refgen as single source of truth with a CI drift check; delete
   the tautological and vestigial tests; one bounds table, one reader; split
   internals tests from the public acceptance suite (11).
7. Docs: ADRs committed or references stripped (1.6); provenance statement
   (2.5); re-grep every file path in CONTEXT.md, AGENTS.md, skills; contract
   comments on the public prototypes.
8. Second review of the resulting diff.
