# ADR-0006: Operation-count cost estimation in `plan_ng`

## Status
Accepted (2026-09-12).

## Context
Every plan carried one cost scalar, `double pcost`, and it served two masters.
Each solver's `mkplan` wrote an analytic estimate into it; the measured race
then overwrote it with seconds read off the clock. Three defects followed
from that one field:

- **Estimates could not compose.** The fast solver's cost function carried a
  `// TODO: Defer to cost estimates for DECONV and CONV` — it could not add
  its children's costs, because a child's scalar described a plan that was
  never measured and the parent had no way to combine such scalars into
  anything meaningful.
- **One field, two units.** The old field comment read "analytic estimate,
  may get overwritten with measured seconds". Nothing in the type said which
  unit a given `pcost` held; only the call path knew.
- **Machine structure had nowhere to live.** A single double cannot say
  "eight multiplies, six additions, two hundred loads". The old models fitted
  that structure into fudge constants instead, and the fudge was per solver.

FFTW had already solved this, and its solution is what we adopted.

## Decision
Adopt FFTW's operation record, central estimator and evaluation split.

- **`opcnt {add, mul, fma, other}`** (`include/iplanner.h`, after FFTW's
  `opcnt`, `kernel/ifftw.h:350`), doubles rather than `INT` because a direct
  NDFT's count exceeds 2^31 at moderate size. `other` is loads, stores and
  index arithmetic — the column no previous model counted. Seven helpers
  (`Y(ops_zero)`, `Y(ops_other)`, `Y(ops_cpy)`, `Y(ops_add)`, `Y(ops_add2)`,
  `Y(ops_madd)`, `Y(ops_madd2)`) live in `kernel/planner/ops.c`, after FFTW's
  `kernel/ops.c`.
- **Solvers fill `ops` at `mkplan`; the parent composes.** A composite plan's
  count is the sum of its children's, exactly as FFTW's
  `dft/ct.c:222` adds its two halves. The fast solver's count is its DECONV
  child plus its CONV child plus `FFTW(flops)(pln->pfwd, ...)` — the FFT term
  comes from FFTW's own flops counter on the plan it actually built, not from
  a formula. `other` is the two children's alone, since `FFTW(flops)` reports
  none.
- **One central estimator.** `Y(iestimate_cost)`
  (`kernel/planner/planner.c`, after FFTW's `planner.c:426`) collapses the
  record to one comparable scalar with FFTW's weights unchanged:
  `add + mul + fma + other`, where `fma` counts once on hardware that fuses
  and twice on hardware that does not. `HAVE_FMA` comes from C99's
  `FP_FAST_FMA` (`include/iplanner.h`) — that macro already means "the
  hardware fuses", which is the distinction FFTW's own macro draws, so no
  configure probe was added.
- **`pcost` is zero until evaluated** (FFTW's `plan.c:35`). `mkplan` leaves
  it at 0; `Y(evaluate_plan)` writes it — the estimate in estimate mode, the
  measured cost otherwise, falling back to the estimate per candidate when
  the clock cannot resolve the plan (FFTW's `planner.c:446` condition and
  `planner.c:462` fallback). Selection never reads `pcost`:
  `Y(planner_mkplan)` and the race's prune gate compare
  `Y(iestimate_cost)` values only. This gives `PLNR_BELIEVE_PCOST` — set by
  the mapping stage below `NFFT_PATIENT`, and recorded by ADR-0005 as read
  by nothing — its first consumer: with the bit set, `Y(evaluate_plan)`
  skips re-evaluating a plan that already carries a cost. It only suppresses
  re-evaluation; it never changes what a fresh plan costs to learn.

### The operation counts, as the kernels actually run

| solver | add | mul | other |
|---|---|---|---|
| NDFT (rank 1 and nd) | `6*taps` | `8*taps` | `2*taps + 3*taps/B`, `taps = Ntot*M`, `B = NDFT_RECURRENCE_BLOCK` |
| rank-0 solver | 0 | 0 | `M` |
| DECONV rank 1 | 0 | `2*Ntot` | `(n-N) + 2*Ntot` |
| DECONV ranks 2, 3, nd | 0 | `3*Ntot` | `ntot + 2*Ntot + 2*(d-1)*rows` |
| CONV (all ranks) | `2*taps` | `2*taps` | `2*taps + 2*M*d`, `taps = M*(2m+2)^d` |
| fast native | children + `FFTW(flops)` | children + `FFTW(flops)` | children's alone |

### Some terms are estimators, not traces
Three `other` terms are approximations, each commented at its fill site in the
kernel. DECONV's row-setup term undercounts the 2D/3D box-corner code, whose
row body re-executes once per last-axis corner, and overcounts the nd flat
loop, whose incremental multi-index amortises the same work. CONV at rank >= 2
omits a strictly lower-order axis-combine term above the leaf accumulate. The
NDFT's trig term amortises the block-start phase setup over
`NDFT_RECURRENCE_BLOCK` taps. The counts exist to rank candidates, not to
profile cycles; the tests pin them against the kernel's loops so a change to
a loop forces a change to the count.

### The estimate fallback, amended from FFTW
FFTW's contract is that a negative timing means the machine has no cycle
counter — a property of the machine, so a per-candidate fallback is safe
there. `Y(plan_measure_cost)` (`kernel/planner/timer.c`) returns negative for
four reasons and only one, the missing clock at `timer.c:151`, is
machine-wide; the other three are properties of this plan's `apply` against
the clock resolution — a zero-reading level, a frozen clock, or the
per-candidate budget expiring before any strictly-positive batch reading was
taken. Left as a pure per-candidate fallback, a race could
compare an op count against seconds and discard the candidate whose `apply`
was too fast to resolve, and a clockless machine would bless a selection
backed by zero measurement.

So `Y(evaluate_plan)` returns a grade — 1 measured, 0 fell back to the
estimate — and a direction in which any candidate fell back is graded
estimate: the race stops, selects the cheapest candidate by estimate, and
memoises the winner **unblessed**, which is the branch the timeout path
already took. No teardown and no replan. `pcost` on that path carries the
estimate, so `Y(fprint_plan)`'s `pcost=` stays meaningful and the calibration
sweep can still read it.

### The untimed adopts carry the estimate
Three paths hand a plan back without measuring it: the wisdom hit, the lone
survivor of the prune gate, and the timed-out or graded adopt. All three now
write `Y(iestimate_cost)` into `pcost` directly rather than leaving it zero.
`Y(evaluate_plan)` is deliberately not used there: in measured mode it
measures, and those paths exist precisely to spend no measurement —
evaluating a wisdom hit would make every cached plan cost one full timing.

### The deleted restart made the directions independent
`need_restart` and the shared estimate-grade replan are gone. Previously a
bail in the forward race tore down the whole bundle, including a forward
winner that had already been blessed; now each direction resolves on its own.
This retires a property nothing ever wrote down as a contract: that both
directions carry the same evidence grade.

### `NDFT_COST_UNIT` and the fitted constants are gone
`NDFT_COST_UNIT = 2.31` was a fitted scale correction over a model that
counted no memory traffic, living to pull the direct NDFT's fitted cost onto
the fast solver's fitted scale. The `other` column replaces it by counting
the traffic itself. It went together with `NDFT_COST_TRIG`, `NDFT_COST_MUL`,
the three per-solver cost functions and `nfft-nd.c`'s own static `pcost()`.
The calibration re-run (below) found the scale gap the constant existed to
correct shrank from 2.31x to 1.43x but did not vanish — the residual is
reported honestly in `docs/agents/pcost-calibration.md`; the deleted constant
is neither simply vindicated nor simply refuted by it.

### `PLNR_PRUNE_RATIO` recalibrated
The sweep in `docs/agents/pcost-calibration.md` re-ran against
`Y(iestimate_cost)` — 976 geometries over ranks 1-3, the same grid as before.
The ratio is confirmed at 4.0, not carried over: zero harmful prunes across
the grid, a 2.5x margin over the observed unsafe floor of 1.6, and 316 ms of
wasted timing summed over all 976 geometries, which is not worth tightening
against.

## Divergences from the plan's text found during implementation
- **Two rows of the plan's op-count table were wrong, and the kernel won.**
  The NDFT is 8 `mul` and 6 `add` per tap, not 4 and 4: the `w *= dw`
  recurrence step is a second complex multiply the plan did not count, and it
  runs on every tap. The rank-0 solver's forward is a pure broadcast — `add`
  0, `other` `M` — where the plan had `2*M` and `2*M`, numbers that bundled
  in the adjoint's reduce; `ops` records the forward, which is what
  `Y(plan_measure_cost)` times.
- **The estimate terms marked as estimators, not silently exact** — see the
  section above.
- **The FFTW fallback contract was not safe as specified** — see the amended
  fallback section. The plan's deletions also left one live path uncovered:
  a budget expiring before the first candidate is measured leaves no winner,
  so the prune gate's argmin is adopted untimed and memoised unblessed.
- **The test hook is one accessor**, `Y(plan_ng_test_native_parts)`
  (`include/plan_ng_test.h`), not the four the plan sketched. The opaque
  structs it must cross (`plan_ng_s` in `plan.c`, `native_fast_plan` in
  `nfft-nd.c`, which stays file-static) need exactly one crossing per file:
  `plan.c`'s accessor delegates to `nfft-nd.c`'s
  `Y(nfft_fast_native_test_parts)`.
- **No winner moved.** Every `expect_winner` case in the solver-dispatch
  tests names the same solver under the ops estimate as under the deleted
  models, and the accuracy records are identical, line for line, to the merge
  base: selection changed, arithmetic did not.

## Deliberate divergences from FFTW
1. **No `cost_hook`.** FFTW lets a host redirect estimation; no NFFT caller
   has asked for it.
2. **No `nplan` / `epcost` statistics.** FFTW counts plans produced and
   estimates spent for its own diagnostics; nothing here reads them.
3. **`other` is weighted at one**, pending evidence. Our kernels are more
   bandwidth-bound than FFTW's codelets, so a higher weight is plausible, but
   moving it now would make any movement in the calibration numbers
   impossible to attribute.
4. **FFTW measures every fresh candidate at every node from `MEASURE`
   upward; our recursion still does not.** In FFTW, `evaluate_plan`'s
   `pcost == 0.0` disjunct is always true for a freshly built plan, so every
   fresh candidate is measured, and `BELIEVE_PCOST` only suppresses
   re-timing one that already carries a cost. Our `Y(planner_mkplan)` still
   selects purely on the estimate and measures nothing; only the top-level
   race measures. Measured children land with the threaded solvers.

## Consequences
- **`pcost` is two-valued and callers must not mix the units.** An estimate
  is an op count; a measurement is seconds or ticks. `Y(evaluate_plan)`'s
  return value is the grade, and the race grades a whole direction rather
  than comparing across units.
- **`pcost=` in `Y(fprint_plan)` is meaningful on every path**, measured or
  not — the untimed adopts carry their estimate.
- **The mock solvers in the tests moved from `pcost` to `ops.other`** — the
  natural field for a mock that does no arithmetic, and the only field whose
  weight is 1 regardless of `HAVE_FMA`, so the magnitudes are stable across
  hosts.
- **No wisdom vocabulary bump.** The `"plnr-flags"` tag in the configuration
  signature is unchanged: no `PLNR_*` bit was added, removed or renumbered,
  and no problem-hash field changed. `pcost` was never in the wisdom key —
  an entry stores a solver index, not a cost — so stored wisdom imports
  unchanged. The tag stays at its first value because no wisdom file has
  been exported by a released build.
- `docs/agents/pcost-calibration.md` is the calibration's document and was
  rewritten for the new estimator.
- **A clockless host now exports no wisdom from measured planning.** Every
  candidate there returns a negative `Y(plan_measure_cost)`, so the direction
  grades estimate and memoises unblessed — where the deleted restart used to
  rerun in estimate mode and bless that selection, nothing is blessed now.

## Known gaps
- **The estimate is precision-dependent through its FFT term — a permanent
  property, not a defect.** `FFTW(flops)` counts the plan the precision's own
  FFTW library built, and on this host `fftw3f` plans `n = 128` roughly 3.5x
  cheaper than `fftw3`/`fftw3l` (808 against 2880 counted operations), so
  the direct/fast estimate ratio at one geometry differs by precision. This
  is the estimate tracking what FFTW actually planned, which is the design.
  The consequence is that the race fixtures in `tests/nplan.c` are a
  per-precision calibration: their first choice (`M = 16`) held the ratio
  under the gate in double and long-double (3.60) but not float (6.58), and a
  float build failed three race-dependent tests. The fixture is now `M = 2`,
  which holds the worst precision at 1.58 — 2.5x under the 4.0 gate; the
  re-check procedure, with per-precision numbers, is in
  `docs/agents/pcost-calibration.md`.
- **`check_nplan_ndft_accuracy` fails in a float OpenMP build — pre-existing
  on the base branch, not caused by anything here.** It fails in float
  `checkall_ng_threads` only, on `err_planner <= err_legacy`. Measured on
  this host (N=8192, M=64, nodes near +1/2): the planner's direct NDFT
  returns the same error in every build (2.00e-06), but the legacy
  `X(trafo_direct)` compiles to a more accurate result under `-fopenmp`
  (1.60e-06 against 2.00e-06, thread count irrelevant — a compile-time
  codegen effect on the blocked phase recurrence, not a threading race), so
  the comparison flips. `kernel/nfft/nfft.c` is identical to the merge base
  and the assertion's text is unchanged there, so the base fails it
  identically; it surfaced only because this branch built float for the first
  time. Out of scope here: it needs its own ruling — a tolerance in the test
  or an investigation of the `-fopenmp` codegen difference.
- The measured-mode fallback leg of `Y(evaluate_plan)` has no deterministic
  unit test: reaching it needs a machine without a fine clock or an `apply`
  the clock cannot resolve. It is covered indirectly on clockless hosts by
  the `check_nplan_measured` family's known-skip path.
- The estimator's residual 1.43x direct/fast unit-scale gap is a property of
  this host's operation mix, not of the model; re-run the sweep on a
  different host before trusting the ratio there.
- `prob[ADJ]` is still never built: `Y(evaluate_plan)`'s estimate and
  measurement are both forward quantities, commensurable with each other and
  with the race. Separate per-direction solvers remain the way out, as in
  ADR-0005's accepted limitation.

## Not decided here
- Weighting `other` above one (the sweep's evidence, not judgement).
- Measuring children in the recursive search, as FFTW does.
- Per-direction solvers, which would give the adjoint its own race.
