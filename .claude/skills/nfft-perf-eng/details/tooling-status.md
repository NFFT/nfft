# Tooling status (agent-operability)

*[← Overview & map](../REFERENCE.md) — cross-cutting reference, consult from any phase.*

Verified end to end in the dev container, using the CMake tree throughout. The full
five-phase loop was run on the worked example (`X(trafo_direct)`, optimised by phase
recurrence): the fault enumerated a 149-case net, the slowdown isolated the
`forward_direct_1d` metric, and the change landed at ~1.3× wall-clock with the net
green — the gaps that surfaced are folded into the notes above.

- **Phase B is fully agent-operable.** `cmake --build` + `ctest` / direct
  `build-cmake/tests/checkall`, the `-> FAIL` stdout signal, exit codes, and the CUnit
  XML all work with no human step. Fault-inject → observe → revert was exercised end to
  end (67 `-> FAIL` lines on the seeded fault; 0 after revert).
- **Phase C (walltime) is fully agent-operable, no account.** `cmake -B build-cmake
  -DNFFT_BENCHMARK_MODE=walltime` built and ran offline end to end: `FetchContent`
  fetched codspeed-cpp `v2.3.0` (submodules auto-recursed), and the binary wrote a
  per-case stats JSON (`median_ns`, …) with no runner, token, or upload. This is the
  metric the whole loop uses, and the same mode CI runs on its CodSpeed Macro Runners —
  so local results preview the CI gate directly.
- **The callgrind tie-breaker needs CodSpeed's Valgrind fork**, not stock valgrind
  (`codspeed setup status` otherwise reports *"not a CodSpeed build"*). The dev
  container installs it (`.devcontainer/Dockerfile` runs `codspeed setup --mode
  simulation`, which apt-installs `valgrind-…codspeedN` over `/usr/bin/valgrind`); **no
  account** is needed for `codspeed setup` or for a raw `valgrind --tool=callgrind` run.
  Its only role is the optional deterministic tie-breaker for an ambiguous untouched
  control case (raw single-case `I refs` delta, layout vs real work — see
  [measurement-modes.md](measurement-modes.md#deterministic-tie-breaker-optional--is-a-stuck-control-case-real-work-or-layout)).
  It is **not** a routine metric and **not** a gate; clean per-case counts (`codspeed
  run`) would need a token and are out of scope.
- **The Autotools benchmark path is legacy** (`./configure --enable-benchmarks
  --with-codspeed=<path>` + `make bench`): it needs a hand-built codspeed-cpp and
  prints no measurements. AGENTS.md §4 and the loop above use the CMake build instead.
- **No benchmark covers the fast path** — the substantive gap for *real* work. See
  [caveats](caveats.md): only the `*_direct` transforms are benchmarked (the worked
  example, `trafo_direct`, is among them). The *fast* path (`trafo` / `adjoint` and the
  `precompute_one_psi` strategies) and `init`-only helpers like `intprod` sit outside
  every measured region; closing perf work there requires *adding* benchmarks first.
