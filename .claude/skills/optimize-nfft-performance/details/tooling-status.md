# Tooling status (agent-operability)

*[← Overview & map](../REFERENCE.md) — cross-cutting reference, consult from any phase.*

Verified end to end in the dev container, using the CMake tree throughout. The full
five-phase loop was run on the worked example (`X(trafo_direct)`, optimised by phase
recurrence): the fault enumerated a 149-case net, the slowdown isolated the
`forward_direct_1d` metric, and the change landed at ~1.3× wall-clock / ~2.3×
instructions with the net green — the gaps that surfaced are folded into the notes
above.

- **Phase A is fully agent-operable.** `cmake --build` + `ctest` / direct
  `build-cmake/tests/checkall`, the `-> FAIL` stdout signal, exit codes, and the CUnit
  XML all work with no human step. Fault-inject → observe → revert was exercised end to
  end (67 `-> FAIL` lines on the seeded fault; 0 after revert).
- **Phase B (local, walltime) is agent-operable.** `cmake -B build-cmake
  -DNFFT_BENCHMARK_MODE=walltime` built and ran offline end to end: `FetchContent`
  fetched codspeed-cpp `v2.3.0` (submodules auto-recursed), and the binary wrote a
  per-case stats JSON (`median_ns`, …) with no runner, token, or upload. The
  `simulation` build also measures under `valgrind --tool=callgrind`. `valgrind` and
  the `codspeed` CLI ship in the dev container (`.devcontainer/Dockerfile`).
- **Clean per-case simulation needs a CodSpeed account.** The CodSpeed CLI (`codspeed
  run`/`exec`, in the dev container) validates a token before running — there is **no
  offline mode** in 4.17.5 (`codspeed run` with no/empty token fails immediately;
  verified). Account-free, the only local simulation is a raw `valgrind` single-case
  **process total** (startup-contaminated, rough). So with no CodSpeed access the loop
  runs on **walltime alone** (see [Working without CodSpeed](measurement-modes.md#working-without-codspeed))
  — the supported baseline, accepting timing noise.
- **The simulation instrument needs CodSpeed's Valgrind fork**, not stock valgrind
  (`codspeed setup status` otherwise reports *"not a CodSpeed build"*). The dev
  container installs it (`.devcontainer/Dockerfile` runs `codspeed setup --mode
  simulation`, which apt-installs `valgrind-…codspeedN` over `/usr/bin/valgrind`); no
  account is needed for `codspeed setup` itself, only for `codspeed run`.
- **Optional — CodSpeed MCP for CI-history parity.** To compare against the base
  branch's CI numbers in Phase D, register the **CodSpeed MCP server**
  (`npx add-mcp https://mcp.codspeed.io/mcp --name CodSpeed`, or the Claude Code plugin
  `CodSpeedHQ/codspeed`). It exposes tools to list/compare runs and read flamegraphs.
  **Prerequisites the user must provide:** a CodSpeed account, the repo onboarded to
  CodSpeed, and CI uploading results — so this is an *enhancement*, not a requirement,
  and is **not wired by default**. Without it the Phase-D simulation check falls back to
  a rough raw-`valgrind` single-case cross-check, or is skipped (walltime-only).
- **The Autotools benchmark path is legacy** (`./configure --enable-benchmarks
  --with-codspeed=<path>` + `make bench`): it needs a hand-built codspeed-cpp and
  prints no measurements. AGENTS.md §4 and the loop above use the CMake build instead.
- **No benchmark covers the fast path** — the substantive gap for *real* work. See
  [caveats](caveats.md): only the `*_direct` transforms are benchmarked (the worked
  example, `trafo_direct`, is among them). The *fast* path (`trafo` / `adjoint` and the
  `precompute_one_psi` strategies) and `init`-only helpers like `intprod` sit outside
  every measured region; closing perf work there requires *adding* benchmarks first.
