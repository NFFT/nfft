# Preflight — pick the measurement track

*[← Overview & map](../REFERENCE.md) · Next: [Phase A — baseline](phase-a-baseline.md)*

**Do this first, once, and record the result** — it decides the build mode in Phase A
and the performance metric in Phases C/D/E. The correctness signal (tests) is the same
either way; only how *performance* is measured changes.

Check, in order, whether the **CodSpeed track** is fully available *for this repo*:

1. **MCP present** — are the CodSpeed MCP tools (e.g. list/compare runs) actually in
   your toolset this session?
2. **Authenticated** — `codspeed status` reports logged in (not *"Not logged in"*).
3. **Repo onboarded with data** — an MCP *list-runs* / *list-repositories* call shows
   **this** repository with at least one run on the base branch to compare against.

- **All three yes → CodSpeed (simulation) track.** Build `-DNFFT_BENCHMARK_MODE=simulation`;
  the metric is the deterministic instruction count (per-case via `codspeed run`, and
  the base branch via the MCP). This is what CI gates on, so local and CI agree.
- **Any no → local (walltime) track.** Build `-DNFFT_BENCHMARK_MODE=walltime`; the
  metric is `median_ns` under the noise rule in
  [Working without CodSpeed](measurement-modes.md#working-without-codspeed). Fully
  offline, no account. Less precise, but a complete, supported loop.

State which track you're on at the start of the report; everything below adapts to it.
The default the rest of these docs show is the **local track** (the always-available
one); where the CodSpeed track differs, it is called out inline.

## Deliverables (exit criteria)

The task dir and its tracker `README.md` already exist — Step 0 created
`docs/perfeng/NNNN-<target-slug>/` (worked example `0001-trafo-direct`) before
Preflight. This phase only records the track decision; it writes no `artifacts/`.

Copy [`../templates/preflight.md`](../templates/preflight.md) to
`docs/perfeng/0001-trafo-direct/preflight.md` — a short record (≈3 lines, don't pad
it) stating the chosen track (`walltime` | `simulation`) and the **evidence** —
the three checks in order: MCP present? `codspeed status` logged in? repo onboarded
with a base-branch run? Three yes ⇒ `simulation`; any no ⇒ `walltime`, naming the
failing check.

Then fill the tracker's **Track** field (header in
[deliverables.md](deliverables.md)) with that track, and flip the **Preflight** row
to ✅.

*Deliverable = exit gate:* Preflight is not exitable until `preflight.md` exists
**and** the tracker's **Track** field is set and its Preflight row is ✅.

*[← Overview & map](../REFERENCE.md) · Next: [Phase A — baseline](phase-a-baseline.md)*
