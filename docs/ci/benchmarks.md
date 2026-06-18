# Running benchmarks in CI

CodSpeed benchmarks (`.github/workflows/bench-linux.yml`) run on expensive
macro runners with a limited monthly allowance, so on pull requests they are
**opt-in behind a maintainer approval**.

## When benchmarks run automatically

- On every push to `main`, `develop`, or `master`. These keep CodSpeed's
  comparison **baseline** current and need **no** approval. They happen on
  merge, so they are infrequent.

## How to run benchmarks on a pull request

1. Every PR push creates a benchmark run that immediately enters a **"Waiting"**
   state. It consumes **no** runner minutes while waiting.
2. A listed reviewer opens the run (Checks / Actions tab) and clicks
   **"Review deployments" → `benchmarks` → "Approve and deploy."**
3. The benchmark runs against the PR's merge commit, and the **CodSpeed** check
   (regressions/improvements vs. the base-branch baseline) posts automatically.

Only the maintainers listed as **required reviewers** on the `benchmarks`
environment can approve — see *Settings → Environments → benchmarks*.

## What happens to older, un-approved runs

Each new push **supersedes** the PR's previous not-yet-approved run (per-ref
concurrency with `cancel-in-progress`), so there is always exactly **one** run
to approve: the latest commit's. You don't need to clean up stale waiting runs.

## Benchmarks are not a merge gate

Benchmarks are advisory. If no one approves, the run never executes and there is
no benchmark status on the PR, so it never blocks merging. Maintainers decide
per PR whether to run benchmarks and whether to merge regardless of the result.

## Security notes

- The workflow uses `pull_request` (never `pull_request_target`): fork PR code
  runs with a **read-only token and no secrets**. The approval gate means a
  named maintainer explicitly authorizes each fork-code execution.
- The workflow holds **`contents: read`** only — no write token.
- No `CODSPEED_TOKEN` secret is used (tokenless CodSpeed upload on this public repo).

## Manual / off-PR runs

Use **Actions → Bench Linux → Run workflow** (`workflow_dispatch`) to benchmark
any branch without a PR; these are not gated.
