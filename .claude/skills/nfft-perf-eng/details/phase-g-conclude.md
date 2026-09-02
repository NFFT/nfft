# Phase G — conclude (squash, push, PR, attach deliverables)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase F — exit gate](phase-f-exit-gate.md)*

The optimization is finished — Phase F passed (or ended `reverted`/blocked) and every
deliverable is written under the gitignored `.perfeng/`. Phase G hands the work off: it turns
the run's *source* history into one clean commit and offers to publish it, carrying the
deliverables to the reviewer as an attachment rather than as repo files.

Two things are true by now and shape this phase:

- **Deliverables were never committed.** `.perfeng/` is gitignored, so the run's commits only
  ever captured *source* changes (see [deliverables.md](deliverables.md#where-it-lives)). The
  deliverables travel to the reviewer as the **zip attached to the PR**, not in the tree.
- **You committed as you went.** During the loop each self-contained unit of work landed as
  its own commit ([the inner loop](phase-e-inner-loop.md) — a kept optimization step, a
  permanent test addition). Phase G collapses those into one.

The archive is named for the PR, and the PR number only exists once the PR is created — after
the squash. So `perf-conclude.sh` has **two steps**: `squash` (before the PR) and `package
<pr>` (after it), with the interactive push/PR in between.

## G1. Squash the run into one commit (deterministic)

Collapse every commit since the [recorded base](deliverables.md#step-0--open-the-perfeng-directory-gated)
(`.perfeng/BASE`, stamped at Step 0) into a single commit. This is scripted — it preflights,
then squashes non-interactively (`git reset --soft <base>` + one commit; interactive rebase is
unavailable in this environment):

```bash
SCR=.claude/skills/nfft-perf-eng/scripts
$SCR/perf-conclude.sh squash -m "Optimize <target>: <one-line what + speedup>"
```

`perf-conclude.sh squash` **refuses** unless the working tree is clean (all source work
committed — the gitignored `.perfeng/` is never in the diff), the current branch is **not** the
default branch, and the recorded base is an ancestor of HEAD. It squashes to one commit and
prints the exact push/label/PR/package commands.

If you are on the default branch, create a feature branch first (`git switch -c perf/<slug>`)
and re-run — the squash and PR belong on a branch.

## G2. Offer to push, then optionally open a PR

Both steps are **opt-in — confirm with the user before running them.** The loop is autonomous;
publishing is not.

1. **Push** the branch: `git push -u origin <branch>` (add `--force-with-lease` if the branch
   was pushed *before* the squash rewrote its history).
2. **Open a PR** to the default branch (`develop`), labelled `perf-eng`. Create the label once
   if the repo doesn't have it yet:
   ```bash
   gh label list | grep -q '^perf-eng' || \
     gh label create perf-eng --description "nfft-perf-eng optimization run" --color 0e8a16
   gh pr create --base develop --label perf-eng \
     --title "Optimize <target>: <speedup>" --body "<summary — link the phases, cite the numbers>"
   ```
   The PR body is the reviewer's entry point: state what changed and why, the measured speedup
   (per precision), the Phase-D accuracy verdict, and any `residual` risk that made the run
   `partial`. `summary.html` (in the attached zip) is the full walkthrough.

## G3. Package the deliverables for the PR, and attach

Once the PR exists and you know its number `N`, package the deliverables and attach them so the
reviewer gets the complete evidence bundle without it living in the repo:

```bash
$SCR/perf-conclude.sh package <N>          # e.g. package 231
```

`package` **renames `.perfeng/` → `.perfeng-pr-<N>/`** (no leading zeros) and zips *that* to
`perfeng-pr-<N>.zip` outside the tree, as a **standard ZIP** (default DEFLATE — any
unzip/Explorer/Finder opens it, no exotic codecs). Naming the top-level directory for the PR
means a downloaded archive **unpacks into a PR-unique directory** — unambiguous within this repo
even when several runs' result zips are downloaded side by side. (`.perfeng-pr-*` is gitignored too.)

> **`gh` cannot upload a binary attachment to a PR or a PR comment** (its `--body`/`--body-file`
> are text only). Attach the printed zip by **drag-dropping it into a PR comment in the web UI**
> — a short manual step inside the already-interactive publish flow. Give the user the archive
> path and the PR URL to complete it. (Don't work around this by committing the zip: the
> deliverables are gitignored on purpose.)

## G4. Follow-up after conclude — keep the attached archive current

Concluding is not always the last word: a reviewer's follow-up question can trigger more code
changes *and* deliverable edits after the PR is open. When that happens, the attached archive
must not go stale — **the zip on the PR always reflects the final state.**

After any post-conclude follow-up work:

1. Commit the follow-up **source** change on the branch (deliverables stay gitignored, as always).
   These are ordinary follow-up commits — don't re-squash a pushed, under-review branch.
2. Update the deliverables in `.perfeng-pr-<N>/` (the run directory, already renamed by the first
   `package`) — the summary/phase docs/artifacts that the change affects.
3. **Re-package:** `perf-conclude.sh package <N>` again. With the dir already renamed it
   re-packages in place and rebuilds `perfeng-pr-<N>.zip` (same name, standard ZIP).
4. **Overwrite the archive on the PR** — because attaching is manual: edit the existing PR
   comment that holds the old `perfeng-pr-<N>.zip`, remove that attachment, and drag-drop the
   refreshed zip into the **same** comment. One comment, one current archive — not a pile of
   stale versions. (`perf-conclude.sh package` prints exactly this instruction on a re-package.)

## Exit — the run is fully concluded when

- the run is **one commit** on a feature branch (`perf-conclude.sh squash` ran clean);
- if the user opted in: the branch is **pushed**, a **PR to `develop`** exists with the
  **`perf-eng`** label, and the **`perfeng-pr-<N>.zip` archive is attached** to it (from
  `perf-conclude.sh package <N>`);
- if the user declined publishing: the squashed commit is reported and the deliverables stay in
  `.perfeng/`, so they can push/PR and package later;
- if follow-up work landed after conclude: the code is committed, the deliverables refreshed, and
  the attached `perfeng-pr-<N>.zip` re-packaged and **overwritten in its PR comment** (G4) — the
  archive on the PR matches the branch's final state.

A `reverted` or hard-gate-blocked run still concludes here — there may be nothing to squash
(no source change landed), but `summary.html` and the zip still document why it stalled; offer
the PR only if there is a change worth reviewing (e.g. a permanent test that closed a coverage
gap).

*[← Overview & map](../REFERENCE.md) · Prev: [Phase F — exit gate](phase-f-exit-gate.md)*
