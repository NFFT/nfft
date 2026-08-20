---
name: deslop
description: Use when writing, reviewing, or cleaning up code commentary and file-local static names in this repo, and always before committing agent-generated C code. Triggers - "deslop", "trim the comments", "too verbose", comments that restate the code, carry task/item/spec/review numbers, porting provenance, design speculation, banners, or explain C basics to experienced developers.
---

# Deslop

## Overview

Comments in this repo are read by two audiences: experienced C developers, and
coding agents that must find relevant code and trace calls fast. Neither needs
the code explained. Both need the facts that are NOT in the code.

**Rule: the shortest comment that keeps every non-obvious fact wins.**

Agent-written C is reliably too long. Deslopping is a mechanical pass: delete
the categories in Cut, keep the categories in Keep, then check the length
budget.

## Cut

Each category below is a real before/after from this repo.

### 1. Process archaeology

Task numbers, item numbers, spec sections, review IDs, phases, dates,
decision-meeting references, tracker filenames. The reader has git and the
issue tracker.

```c
/* Task 9: adjoint entry. Forward-only race selects one plan (dir[FWD]) ... */
/* Item 15: awaken the candidate to AWAKE_ZERO ... */
/* ESTIMATE mode: spec 08 section 7 / spec 07 amended. Bounds ... */
/* closes, just for this one bit (2026-07-30 /grilling decision). */
/* Review D2: the quadrant split hinges on N/2 ... */
/* Fast-NFFT decomposition problem types (IMPROVEMENTS.md item 13). */
```
becomes
```c
/* Adjoint entry. Forward-only race selects one plan (dir[FWD]) ... */
/* Awaken the candidate to AWAKE_ZERO ... */
/* ESTIMATE mode. Bounds ... */
/* closes, just for this one bit. */
/* The quadrant split hinges on N/2 ... */
/* Fast-NFFT decomposition problem types. */
```

Also cut the schedule words that go with them: "this item", "yet", "today",
"a future task may", "deferred (not ported here)". Keep only a real
`/* TODO: ... */` for work someone must still do.

### 2. Porting provenance

Where code came from, which lines, and that it was reformatted.

```c
/* Forward B (g -> f[j]): ported from nfft_trafo_2d_compute
 * (kernel/nfft/nfft.c:2931-3008), then reformatted to the repo BSD style;
 * the 4-case wrap-split accumulation. */
```
becomes
```c
/* Forward B (g -> f[j]) */
```

Line-number citations rot. Name the function only when a reader must still
consult it.

### 3. Restating the code

```c
/* precompute the psi table from ego->x (borrowed, node-dependent); guarded by
 * `precomputed`, invalidated back to SLEEPY -- mirrors conv_awake (1D). */
```
becomes
```c
/* precompute the psi table from ego->x. */
```
The guard flag, the invalidation, and the sibling implementation are all
visible in the ten lines below.

### 4. Explaining C, or the local API, to an expert

```c
/* Integer floor(log2(v)) for v >= 1 (M-bucketing in the key; local copy of the
 * helper in nproblem.c -- only CONV needs it). */
/* the virtual count register (fixed-frequency, frequency-scaling invariant) */
```
Delete outright, or reduce to the name of the thing.

### 5. Justification aimed at a reviewer

Defence of a choice, "deliberately", "the CAPSTONE", "no exceptions",
"(Defense-in-depth; A() is debug-only)", enumerated alternatives not taken.
The decision is made. State the fact, not the argument.

### 6. Parenthetical qualifiers that add nothing

```c
/** forward transform (planner-chosen algorithm) */   -> /** forward transform */
/* d == 2 only; KB only (only implemented window). */ -> /* d == 2 only */
/* x is a borrowed caller array -- never freed here. */ -> /* x is a borrowed caller array */
/* kind stays PLAN_KIND_NATIVE (Y(plan_create) default): self-contained. */ -> delete
```

### 7. The same fact stated twice

State a fact once, at the highest site that covers every reader who needs it,
then delete the copies below it.

- Repeated across a group of accessors or callbacks -> one section header above
  them, which must then carry the shared fact and not merely name the group.
- Repeated in a helper and in its only caller -> the helper.
- Repeated in a caller and in the constraint it enforces -> the constraint.

```c
/* direction-aware geometry accessors */
INT Y(problem_nfft_N)(const problem *p, int t) {
  /* t indexes the axis list in caller order; valid range is [0, sz->rnk).
   * On the top-level copy_x=1 path this is the unit-elided (compressed)
   * live-axis list (rnk <= d); on the borrowed copy_x=0 path all d axes are
   * retained. */
INT Y(problem_nfft_n)(const problem *p, int t) {
  /* ... the same four lines again ... */
```
becomes
```c
/* Direction-aware geometry accessors. Forward: n_in = N_t, n_out = n_t;
 * adjoint (stored after tensor_adjoint): n_in = n_t, n_out = N_t. t indexes
 * the axis list in caller order, range [0, sz->rnk) -- the unit-elided live
 * axes on the copy_x=1 path, all d axes on the borrowed path. */
```

Do not compromise by keeping a short restatement at the lower site. Either the
higher site covers it, or it does not belong there.

### 8. Duplicated prose across sibling files

Rank-specialised or otherwise parallel files get ONE shared header text,
identical word for word except the dimension. Do not re-derive the shared
explanation per file.

```c
/* 3D CONV solver: Step C of the fast NFFT decomposition -- the node convolution
 * (matrix B). Sums the oversampled grid g against the window psi at each
 * nonequispaced node x_j (forward), or scatter-adds f onto g with the same psi
 * weights (adjoint). psi depends on x/window/n/N/m, precomputed once at awake
 * (sparse PRE_PSI strategy in legacy code). */
```

### 9. Shouting and asides

ALL-CAPS for emphasis, stacked em-dash clauses, rhetorical contrast
("never one transform against another"). Keep capitals only where they mark a
real term (`AWAKE_ZERO`, `PRE_PSI`, `NULL`). Emphasis capitals survive a first
sweep easily because they sit inside otherwise-good sentences - grep for them:
`grep -nE "\b[A-Z]{2,}\b" <files>`, then keep only the real terms.

```c
/* Refresh thread count BEFORE any keying. */   -> before
/* measure each candidate on its OWN core */    -> own
/* the RAW fftw_flags still reaches ... */      -> raw
```

### 10. Narrating what does not happen

A comment whose whole content is that this site is uneventful. It reads as
information and carries none; every one of these is deleted outright, never
shortened.

```c
/* Every winner is a coreless native plan. Nothing to build or install here. */
/* The estimate winners already own their cores. */
/* no shared bundle core to finalize here (every plan is coreless native) */
```

Same for the reviewer-facing defence of a defensive check:

```c
/* Required args. A() is a no-op in release builds, so validate explicitly and
 * return NULL (FFTW's contract) instead of dereferencing. Checked before any
 * allocation, so a NULL return leaks nothing. */
```
becomes nothing - the `if (... == 0) return 0;` below states it.

### 11. Design speculation

What a future variant would do, what the current shape leaves room for, what
would not be a rewrite. State the dormant thing once, at the declaration that
is dormant, and never again.

```c
/* [ADJ] is left dormant (NULL) so a future two-internal-plan option is
 * "populate prob[ADJ] + race it", not a rewrite. */
/* sign is RETAINED as the direction-dependent-wisdom hook for that future
 * two-problem split. */
/* The (fwd)/(adj) framing is kept -- the [ADJ] substrate is dormant, not
 * removed -- so a future two-internal-plan option fills in (adj %p) without a
 * format change. */
```
becomes, once, on the struct member:
```c
plan *dir[2]; /* [FWD] the winning plan; [ADJ] dormant (NULL) */
```

### 12. Banners and step numbers

The `/*---*/` sandwich around a section title, column-padded to the right
margin, and ordinal prefixes inside a linear function. A section title earns a
line only when it states a fact the function names below do not.

```c
/*-----------------------------------------------------------------------*/
/* Destroy                                                                */
/*-----------------------------------------------------------------------*/
void Y(plan_ng_destroy)(...)          -> delete all three lines

/* Guru constructor                                                       */
                                      -> /* Guru constructor */

/* 1. Ensure the roster is registered into the process-global planner. */
                                      -> /* Ensure the roster is registered. */
```

## Names

A comment that exists to explain a naming convention is a defect in the names.
Fix the names and the comment disappears with it.

```c
/* ADT callbacks use the pk_ prefix: Y(problem_*) expands to nfft_problem_* in
 * this translation unit. */
static void pk_hash(const problem *p, md5 *ctx) {
```
becomes
```c
static void hash(const problem *p, md5 *ctx) {
```

A `static` function is already scoped to its file. Name it for its role, not
its file, its type, or its adt slot: `apply`, `apply_adjoint`, `print`,
`destroy`, `live_axes` - not `rnk0_apply`, `nfft_live_axes`. Exported names
keep their `Y(...)` / `X(...)` mangling as `CONVENTIONS.md` requires; this rule
is about file-local statics only. Renaming a static is a mechanical, in-file
change - do it in the same pass, then rebuild.

## Mechanics

- One space after a full stop, never two.
- No column padding inside or after a comment; no trailing whitespace.
- Do not leave a hand-trimmed comment mid-sentence or split across a stale
  line break. Rewrap the block after every cut.

## Keep

Trimming is not blanket deletion. These earn their lines, and if they are
missing you ADD them:

- **Non-obvious invariants, with the reason.** Added during the cleanup:
  `/* The scatter accumulates (+=) into an overlapping, node-dependent set that
  does not cover the grid, so the whole grid must start zeroed. */`
- **Traps.** `/* x is the COMPRESSED copy (rnk*M reals, not d*M). Sizing the
  guard from the original d would over-read out of bounds in debug builds. */`
- **Ownership and lifetime.** Who allocates, who frees, what is borrowed.
- **Applicability gates.** `/* d >= 4 only. */`
- **Memory layout and index conventions.** `psi[(j*d+t)*(2m+2)+lj]`, frequency
  range of a table, row-major vs adjoint-flipped.
- **Architecture in one block at the top of a file**: what the unit is, its
  place in the decomposition, its dataflow line.
- **Trace anchors for agents.** A name the reader must actually go and read:
  `/* guards_ok (nfast.c) declines odd N upstream */`, a roster comment listing
  the leaf files it registers, the header a prototype lives in. Prefer a name
  over a line number. An anchor naming something merely related, that no reader
  needs to open, is category 6 - cut it.
- **`TODO:`** for genuine open work.

## Length budget

| Site | Default | Hard cap |
|---|---|---|
| File header | 3-5 lines | 10, only with dataflow or ownership content |
| Function | 1 line | 3 |
| Inline | 1 line | 3, only for a trap or an invariant |
| Doxygen in a public header | 1 clause | 2 lines |

Over budget is allowed only when every extra line carries a fact you cannot
get by reading the function.

Public-header doxygen states the contract, not the mechanism:
```c
/** import wisdom from a file; returns 1 on success, 0 on failure */
/** build the node-dependent tables, call after filling x */
```

## Procedure

1. Grep the target files for the markers:
   `grep -nE "item [0-9]|Task [0-9]|spec [0-9]|specs/|[Rr]eview [A-Z][0-9]|Phase [0-9]|ported from|IMPROVEMENTS|deliberately|mirrors |[Nn]othing to |already own|future |would be" <files>`
2. Delete every hit's clause; keep the sentence if a real fact survives.
3. Re-read each remaining comment and ask: does the code below already say
   this? If yes, delete.
4. Apply the length budget. Merge what is left into one block per site.
5. Add the missing Keep items - invariants, traps, ownership - that the long
   version buried.
6. For every fact stated more than once, keep the copy at the highest covering
   site and delete the rest. Check sibling files use the same header text.
7. Rename the file-local statics whose names carry a redundant prefix, and
   drop the comments that explained those prefixes.
8. Rebuild the touched directories. The pass is comment- and static-name-only,
   so a compile error means you broke something.
9. Review the result (below), then sweep again on the confirmed findings.

## Review the sweep

A first sweep is reliably not enough. It cuts the loud slop - the item numbers,
the ported-from lines - and leaves the quiet slop: emphasis capitals inside
good sentences, "nothing happens here", the same design speculation restated
at four sites, a section banner the eye reads as structure.

So after the sweep, dispatch a fresh subagent to review it. A fresh context is
the point: you have spent the whole sweep reading these comments and can no
longer see them.

```
Agent(subagent_type: "Explore", description: "Review deslop pass", prompt: """
Read .claude/skills/deslop/SKILL.md, then `git diff <base>..HEAD -- <files>`
and the current text of <files>.

Report every comment that still violates the skill, and every comment the
sweep deleted that carried a fact not recoverable from the code (an
over-cut). For each: file, line, the comment text, the Cut category number or
Keep bullet it offends, and the one-line replacement you propose. Cite the
category; do not report a comment you cannot name a rule for.

Rank by confidence. Report nothing you are unsure about.
""")
```

Then act on the findings yourself: apply the ones that hold, and say which you
rejected and why. A finding is not a mandate - the reviewer has not read the
surrounding module, so an "over-cut" it flags is often a fact that lives at a
higher site on purpose.

Repeat only if the review found something substantial. Two rounds is normally
the end of it.

## Red flags

Any of these means the comment is slop:

- A number that only means something inside a planning document.
- A file path with a line range.
- "mirrors", "exactly as", "just like", followed by a paragraph.
- A sentence that names a variable and then says what the next line does to it.
- A parenthesis that says what the code is NOT, or does NOT do.
- A whole comment whose content is that nothing happens here.
- A word in capitals that is not a real term.
- The word "future", or a description of a variant that does not exist.
- A `/*---*/` banner, or a comment padded out to the right margin.
- An ordinal ("1.", "Step 2") inside a single linear function.
- A comment that explains what a prefix in a name means.
- The same fact in a section header and again in the block below it.
- A section header that only repeats the name of what it sits above.
- Two sibling files with different prose for the same role.
- A comment longer than the function it documents.

## Rationalizations

| Excuse | Reality |
|---|---|
| "The context helps a future maintainer" | Git and the tracker hold context. Code holds facts. |
| "An agent reading this needs the background" | An agent needs names to grep and invariants to respect. Prose slows it. |
| "It documents why I chose this" | Ship the decision, not the deliberation. |
| "Line numbers make it traceable" | They are wrong at the next commit. Use names. |
| "Each file should stand alone" | Sibling files share one header text. Divergent prose reads as a difference that is not there. |
| "It is only a few extra lines" | Per file, times 178 files. That is the problem being fixed. |
| "That comment is fine, I just wrote it" | You wrote it inside the same context that produced the slop. Have a fresh agent review the sweep. |
| "One sweep was enough" | It never is. The first pass cuts the loud slop and leaves the quiet slop. |
