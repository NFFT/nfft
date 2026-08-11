---
name: understanding-fftw
description: Use when you need a precise FFTW3 C function signature, planner-flag/alignment/execute-function detail, license text, or exact reference-manual wording — instead of relying on memory or a live web fetch.
---

# Understanding FFTW

## Overview

Offline, local mirror of the upstream FFTW3 project: a shallow clone of its
source (`github.com/FFTW/fftw3`) plus its reference manual (`fftw3.pdf`)
converted to Markdown. Grep it for ground truth instead of guessing from
training data or depending on network access.

The mirror is cached **inside the current project**, at
`<repo-root>/.fftw3-reference/`, not inside this skill directory —
run the fetch scripts below the first time you need it in a given repo. The
scripts also add that path to the repo's `.git/info/exclude` (never to a
tracked `.gitignore`), so it's ignored locally without affecting anyone
else's clone.

## When to use

- Confirming an exact FFTW3 public API signature, macro, or type
- Checking FFTW's license text or source layout
- Looking up reference-manual sections (planner flags, alignment rules,
  new-array execute functions, guru interface, MPI/threads, ...)
- Any FFTW-related question — e.g. while working on NFFT3, which depends on
  and is influenced by FFTW — where source- or manual-level precision matters
  more than a paraphrase

## Contents

Under `<repo-root>/.fftw3-reference/`:

- `fftw3-src/` — depth-1 clone of `https://github.com/FFTW/fftw3.git`.
  Public API lives in `api/fftw3.h`, generated via `X()`/`R`/`C` mangling
  macros (precision-agnostic, one source per float/double/long-double —
  analogous to NFFT3's own `Y()`/`X()` convention).
- `fftw3-docs.md` — `fftw3.pdf` converted to Markdown (headings and code
  blocks preserved; not pixel-perfect, see Common mistakes below).

## How to use

First time in a given repo (or to check it exists), fetch/update:

```bash
.claude/skills/understanding-fftw/scripts/update.sh          # both
.claude/skills/understanding-fftw/scripts/fetch-source.sh    # source only
.claude/skills/understanding-fftw/scripts/fetch-docs.sh      # docs only
```

Run these from anywhere inside the target project's working tree — they
resolve the project root via `git rev-parse --show-toplevel`. Not being
inside a git repo is an error (the script has nowhere sensible to cache to).

Then grep, don't read whole files, from `<repo-root>/.fftw3-reference/`:

```bash
grep -n "plan_many_dft_r2c" fftw3-src/api/fftw3.h
grep -n -i "new-array execute" fftw3-docs.md
```

The header uses unexpanded macros, so search for the name **without** the
`fftw_` prefix (e.g. `plan_many_dft_r2c`, not `fftw_plan_many_dft_r2c`).

## Quick reference

| Need | Where |
|---|---|
| Exact function signature | `.fftw3-reference/fftw3-src/api/fftw3.h` |
| License | `.fftw3-reference/fftw3-src/COPYING` |
| Build/config options | `.fftw3-reference/fftw3-src/configure.ac` |
| Manual section (planner flags, alignment, execute functions, ...) | `.fftw3-reference/fftw3-docs.md` |
| Precision-mangling macros | `.fftw3-reference/fftw3-src/api/fftw3.h` (search `X(`, `R`, `C`) |

## Refreshing

Cached content is a point-in-time snapshot. Re-run the fetch scripts above
(safe to re-run any time — they update in place) when you need the latest
upstream state. `fetch-docs.sh` tries `pandoc` first, falls back to
`uvx`+`pymupdf4llm`, then `pdftotext`.

## Common mistakes

- Grepping `api/fftw3.h` for `fftw_`-prefixed names — the header uses
  unexpanded `X()` macros; search the unprefixed name.
- Assuming the Markdown conversion is typographically perfect — pymupdf4llm
  occasionally merges hyphenated words split across a line break (e.g.
  "new-array" → "newarray"). If a grep misses, retry with a shorter
  substring or without the hyphen.
- Treating `fftw3-src` as having full history — it's a depth-1 shallow
  clone (latest default-branch commit only). Use the upstream GitHub repo
  directly for blame/history questions.
- Looking for `.fftw3-reference/` inside this skill directory — it lives in
  the *project* you're working in, not here. If it's missing, run the fetch
  scripts first.
- Committing it by hand (e.g. `git add -A`) — it's excluded via
  `.git/info/exclude`, so normal `git add`/`git status` never see it; no
  extra care needed, but don't force-add it either.
