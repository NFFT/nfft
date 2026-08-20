# Planner internals (the kind-agnostic core)

The `kernel/planner/` trinity is **kind-agnostic** — it knows nothing about
NFFT. NFFT-specific problems/solvers live in `kernel/nfft/`, `kernel/deconv/`,
`kernel/conv/` and register into the global planner. Ground truth:
`include/iplanner.h` + `kernel/planner/*.c`.

## The trinity ADTs (`include/iplanner.h`)

Each family is a struct whose first member is a const vtable pointer; concrete
instances embed the base struct first.

```c
/* problem — hashable "what" */
typedef struct {
  int kind;
  void (*hash)(const problem *p, md5 *ctx);
  void (*print)(const problem *p, printer *pr);
  void (*destroy)(problem *p);
} problem_adt;
struct problem_s { const problem_adt *adt; };

/* plan — executable "how" */
typedef struct {
  void (*apply)(const plan *ego, const problem *p);          /* forward  */
  void (*awake)(plan *ego, int wakefulness);                 /* may be NULL */
  void (*print)(const plan *ego, printer *pr);
  void (*destroy)(plan *ego);                                /* may be NULL */
  void (*apply_adjoint)(const plan *ego, const problem *p);  /* may be NULL */
} plan_adt;
struct plan_s {
  const plan_adt *adt;
  double pcost;      /* analytic estimate; overwritten with measured cost in a race */
  int awake_state;   /* PLNR_SLEEPY | PLNR_AWAKE_ZERO | PLNR_AWAKE */
  int uses_x;        /* whether apply/apply_adjoint read x from the problem */
};

/* solver — stateless refcounted recipe */
typedef struct {
  int problem_kind;
  void (*destroy)(solver *ego);                                    /* may be NULL */
  plan *(*mkplan)(const solver *ego, const problem *p, planner *pl);
} solver_adt;
struct solver_s { const solver_adt *adt; int refcnt; };
```

**`mkplan` contract:** return a plan with `pcost` set, or `NULL` when not
applicable. It must be **cheap** and must **not** build node-dependent state —
that is `awake()`'s job. Impatience gates (`PLNR_L(pl)`) are honoured here by
returning `NULL`.

Problem kinds (current):
`enum { NFFT_PROBLEM_UNSOLVABLE, NFFT_PROBLEM_NFFT, NFFT_PROBLEM_DECONV,
NFFT_PROBLEM_CONV, NFFT_PROBLEM_LAST }`. (Note: **no** NFCT/NFST kinds — they
were removed; the DECONV/CONV kinds are the fast-NFFT decomposition children.)

## The native architecture (important)

Every candidate plan — direct or fast reads/writes exclusively through its `problem` 
(which holds the borrowed `x`/`f_hat`/`f` pointers). There is:

- **no** "legacy core", **no** `set_core`/`provisional core`/`psi_valid` coordination;
- **no** wrapping of the legacy `nfft.c` `trafo`/`adjoint` entry points.

The fast-native solver instead *decomposes* the problem into DECONV + FFTW +
CONV children (see [solvers-problems-windows.md](solvers-problems-windows.md)).
If any doc or comment mentions a shared/provisional/legacy core, it describes a
removed design — see [history-and-drift.md](history-and-drift.md).

## The search (`kernel/planner/planner.c`)

The global planner is `Y(the_planner)()` — lazily created, process-global, and
**not thread-safe** for planning (FFTW parity). Executing plans is fine.

```c
plan *Y(planner_mkplan)(planner *pl, const problem *p);   /* estimate search + memo */
int   Y(planner_candidates)(planner *pl, const problem *p,
                            plan **plans, unsigned *slvndx, int cap);
void  Y(planner_bless)(planner *pl, const problem *p, unsigned slvndx);
solution *Y(planner_hlookup)(planner *pl, const md5sig s, const flags_t *q);
void  Y(planner_hinsert)(planner *pl, const md5sig s, const flags_t *f, unsigned slvndx);
```

- **`planner_mkplan`** = the estimate path: hash the problem, consult wisdom;
  on a hit re-run the winning solver's `mkplan`; on a miss iterate
  `FORALL_SOLVERS_OF_KIND`, keep the **strictly smallest `pcost`** (ties keep the
  earlier-encountered), destroy losers immediately, and memoise the outcome
  **unblessed** (or `INFEASIBLE_SLVNDX`).
- **`planner_candidates`** = the measured path's enumerator: run every
  applicable solver of the kind under current bounds, store up to `cap` (plan,
  descriptor-index) pairs, return the count. **No** wisdom lookup, memoise, or
  store mutation — the guru drives the race and blessing itself.
- **`planner_bless`** persists a measured winner into the blessed table.

### Registration order determines tie-breaks

Solvers register push-front onto a per-kind chain (`kind_head[kind]`,
`next_same_kind`), so **iteration order is reverse registration order**. On an
exact-`pcost` tie the earlier-*encountered* (i.e. later-*registered*) plan wins.
Rosters (`kernel/nfft/conf.c`) order registrations deliberately — e.g. the
native fast solver is registered *first* so it is iterated *last* and loses
exact ties to more-proven candidates.

### The generation / ABA guard

Lazy registration keys on `Y(the_planner_generation)()`, which strictly
increases on each planner (re)creation. Registration hooks must compare
**generations, never planner pointers** — a freed-and-recreated planner can
reuse the same malloc address (ABA). This was a real, adversarially-found bug.

## Geometry: tensor / mvdim (`kernel/planner/tensor.c`)

NFFT operators are **rectangular** (output length ≠ input length), unlike FFTW's
square `iodim`. The vocabulary:

```c
typedef struct { INT n_in, is, n_out, os; } mvdim;   /* one n_out x n_in strided factor */
typedef struct { int rnk; mvdim *dims; } tensor;     /* Kronecker product of rnk factors */
```

- Square case (`n_in == n_out`) is a *state* (`Y(mvdim_square)`,
  `Y(tensor_squarep)`), not a separate type. `rnk == 0` is the scalar identity.
- The scattered-node side of NFFT is **not** a tensor factor — it stays a flat
  count `M` in the problem struct.
- `Y(tensor_adjoint)` swaps `n_in↔n_out` and `is↔os` per factor (an involution).
- `Y(tensor_compress)` / `_compress_contiguous` canonicalise (drop unit factors,
  sort by a total canonical order, merge contiguous factors) — behavior-matched
  to FFTW's `tensor_compress`. Problem constructors canonicalise at birth.
- `Y(tensor_md5)` hashes rank + every factor's `{n_in,is,n_out,os}` (strides
  included; direction-sensitive); it does **not** compress internally.

## Supporting infrastructure

- **md5** (`md5.c`) — RFC-1321 MD5, used as a fast non-cryptographic key. Feed
  bytes with `Y(md5_put_byte/bytes/str/int/INT/unsigned)`, then `Y(md5_end)`.
  `put_str` includes the terminating `'\0'` (so `"ab"+"c" ≠ "a"+"bc"`). Not
  endianness-portable (per-machine cache — a non-goal).
- **printer / scanner** (`print.c`, `printers.c`, `scan.c`, `scanners.c`) — a
  tiny locale-independent printf/scanf pair with NFFT's own directive language,
  used for plan trees and wisdom I/O. Directives include `%D` (INT), `%w` (an
  8-hex-digit wisdom word), `%p`/`%P` (plan/problem tree, `(null)` for NULL),
  `%(`/`%)` (indent/unindent). Backends: file / string / counter (for two-pass
  sizing).
- **timer** (`timer.c`) — `Y(plan_measure_cost)`; see
  [planning-modes-and-flags.md](planning-modes-and-flags.md#the-timer-kernelplannertimerc).
- **hash / primes** (`hash.c`, `primes.c`) — string hash for solver-name
  resolution; `Y(is_prime)`/`Y(next_prime)` for table sizing.

## Iteration macros

```c
FORALL_SOLVERS(pl, s, d, what)              /* every registered solver */
FORALL_SOLVERS_OF_KIND(kind, pl, s, d, what) /* one kind's chain */
```

Do **not** register solvers while one of these is iterating a raw
`slvdescs` pointer — a realloc mid-iteration is a use-after-free. This is why
cross-kind children are registered eagerly; see
[solvers-problems-windows.md](solvers-problems-windows.md#adding-a-solver-or-a-new-kind).
