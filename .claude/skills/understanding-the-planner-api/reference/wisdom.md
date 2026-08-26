# Wisdom: the size-class key & the cache

Wisdom is the planner's **memo**: for a problem *size class* (a 128-bit MD5
signature), which solver won a search conducted under given impatience bounds.
It is a **cache, never a source of truth** — any mismatch on lookup or import
degrades to a miss, never a wrong plan. Ground truth: `kernel/planner/planner.c`
(store), `kernel/planner/md5.c`, `kernel/nfft/api_ng.c` (public wrappers),
`include/iplanner.h`.

## What the key is (the "size class")

`Y(problem_md5)(pl, problem, out)` builds the key. It folds in, structurally
(never node coordinates):

- `sizeof(R)` (precision) and `pl->nthr` (thread count) — added by
  `problem_md5` itself, so concrete problem hashes must not repeat them.
- the problem **kind**;
- the frequency-tensor geometry `sz` and batch geometry `vecsz` via
  `Y(tensor_md5)` — **strides included** (layout changes which solver wins);
- `M` **bucketed** to `floor(log2 M)` (node count in coarse power-of-two
  buckets; everything else exact);
- window cutoff `m`, the **runtime window ordinal**, per-axis NDFT `variant`;
- `sign` (direction — retained even though only forward is raced today);
- the **keyable** `fftw_flags` (planning-rigor bits only; preservation bits and
  the removed `nfft_flags` do not participate).

Consequences you can rely on:

- Unit axes are elided at construction, so a top-level `{64,64,1}` problem hashes
  **identically** to `{64,64}`.
- `M = 900` and `M = 1000` share a bucket (`floor(log2 M) = 9`) → same key.
- Node coordinates never enter the key; wisdom means "for this *size class* on
  *this machine*, solver X won." Accepting that node distribution can
  occasionally make the cached choice suboptimal is deliberate.

## Blessed vs unblessed

Two hash tables in `struct planner_s`: `htab_blessed` and `htab_unblessed`.

- **Blessed** (`PLNR_BLESSING` set): worth persisting. Exported to wisdom files,
  survives `PLNR_FORGET_UNBLESSED`. Only *measured* winners are blessed.
- **Unblessed**: session-only search memoisation (estimate memos, stale-search
  outcomes). Dropped on forget; never exported.

Infeasibility (no applicable solver under the bounds) is recorded with the
sentinel `INFEASIBLE_SLVNDX` (`0xFFFF`) so a repeated query short-circuits.

The store is open-addressed with double hashing, prime table sizes, and
grow-before-load-8/9 rehashing (drops dead/invalid slots). `flags_t` packs
`{l:16, info:3, timelimit_imp:12, u:16, slvndx:16}` into 64 bits.

## Export / import

```c
int   nfft_export_wisdom_to_filename(const char *filename);
char *nfft_export_wisdom_to_string(void);        /* caller frees with nfft_free */
int   nfft_import_wisdom_from_filename(const char *filename);
int   nfft_import_wisdom_from_string(const char *s);
void  nfft_forget_wisdom(void);                  /* drops blessed AND unblessed */
```

All call `Y(nfft_ensure_registered)()` first so solver names resolve. Export is
two-pass (a counter printer sizes the buffer, a string printer fills it). Import
returns `1` on success, `0` on any failure, and is **atomic**: `planner_import`
snapshots the blessed table (deep-copying its entries) and restores it on any
malformed input, so a failed import never corrupts the live store.

### The wisdom file grammar (NFFT's own design)

Shape (verify exact spelling in `kernel/planner/printers.c`/`scanners.c` if you
touch it):

```
(nfft_wisdom-<PACKAGE_VERSION> #x<c0> #x<c1> #x<c2> #x<c3>
  (<registrar-name> <reg-id> #x<l> #x<u> #x<tli> #x<s0> #x<s1> #x<s2> #x<s3>)
  ...
)
```

- The leading token binds precision + release: `nfft_wisdom-` for double,
  `nfftf_wisdom-` / `nfftl_wisdom-` for float/long-double, then the exact
  `PACKAGE_VERSION`. **Wisdom never crosses precision or release** — a mismatch
  is a clean rejection.
- `#x<c0..c3>` is the **configuration signature** (see below).
- One line per live blessed entry, each identifying its winning solver by
  `(registrar-name, reg-id)` plus its bounds and 128-bit key words. An
  infeasible entry uses the reserved name `!` with id `0`.

### The configuration signature (the roster fingerprint)

An MD5 over `sizeof(R)` followed by, in registration order, **every** registered
solver's `reg_id` and registrar name — computed across *all* kinds, not
per-kind. Written on export, checked on import **before** the store is touched.

Because it covers the whole roster, **adding or removing any solver changes the
signature**, so wisdom produced by a process with a different solver set (a
different precision, or a build that registers extra kinds) fails import — safe
by design. This is why re-adding NFCT/NFST solvers later would invalidate old
wisdom (harmless: it is a cache).

## Three import guards (defense in depth)

1. **Preamble** — version + precision token must match exactly.
2. **Configuration signature** — the whole-roster MD5 must match.
3. **Per-entry solver resolution** — each `(name, id)` must resolve to a
   registered solver (via `nam_hash` prefilter + `strcmp` + `reg_id`);
   unresolvable → treated as malformed. Feasible entries must have
   `timelimit_imp == 0`; every entry must satisfy `LEQ(l, u)`.

Any failure at any guard → clean rejection / cache miss, never a wrong plan.
