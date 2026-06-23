# Floating-point error analysis of the fast 1D NFFT (forward & adjoint)

Status: analysis (companion to PR #222, which covers the *direct* NDFT path)

This note analyses the **fast** one-dimensional transforms `nfft_trafo_1d` and
`nfft_adjoint_1d` in `kernel/nfft/nfft.c`. PR #222 ("NDFT blocked recurrence")
showed that the *direct* transforms accumulate rounding error that grows like
`O(√N)` (forward) / `O(N)` (adjoint), and fixed it with FMA-based argument
reduction. The natural question for the fast path is whether it has an analogous
problem. The short answer is **no** — and this note explains why, quantifies the
real rounding sources, validates the result empirically, and lists concrete
accuracy and performance improvements.

---

## 1. The three-step algorithm

The fast NFFT factorises the nonequispaced Fourier matrix `A ≈ B F D`:

Forward `nfft_trafo_1d` (`nfft.c:2779`):

1. **D — deconvolution** (`nfft.c:2808`): `ĝ_k = f̂_k · c_k`, where
   `c_k = 1/φ̂(k)` (precomputed reciprocal of the window's Fourier coefficient,
   `PRE_PHI_HUT`), written into the oversampled array `g_hat` of length `n = σN`
   in FFT layout (negative frequencies at the top).
2. **F — FFT** (`nfft.c:2838`): `g = FFT_n(ĝ)` via FFTW (`my_fftw_plan1`).
3. **B — convolution / gather** (`nfft.c:2842`, `nfft_trafo_1d_compute`,
   `nfft.c:2135`): for each node, `f_j = Σ_{l=0}^{2m+1} ψ_{j,l} · g_{(u_j+l) mod n}`,
   a dot product of length `2m+2` with window samples `ψ_{j,l} = φ(x_j-(u_j+l)/n)`.

Adjoint `nfft_adjoint_1d` (`nfft.c:2847`) runs the transpose in reverse order:
**Bᵀ — spread** (`nfft_adjoint_1d_B`, `nfft.c:2589`) → **Fᴴ — inverse FFT**
(`my_fftw_plan2`) → **Dᴴ — deconvolution** (`nfft.c:2878`).

The total error of the computed result splits as

```
   f̃ − f_exact  =  E_approx  +  E_round
```

* `E_approx` — the factorisation `B F D ≠ A` is only approximate. This is the
  classical NFFT truncation/aliasing error, **independent of N**, decaying like
  `e^{-2πm√(1−1/σ)}` for Kaiser–Bessel (exactly the bound `err_trafo` checks in
  `tests/nfft.c:325`). It is governed by `m` and `σ` and is the deliberate
  accuracy knob.
* `E_round` — finite-precision evaluation of `D`, `F`, `B`. This is what the
  analysis below is about.

Throughout, `u = 2⁻⁵³ ≈ 1.1·10⁻¹⁶` is the double-precision unit roundoff and the
standard model `fl(a∘b) = (a∘b)(1+δ), |δ|≤u` is used.

---

## 2. Rounding analysis — forward

### Step D (deconvolution) — `O(u)`, elementwise

`ĝ_k = fl(f̂_k · c_k)` with `c_k = fl(1/φ̂(k))` precomputed
(`precompute_phi_hut`, `nfft.c:5758`). Each output carries the relative error of
the stored reciprocal (window evaluation `bessel_i0`/`sqrt`, a few ulp) plus one
multiply:

```
   ĝ_k = (f̂_k / φ̂(k)) (1 + e_k),   |e_k| ≲ (ρ + 2) u
```

with `ρ` the few-ulp window-evaluation constant. **No summation, no N-dependence.**

Conditioning note: `|c_k| = 1/φ̂(k)` grows toward the band edges `|k|→N/2` by a
factor `κ_D = φ̂(0)/φ̂(N/2)` (modest for `σ=2`, larger as `σ→1`). This amplifies
high-frequency content of `f̂` — but `κ_D` depends only on `σ, m`, **not on N**,
so it does not introduce N-growth. It is the mathematically-required deconvolution,
not a numerical defect.

### Step F (FFT) — `O(u·log n)`, the only n-dependent term

FFTW is backward stable; the standard bound (Higham, *ASNA* §24) is

```
   ‖g̃ − FFT_n(ĝ)‖₂  ≤  c · u · log₂(n) · ‖ĝ‖₂ ,   c small.
```

So the FFT contributes a **normwise** relative rounding of size `O(u·log n)` —
logarithmic in the problem size and irreducible (internal to FFTW). With
`n = σN`, this is the closest thing to N-dependence in the whole pipeline, and it
is only `log N`.

### Step B (gather) — `O(m·u)`, short dot product

`f_j` is a dot product of length `L = 2m+2`. With precomputed window samples
(`PRE_PSI`) the dot-product bound (Higham §3.1) gives

```
   |f̃_j − f_j|  ≤  γ_L · Σ_l |ψ_{j,l}|·|g_{·}| ,   γ_L = Lu/(1−Lu) ≈ (2m+2)u.
```

`L = 2m+2` is a **small constant** (`m ∈ [2,12]`), and the window samples are
non-negative so there is no cancellation among the `2m+2` terms (local condition
number `Σ|ψ g| / |Σ ψ g| = O(1)`). Hence **`O(m·u)`, N-independent.**

A subtlety worth recording — the window **argument** `x_j − (u_j+l)/n` is formed
by an explicit subtraction (`nfft.c:2444`, `precompute_psi` `nfft.c:5844`).
Naively this is a difference of two `O(½)` quantities producing an `O(m/n)`
result — apparent cancellation. It is saved by the fact that **`n` is a power of
two** for every default plan (`n = 2·next_power_of_2(N)`, `nfft.c:6068`):

* `(u_j+l)/n` is an integer divided by a power of two ⇒ **exact**;
* `x_j − (u_j+l)/n` is then exact by Sterbenz (operands within a factor 2);
* the window's internal `×n` (`(x)*(R)n` in the Kaiser–Bessel `PHI`,
  `infft.h:259`) multiplies by a power of two ⇒ **exact**.

So for power-of-two `n` the scaled window argument carries **no** rounding error,
and the only B-step error is the few-ulp window evaluation. (See §5.1 for the
non-power-of-two case, which is *not* protected and is the one genuine accuracy
recommendation.)

### Forward total

```
   E_round(forward)  ≲  u · ( (ρ+2)      [D]
                            +  c·log₂ n   [F]
                            + (2m+2)·κ    [B] )
                     =  O( u · (m + log n) )      — N-independent up to log N.
```

---

## 3. Rounding analysis — adjoint

The adjoint reverses the order; the per-step bounds are the transpose of the above
with one important difference in the spread step.

* **Dᴴ** (`nfft.c:2878`): elementwise, `O(u)`. Same as D.
* **Fᴴ** (`nfft.c:2875`): `O(u·log n)`. Same as F.
* **Bᴴ — spread** (`nfft_adjoint_1d_compute_serial`, `nfft.c:2160`): each grid
  point `g_l` *accumulates* `g_l += ψ · f_j` over **every node whose stencil
  covers `l`**. The accumulation length is therefore the per-cell node
  multiplicity `p_l`, not a fixed `2m+2`:
  ```
     |g̃_l − g_l| ≤ γ_{p_l} · Σ |ψ·f| ,   p_l = #{ j : l ∈ supp ψ_j }.
  ```
  For quasi-uniform nodes with `M=N`, `p_l ≈ (2m+1)/σ = O(m)` ⇒ `O(m·u)`,
  N-independent — the common case. **But under node clustering `p_l` can reach
  `O(M)=O(N)`**, and the adjoint spread then degrades to `O(N·u)`. This is the
  one place where the fast adjoint can, in principle, inherit N-growth — though
  driven by node distribution, not by `N` per se. (The OpenMP atomic/blockwise
  variants change the summation *order*, not this bound.)

### Adjoint total

```
   E_round(adjoint)  ≲  O( u · (p̄ + log n) ),   p̄ = max_l p_l.
```

For uniform nodes `p̄ = O(m)` ⇒ N-independent; pathological clustering can lift it.

---

## 4. Main result and empirical validation

> **The fast 1D forward and adjoint transforms have rounding error
> `O(u·(m + log n))` — essentially independent of `N` (logarithmic at worst).
> The `O(eps)`, N-independent error bound the tests use (`err_trafo`,
> `tests/nfft.c:219`) is therefore justified for the fast path, which is exactly
> why the fast path needed no analogue of the PR #222 fix.**

Validation experiment (`.scratch/fast-nfft-1d/fp_experiment.c`): random `f̂`/`f`, random
uniform nodes, `m=12` (so `E_approx ≈ e^{−2π·12·√½} ≈ 10⁻²² ≪ E_round`,
isolating rounding), `σ=2`, relative ℓ₂ error against an **independent
long-double reference** with exact phase reduction. Errors of the fast transforms
vs. the library's own direct transforms:

```
        N        fast_fwd      direct_fwd        fast_adj      direct_adj
        16      1.0524e-14      7.2706e-16      7.6050e-15      4.1005e-16
        64      1.0318e-14      2.6082e-15      9.4946e-15      2.9728e-15
       256      1.1247e-14      1.1280e-14      1.1581e-14      1.0522e-14
      1024      1.0438e-14      4.1841e-14      1.0603e-14      4.3670e-14
      4096      1.0587e-14      1.6509e-13      1.0705e-14      1.6723e-13
      8192      1.0341e-14      3.3422e-13      1.0291e-14      3.3009e-13
```

* `fast_fwd`, `fast_adj`: **flat at ≈1.0·10⁻¹⁴** across three decades of `N`
  (no trend — the `O(m·u)` window/deconvolution floor dominates; even the `log n`
  term, a 2.8× span here, is not visible). `≈1.0·10⁻¹⁴ ≈ 95u`, consistent with a
  modest `(m + log n)` constant at `m=12, n≤2¹⁴`.
* `direct_fwd`, `direct_adj`: **grow with `N`** (≈ linearly here for random data),
  crossing above the fast error near `N≈256`. This is the pre-PR-#222 behaviour.

Two takeaways:
1. The fast transform is the numerically *robust* path at scale; for large `N` it
   is **more accurate** than the unpatched direct transform.
2. For *small* `N` (≤128) the direct transform is more accurate, because its sum
   is short and its arguments small, while the fast path pays a fixed `~10⁻¹⁴`
   FFT+window floor regardless of `N`. (The library already dispatches to
   `trafo_direct` for `N ≤ m` or `n ≤ 2m+2`, `nfft.c:2781`.)

---

## 5. Recommendations

### Accuracy

**5.1 — Robust window argument for non-power-of-two `n` (the one real gap).**
The exactness argument in §2 relies on `n` being a power of two. `nfft_init`
always picks such `n`, but `nfft_init_guru` lets callers pass an arbitrary `n`
(e.g. `σ=1.25` to save memory). For non-power-of-two `n`:

* `(u_j+l)/n` is no longer exact, and `x_j−(u_j+l)/n` then loses absolute
  precision `~u`, which the window's internal `×n` amplifies to a scaled-argument
  error `~n·u`. Because the transform is an `O(N)`-bandlimited function of the
  node, this surfaces as an `O(N·u)` *relative* error in `f` — i.e. the fast path
  would acquire exactly the kind of N-growth PR #222 removed from the direct one.

  Fix: form the scaled argument directly and reuse the product `nfft.c`'s
  `uo`/`uo2` already compute. Let `w = x_j·n` (computed once; exact when `n` is a
  power of two, otherwise the single unavoidable rounding), then the per-term
  scaled argument is `s_l = w − (u_j+l)` (a small integer-offset from `w`, formed
  with no catastrophic cancellation). Pass `s_l` to a window variant that does
  **not** re-multiply by `n`. This is exact for power-of-two `n` (no regression)
  and removes the `~n·u` argument error for general `n`. It also saves a division
  and a multiply per stencil tap (see 5.4). The change touches the shared `PHI`
  macro interface (`infft.h`) and all `*_B` call sites, so it is the larger of the
  recommendations and should be scoped carefully to preserve the precision-agnostic
  build.

**5.2 — Adjoint spread under clustering (situational).** If clustered-node
accuracy matters, a 2Sum/Kahan-compensated accumulation in
`nfft_adjoint_1d_compute_serial` (`nfft.c:2160`) caps the spread error at `O(u)`
regardless of `p̄`. Skip it for the uniform-node common case — the `2m+2`-length
forward gather and the typical `O(m)` spread are too short to benefit.

**5.3 — FG recurrence (minor).** `PRE_FG_PSI`/`FG_PSI` build the `2m+2` window
samples by the recurrence `fg_psij2 *= fg_psij1` (`nfft.c:2346`), accumulating
`O(m·u)` over the stencil. Already `O(m·u)`; a centred recurrence from the stencil
midpoint would roughly halve the constant. Low priority.

### Efficiency

**5.4 — Avoid the per-tap division (pairs with 5.1).** The no-`PRE_PSI` and
precompute paths divide by `n` per stencil tap to form `x_j−(u+l)/n`, and the
window then multiplies by `n` again. Passing the scaled argument `s_l = w−(u+l)`
removes both, replacing `2m+2` divisions + multiplies per node by additions.

**5.5 — Zero only the oversampling gap in step D (implemented).** `trafo_1d`
zeroed all `n` entries of `g_hat` and then overwrote the two outer blocks
(`[0,N/2)` and `[n−N/2,n)`, total length `N`) in the deconvolution loop. Only the
middle gap `[N/2, n−N/2)` (length `n−N`) is genuine zero-padding. Zeroing just the
gap removes `N` redundant complex writes (for `σ=2`, halving the zeroing traffic).
Safe because `σ>1` is enforced (`nfft.c:6199`), so `n>N` and the two written
blocks never overlap — the same precondition the existing D-loop indexing already
assumes. Verified numerically identical (bit-for-bit) on the experiment above.
The adjoint zeroing (`nfft_adjoint_1d_B`, `nfft.c:2595`) must stay full: the
spread touches only a node-dependent subset and the inverse FFT reads all `n`.

**5.6 — FFT dominates; expose `FFTW_MEASURE`.** Step F is the asymptotic cost
(`O(n log n)` vs `O(N)+O(Mm)` for D+B). Plans use `FFTW_ESTIMATE` by default
(`nfft.c:6087`); `FFTW_MEASURE` is already user-selectable via `fftw_flags` and is
the highest-leverage performance knob for repeated transforms.

---

## 6. Implemented in this branch

Only the safe, self-contained 5.5 is applied to `kernel/nfft/nfft.c` (`trafo_1d`,
both the OpenMP and serial zeroing). 5.1/5.4 are left as a scoped follow-up
because they touch the shared window-macro interface.

## 7. Reproduction

```bash
sudo apt-get install -y libfftw3-dev libcunit1-dev
./bootstrap.sh && ./configure --enable-all --enable-tests && make -j
gcc -O2 -I include .scratch/fast-nfft-1d/fp_experiment.c .libs/libnfft3.a -lfftw3 -lm -o /tmp/fp_exp
/tmp/fp_exp          # prints the table in §4
make check           # CUnit suite (tests/checkall): 1600 OK / 0 FAIL with this change
```

`.scratch/fast-nfft-1d/fp_experiment.c` is the standalone harness; `m` and the
`N` range are set at the top of `main`.
