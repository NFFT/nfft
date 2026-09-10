# ADR-0005: Runtime-dispatched SIMD kernels for the NFFT convolution

## Status
Accepted; implemented for the plain NFFT module (`kernel/nfft`).

## Context
The cost of an NFFT that is not FFT-bound sits in the convolution step: for each
node, `f_j` gathers a `(2m+2)^d` window of the oversampled grid `g` (forward),
or scatters into it (adjoint). Along the last dimension the grid points a node
touches are contiguous, and the window weights `psi` are **real** while `g` is
complex. Each innermost loop is therefore one of two kernels:

    f += sum_l psi[l] * g[l]        (forward)
    g[l] += psi[l] * f              (adjoint)

Both vectorise trivially once you look at `g` as `2*len` reals and duplicate
each weight across an adjacent lane pair. Compilers do reach some of this on
their own, but only for the instruction set the *build* targets — and NFFT3 is
overwhelmingly consumed as a distribution package built for the baseline
architecture, where that means SSE2 on x86-64 and nothing wider.

Three constraints shape the design:

- The library is precision-agnostic and builds as float, double and long
  double from the same sources. Long double is padded (10 bytes of value in 12
  or 16), so the "complex is two adjacent reals" identity the kernels rest on
  does not hold there.
- NFFT3 builds with GCC, clang and (via MinGW) on Windows, and is expected to
  build with less capable compilers too. Nothing may become mandatory.
- A binary built on one machine must keep running on another. Selecting an
  instruction set at build time and hoping the host has it turns a portability
  question into a SIGILL.

## Decision
Compile **every** instruction set the compiler can express into the library,
side by side with the plain C kernels, and pick between them **at run time**
against the host CPU.

The split is:

- `include/simd.h` — what this build carries (`NFFT_SIMD_HAVE_<isa>`,
  `NFFT_SIMD_MAX`) and what the host supports (`Y(simd_isa)()`, resolved once
  via `__builtin_cpu_supports` and cached; `kernel/util/simd.c`).
- `include/simd_run.h` — the two primitives above, one instantiation per
  instruction set. No include guard: it is included once per variant with
  `NFFT_SIMD_VARIANT` set, and names its output `nfft_simd_cdot_avx2` and so on.
- `kernel/nfft/compute_body.h` — the per-node 1D/2D/3D forward and adjoint
  kernels, written once against those primitives and likewise compiled once per
  variant.
- `kernel/nfft/nfft.c` — a table of function pointers, filled in
  `init_help()`, i.e. once per plan. Where only the scalar variant exists the
  table disappears and the call sites bind directly.

Widening beyond the build's own target rests on GCC's and clang's function
`target` attribute, which lets one object file hold functions compiled for
different instruction sets. That keeps the whole mechanism inside the C
sources: no per-file `-m` flags, no separate objects, and nothing for either
build system to do beyond a single on/off switch (`--disable-simd`,
`-DNFFT_ENABLE_SIMD=OFF`).

Currently compiled: SSE2, AVX and AVX2+FMA on x86/x86-64; Advanced SIMD on
AArch64, where it is baseline and needs neither attribute nor run-time check.
Everything else — other architectures, other compilers, long double,
`--disable-simd` — gets the scalar kernels, which are plain C and are what the
transforms computed before this change.

`NFFT_SIMD=scalar|sse2|avx|avx2|neon` in the environment overrides the
detection, clamped to what the build has; `Y(simd_force_isa)()` does the same
from code. Both exist so that the same binary can be run down the whole ladder
and compared, which is what `tests/simd.c` does.

### Expressing the runs
Getting the innermost loop into a single contiguous call meant reshaping the
kernels. A node's `2m+2` point run wraps around the periodic grid at most once,
so it is one or two contiguous pieces (`nfft_run`); the previous code spelled
out the wrap/no-wrap combinations as a `2^d`-way branch tree with the loop body
repeated in each arm. The runs now carry that structure as data, which collapses
those trees and hoists the outer dimensions' window factors out of the innermost
loop instead of re-multiplying them per grid point. The OpenMP blockwise adjoint
kernels lose their `index_temp` modulo tables to the same descriptor.

This changes the *association* of the sums, not their terms: results move by a
few ulps and, in the 2D and 3D cases, generally toward the more accurate answer,
since the per-dimension partial sums are now accumulated separately.

## Consequences
- The plain C path is never removed and never optional. Every host that could
  run NFFT3 before still runs it, taking the scalar kernels if nothing else
  applies.
- Two kernels, six shapes, one source each. A new instruction set is a block of
  eight macros plus two small helpers in `simd_run.h`; nothing in
  `compute_body.h` or `nfft.c` changes.
- One indirect call per *node*, not per grid point — a single-target call the
  branch predictor resolves, amortised over the whole window.
- Cross-variant agreement is a unit test (`tests/simd.c`), not a review
  argument: every instruction set available on the test host runs the same 1D,
  2D and 3D transforms and has to match the scalar result to rounding.

## Verification evidence
`make check` and `ctest` pass with 0 failures across the double, float,
long-double, `--disable-simd` and OpenMP builds, and under CMake.

Instruction counts (callgrind `I refs`, one forward plus one adjoint transform,
`PRE_PSI`, double precision), `-O3 -ffast-math`:

| build target | case | scalar | AVX2 |
|---|---|---|---|
| baseline x86-64 | 3D, N=32³, M=8192, m=4 | 254.3 M | 177.5 M (−30%) |
| baseline x86-64 | 2D, N=128², M=16384, m=6 | 145.0 M | 113.4 M (−22%) |
| `-march=haswell` | 3D, N=32³, M=8192, m=4 | 226.7 M | 173.3 M (−24%) |
| `-march=haswell` | 2D, N=128², M=16384, m=6 | 125.3 M | 109.1 M (−13%) |

The `-march=haswell` rows are the honest comparison: there the scalar loops are
already auto-vectorised to AVX2 by GCC, and the explicit kernels still win.
1D gains little (2–5%), as expected — a 1D NFFT of this size is FFT-bound.

## Alternatives considered
- **Build-time selection only** (`--enable-avx2` adding `-mavx2` globally, the
  classic FFTW-style knob). Simpler, but it makes the binary non-portable and
  leaves the common distribution build — the one that most needs the speedup —
  on the baseline instruction set.
- **`__attribute__((target_clones))` / ifunc.** Would dispatch for free, but
  clang's support is newer and macOS has no ifunc, so it is not available
  everywhere NFFT3 builds.
- **Dispatching per run rather than per node.** An indirect call every `2m+2`
  points is a real cost at these run lengths; per node it is not.
- **Relying on the auto-vectoriser alone.** It cannot emit code for an
  instruction set the command line did not enable, which is the whole problem.
