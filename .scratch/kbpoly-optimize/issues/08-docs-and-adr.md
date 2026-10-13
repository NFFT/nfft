# 08: Documentation, glossary and ADR

Status: ready-for-agent

Part of: PRD C3. Depends on: issues 02 to 05 for the facts it records.

## Task

- `docs/adr/0005-fit-the-polynomial-window-to-the-transform-error.md`:
  context (the function-fit and its table rule), decision (strict degree
  criterion from the error functional; optional least-squares and Lawson
  correction fitted to the residual in `E`; guard in the table), consequences
  (window is stencil-specific; `tests/window.c`-style comparisons against the
  Kaiser-Bessel function no longer describe the window the transform uses;
  gains are constant factors, the rate is unchanged), alternatives rejected
  (free deconvolution, coefficient table, continuous optimal kernels), with
  the measurements of `analysis.md` sections 6, 9 and 10 summarized.
- `CONTEXT.md`: add glossary entries, with the "avoid" lines the file uses:
  **error functional** (`E(xi, t) = D/C - 1`, section 3), **polynomial part**
  (the polynomial's contribution to the transform error, section 9),
  **correction step** (the optional plan-time fit of the difference to the
  Chebyshev interpolant), **regime guard** (the table threshold that turns the
  step off near the rounding floor).
- `AGENTS.md` section 5 and the quick-reference table: the generator's new
  dependency and command; `tests/kbpolygen/README.md`: the new rule, the
  guard, the certification.
- `include/nfft3.h`: the flag's doc comment; `include/infft.h`: the comment
  block above `Y(kb_poly_degree)` and `Y(kb_poly_fit)`.

## Acceptance

- The ADR is numbered after the last existing one and linked from
  `analysis.md`'s status line.
- `CONTEXT.md` terms are used verbatim in the issues' and the code's comments
  (grep for drift: "optimizer", "optimized window", "LUT" are the synonyms to
  avoid).
