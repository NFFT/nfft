"""Weighted minimax (Remez) approximation at arbitrary precision.

Both approximation forms used for I0 reduce to the same problem: find the
polynomial p of degree n that minimises

    max_{x in [a,b]} |w(x) * (f(x) - p(x))|

for a positive weight w. Choosing w = 1/f turns the weighted error into the
relative error, which is what a floating-point kernel needs.

The polynomial is carried in the Chebyshev basis of [a,b] throughout; the
monomial basis loses digits well before degree 30 at the interval widths used
here.
"""
import mpmath as mp


def cheb_basis(t, n):
    """[T_0(t), ..., T_n(t)] by the three-term recurrence."""
    out = [mp.mpf(1)]
    if n >= 1:
        out.append(t)
    for j in range(2, n + 1):
        out.append(2 * t * out[-1] - out[-2])
    return out


def clenshaw(c, t):
    """Sum_j c[j] T_j(t)."""
    n = len(c)
    if n == 1:
        return c[0]
    b1 = mp.mpf(0)
    b2 = mp.mpf(0)
    for j in range(n - 1, 0, -1):
        b1, b2 = 2 * t * b1 - b2 + c[j], b1
    return t * b1 - b2 + c[0]


def _to_unit(x, a, b):
    return (2 * x - a - b) / (b - a)


def _from_unit(t, a, b):
    return (a + b + (b - a) * t) / 2


def _solve(mat, rhs):
    """Gaussian elimination with partial pivoting on mpf matrices."""
    n = len(mat)
    m = [row[:] + [rhs[i]] for i, row in enumerate(mat)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if m[piv][col] == 0:
            raise ZeroDivisionError("singular Remez system")
        m[col], m[piv] = m[piv], m[col]
        inv = 1 / m[col][col]
        for r in range(col + 1, n):
            fac = m[r][col] * inv
            if fac == 0:
                continue
            for c in range(col, n + 1):
                m[r][c] -= fac * m[col][c]
    x = [mp.mpf(0)] * n
    for r in range(n - 1, -1, -1):
        s = m[r][n] - sum(m[r][c] * x[c] for c in range(r + 1, n))
        x[r] = s / m[r][r]
    return x


def _extrema(err, a, b, count, coarse):
    """Locate `count` alternating extrema of err() on [a,b].

    Scans a dense grid for sign-alternating local maxima of |err|, then
    refines each by ternary search. The grid is fine enough that no
    equioscillation lobe is missed for the degrees used here.
    """
    xs = [_from_unit(-mp.cos(mp.pi * i / coarse), a, b) for i in range(coarse + 1)]
    vs = [err(x) for x in xs]

    # Candidate local extrema of |err|. An endpoint qualifies only if |err|
    # actually rises towards it; where the weight vanishes (branch 1 at y=0)
    # the endpoint is a zero of the error, not a lobe, and including it would
    # pin the equioscillation level to zero.
    cand = []
    if abs(vs[0]) >= abs(vs[1]):
        cand.append(0)
    for i in range(1, coarse):
        if abs(vs[i]) >= abs(vs[i - 1]) and abs(vs[i]) >= abs(vs[i + 1]):
            cand.append(i)
    if abs(vs[coarse]) >= abs(vs[coarse - 1]):
        cand.append(coarse)

    refined = []
    for i in cand:
        lo = xs[max(i - 1, 0)]
        hi = xs[min(i + 1, coarse)]
        sign = 1 if vs[i] >= 0 else -1
        x = _refine(lambda z: sign * err(z), lo, hi)
        v = err(x)
        # Keep the grid point if the refinement wandered downhill.
        if abs(v) < abs(vs[i]):
            x, v = xs[i], vs[i]
        refined.append((x, v))

    # Keep one extremum per sign run, the largest in each run.
    runs = []
    for x, v in refined:
        s = 1 if v >= 0 else -1
        if runs and runs[-1][0] == s:
            if abs(v) > abs(runs[-1][2]):
                runs[-1] = (s, x, v)
        else:
            runs.append((s, x, v))

    # Discard degenerate lobes: a run whose peak is orders below the largest is
    # not an alternation point, it is a zero crossing that survived rounding.
    if not runs:
        return None
    top = max(abs(r[2]) for r in runs)
    if top == 0:
        return None
    runs = [r for r in runs if abs(r[2]) > top * mp.mpf(10) ** -6]

    if len(runs) < count:
        return None
    # Drop the weakest end lobes until the reference has exactly `count` points.
    while len(runs) > count:
        if abs(runs[0][2]) <= abs(runs[-1][2]):
            runs.pop(0)
        else:
            runs.pop()
    return [r[1] for r in runs]


def _refine(g, lo, hi, rounds=6):
    """Maximise g on [lo, hi] by repeated three-point parabolic fitting.

    The error curve is quadratic-flat at an equioscillation peak, so the peak
    *location* barely matters; a handful of rounds pins it far more cheaply
    than a long ternary search, and the objective is evaluated ~20 times
    instead of ~160.
    """
    for _ in range(rounds):
        m = (lo + hi) / 2
        gl, gm, gh = g(lo), g(m), g(hi)
        denom = (gl - 2 * gm + gh)
        if denom != 0:
            step = (hi - lo) * (gl - gh) / (4 * denom)
            cand = m + step
            if lo < cand < hi and g(cand) > gm:
                m, gm = cand, g(cand)
        # Shrink towards the current best.
        q = (hi - lo) / 4
        lo, hi = max(lo, m - q), min(hi, m + q)
    return (lo + hi) / 2


def _memo(fn):
    """Cache fn on exact mpf arguments.

    The grid _extrema scans is identical on every exchange step, so f and w are
    otherwise recomputed from scratch each iteration; besseli at 60 digits is
    by far the dominant cost.
    """
    cache = {}

    def wrapped(x):
        key = mp.nstr(x, 50)
        v = cache.get(key)
        if v is None:
            v = fn(x)
            cache[key] = v
        return v

    return wrapped


def minimax_poly(f, w, a, b, n, iters=30, tol=mp.mpf(10) ** -25, coarse=None):
    """Degree-n weighted minimax polynomial on [a, b].

    f, w: callables taking and returning mpf. w must be >= 0 and positive
    except possibly at an endpoint.

    Returns (cheb_coeffs, err) where the polynomial is
    sum_j cheb_coeffs[j] * T_j((2x-a-b)/(b-a)) and err is the equioscillation
    level of w*(f-p).
    """
    a = mp.mpf(a)
    b = mp.mpf(b)
    f = _memo(f)
    w = _memo(w)
    if coarse is None:
        coarse = max(10 * (n + 2), 160)

    # A weight that vanishes at an endpoint (branch 1 has w(0) = 0) would make
    # that reference row read 0 = +-E and pin the level to zero. The minimax
    # problem is unchanged by excluding a point of zero weight, so the search
    # steps just inside the interval -- but the Chebyshev basis stays on the
    # original [a, b]. Shifting the basis too would leave the returned
    # coefficients describing a slightly different interval than the caller
    # asked for, an error floor of the order of the step itself.
    slo, shi = a, b
    if w(a) == 0:
        slo = a + (b - a) * mp.mpf(10) ** -30
    if w(b) == 0:
        shi = b - (b - a) * mp.mpf(10) ** -30

    # Initial reference: Chebyshev extrema of the search interval.
    ref = [_from_unit(-mp.cos(mp.pi * i / (n + 1)), slo, shi) for i in range(n + 2)]

    def measured(c, samples=None):
        """Max weighted error by dense scan, independent of the solver level.

        The equioscillation level E that falls out of the linear system is only
        trustworthy once the reference has converged; every error this module
        reports is this measured value instead.
        """
        m = samples or max(40 * (n + 2), 800)
        worst = mp.mpf(0)
        for i in range(m + 1):
            x = _from_unit(-mp.cos(mp.pi * i / m), slo, shi)
            worst = max(worst, abs(w(x) * (f(x) - clenshaw(c, _to_unit(x, a, b)))))
        return worst

    c = None
    best = None
    best_err = None
    for _ in range(iters):
        # Solve w(x_i)(f(x_i) - sum_j c_j T_j(t_i)) = (-1)^i E.
        mat = []
        rhs = []
        for i, x in enumerate(ref):
            wi = w(x)
            row = [wi * t for t in cheb_basis(_to_unit(x, a, b), n)]
            row.append(-mp.mpf(-1) ** i)
            mat.append(row)
            rhs.append(wi * f(x))
        sol = _solve(mat, rhs)
        c = sol[:n + 1]

        def errfun(x, c=c):
            return w(x) * (f(x) - clenshaw(c, _to_unit(x, a, b)))

        new = _extrema(errfun, slo, shi, n + 2, coarse)
        if new is None:
            break
        peak = max(abs(errfun(x)) for x in new)
        lo = min(abs(errfun(x)) for x in new)
        cur = peak
        if best_err is None or cur < best_err:
            best, best_err = c, cur
        ref = new
        if peak == 0 or (peak - lo) <= tol * peak:
            break

    if best is None:
        best = c
    # Always report an independently measured error, never the solver level.
    return best, measured(best)


def cheb_series(f, a, b, n):
    """Plain degree-n Chebyshev interpolant of f on [a, b] (no weighting).

    Used as the starting point and as a cross-check on the Remez result.
    """
    a = mp.mpf(a)
    b = mp.mpf(b)
    m = n + 1
    nodes = [mp.cos(mp.pi * (k + mp.mpf(1) / 2) / m) for k in range(m)]
    vals = [f(_from_unit(t, a, b)) for t in nodes]
    c = []
    for j in range(m):
        s = sum(vals[k] * mp.cos(mp.pi * j * (k + mp.mpf(1) / 2) / m)
                for k in range(m))
        c.append(2 * s / m)
    c[0] /= 2
    return c


def cheb_to_monomial(c, a, b):
    """Convert Chebyshev coefficients on [a,b] to monomial coefficients in x.

    Only used for the `1 + y*P(y)` branch, whose Horner form in the raw
    variable is both cheaper and better conditioned than Clenshaw (all
    coefficients there are positive, so Horner cannot cancel).
    """
    n = len(c) - 1
    # T_j in the variable t, then t = (2x - a - b)/(b - a).
    tpoly = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]]
    for j in range(2, n + 1):
        prev = tpoly[j - 1]
        shifted = [mp.mpf(0)] + [2 * v for v in prev]
        cur = shifted[:]
        for i, v in enumerate(tpoly[j - 2]):
            cur[i] -= v
        tpoly.append(cur)

    # p(t) = sum_j c_j T_j(t)
    pt = [mp.mpf(0)] * (n + 1)
    for j in range(n + 1):
        for i, v in enumerate(tpoly[j]):
            pt[i] += c[j] * v

    # substitute t = alpha*x + beta
    alpha = 2 / (b - a)
    beta = -(a + b) / (b - a)
    # Horner in x over the t-polynomial: acc = acc*(alpha x + beta) + pt[k]
    acc = [mp.mpf(0)] * (n + 1)
    for k in range(n, -1, -1):
        new = [mp.mpf(0)] * (n + 1)
        for i, v in enumerate(acc):
            if v == 0:
                continue
            new[i] += v * beta
            if i + 1 <= n:
                new[i + 1] += v * alpha
        new[0] += pt[k]
        acc = new
    return acc
