/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Primal-dual interior-point method for min q^T x subject to h - G x in a
 * product of three-dimensional second-order cones Q^3 = {u : u0 >= |(u1,
 * u2)|}: Nesterov-Todd scaling, Mehrotra predictor-corrector, normal
 * equations by dense Cholesky. */

#include "nfft3.h"
#include "infft.h"

static E e_abs(const E a)
{
  return IF(a < EK(0.0), -a, a);
}

static E dot(const E *a, const E *b, const INT n)
{
  E r = EK(0.0);
  INT k;
  for (k = 0; k < n; k++)
    r += a[k] * b[k];
  return r;
}

static int all_finite(const E *v, const INT n)
{
  INT k;
  for (k = 0; k < n; k++)
    if (!Y(e_finite)(v[k]))
      return 0;
  return 1;
}

/* Smallest u0 - |(u1, u2)| over the cones. */
static E cone_margin(const E *u, const INT M)
{
  E r = u[0] - E_SQRT(u[1] * u[1] + u[2] * u[2]);
  INT i;
  for (i = 1; i < M; i++)
  {
    const E *v = u + 3 * i;
    const E t = v[0] - E_SQRT(v[1] * v[1] + v[2] * v[2]);
    if (t < r)
      r = t;
  }
  return r;
}

/* Whether u + a du lies in the interior of every cone. */
static int interior(const E *u, const E *du, const E a, const INT M)
{
  INT i;
  for (i = 0; i < M; i++)
  {
    const E v0 = u[3 * i] + a * du[3 * i];
    const E v1 = u[3 * i + 1] + a * du[3 * i + 1];
    const E v2 = u[3 * i + 2] + a * du[3 * i + 2];
    if (!(v0 > E_SQRT(v1 * v1 + v2 * v2)))
      return 0;
  }
  return 1;
}

/* Largest alpha with u + alpha du in every cone, capped at 2: callers only
 * use min(1, 0.99 alpha). */
static E step_to_boundary(const E *u, const E *du, const INT M)
{
  E alpha = EK(2.0);
  INT i;
  for (i = 0; i < M; i++)
  {
    const E *s = u + 3 * i, *d = du + 3 * i;
    const E a = d[0] * d[0] - d[1] * d[1] - d[2] * d[2];
    const E b = EK(2.0) * (s[0] * d[0] - s[1] * d[1] - s[2] * d[2]);
    const E c = s[0] * s[0] - s[1] * s[1] - s[2] * s[2];
    const E D = b * b - EK(4.0) * a * c;
    if (D >= EK(0.0))
    {
      const E sq = E_SQRT(D);
      const E qq = EK(-0.5) * (b + IF(b >= EK(0.0), sq, -sq));
      if (a != EK(0.0) && qq / a > EK(0.0) && qq / a < alpha)
        alpha = qq / a;
      if (qq != EK(0.0) && c / qq > EK(0.0) && c / qq < alpha)
        alpha = c / qq;
    }
    if (a == EK(0.0) && b < EK(0.0) && -c / b < alpha)
      alpha = -c / b;
    if (d[0] < EK(0.0) && -s[0] / d[0] < alpha)
      alpha = -s[0] / d[0];
  }
  return alpha;
}

static void mv3(E *out, const E *A, const E *u)
{
  out[0] = A[0] * u[0] + A[1] * u[1] + A[2] * u[2];
  out[1] = A[3] * u[0] + A[4] * u[1] + A[5] * u[2];
  out[2] = A[6] * u[0] + A[7] * u[1] + A[8] * u[2];
}

/* Nesterov-Todd scaling of one cone: W y = W^{-1} s = lambda, W and Winv
 * symmetric. */
static void nt_scaling(E *W, E *Winv, const E *s, const E *y)
{
  const E sJs = s[0] * s[0] - s[1] * s[1] - s[2] * s[2];
  const E yJy = y[0] * y[0] - y[1] * y[1] - y[2] * y[2];
  const E rs = E_SQRT(sJs), ry = E_SQRT(yJy);
  const E beta = E_SQRT(E_SQRT(sJs / yJy));
  const E gamma = E_SQRT((EK(1.0) + (s[0] * y[0] + s[1] * y[1]
      + s[2] * y[2]) / (rs * ry)) / EK(2.0));
  const E J[3] = {EK(1.0), EK(-1.0), EK(-1.0)};
  E v[3], nv;
  INT r, c;

  v[0] = (s[0] / rs + y[0] / ry) / (EK(2.0) * gamma);
  v[1] = (s[1] / rs - y[1] / ry) / (EK(2.0) * gamma);
  v[2] = (s[2] / rs - y[2] / ry) / (EK(2.0) * gamma);
  nv = E_SQRT(EK(2.0) * (v[0] + EK(1.0)));
  v[0] += EK(1.0);
  for (r = 0; r < 3; r++)
    v[r] /= nv;
  for (r = 0; r < 3; r++)
    for (c = 0; c < 3; c++)
    {
      const E Jrc = IF(r == c, J[r], EK(0.0));
      W[3 * r + c] = beta * (EK(2.0) * v[r] * v[c] - Jrc);
      Winv[3 * r + c] = (EK(2.0) * J[r] * v[r] * J[c] * v[c] - Jrc) / beta;
    }
}

/* x with lam o x = w in the Jordan algebra of Q^3; lam strictly interior,
 * checked by the caller. */
static void jsolve(E *x, const E *lam, const E *w)
{
  const E det = lam[0] - (lam[1] * lam[1] + lam[2] * lam[2]) / lam[0];
  x[0] = (w[0] - (lam[1] * w[1] + lam[2] * w[2]) / lam[0]) / det;
  x[1] = (w[1] - lam[1] * x[0]) / lam[0];
  x[2] = (w[2] - lam[2] * x[0]) / lam[0];
}

/* Cholesky factor of the lower triangle of A + reg I into L; 0 if A + reg I
 * is not numerically positive definite. */
static int cholesky(E *L, const E *A, const E reg, const INT n)
{
  INT j, i;
  for (j = 0; j < n; j++)
  {
    E d = A[j * n + j] + reg - dot(L + j * n, L + j * n, j);
    if (!(d > EK(0.0)))
      return 0;
    L[j * n + j] = E_SQRT(d);
    for (i = j + 1; i < n; i++)
      L[i * n + j] = (A[i * n + j] - dot(L + i * n, L + j * n, j))
          / L[j * n + j];
  }
  return 1;
}

/* Newton direction for the scaled complementarity right-hand side t (3M,
 * overwritten); L is the Cholesky factor of the normal matrix. */
static void newton(E *dx, E *dy, E *ds, E *t, const INT M, const INT n,
    const E *G, const E *Gt, const E *Winv, const E *L, const E *rx,
    const E *rz, const E *Wrz)
{
  const INT m3 = 3 * M;
  INT j, k;

  for (j = 0; j < m3; j++)
    t[j] += Wrz[j];
  for (k = 0; k < n; k++)
    dx[k] = -rx[k];
  for (j = 0; j < m3; j++)
    for (k = 0; k < n; k++)
      dx[k] -= Gt[j * n + k] * t[j];
  for (k = 0; k < n; k++)
    dx[k] = (dx[k] - dot(L + k * n, dx, k)) / L[k * n + k];
  for (k = n - 1; k >= 0; k--)
  {
    INT l;
    for (l = k + 1; l < n; l++)
      dx[k] -= L[l * n + k] * dx[l];
    dx[k] /= L[k * n + k];
  }
  for (j = 0; j < m3; j++)
  {
    t[j] += dot(Gt + j * n, dx, n);
    ds[j] = -rz[j] - dot(G + j * n, dx, n);
  }
  for (j = 0; j < M; j++)
    mv3(dy + 3 * j, Winv + 9 * j, t + 3 * j);
}

INT Y(socp3_work)(const INT M, const INT n)
{
  return 2 * n * n + 3 * M * n + 39 * M + 2 * n;
}

int Y(socp3_solve)(const INT M, const INT n, const E *h, const E *G,
    const E *q, E *x, E *y, const int warm, const int early, const E tol,
    const E eps, const INT maxit, E *work, INT *it, E *pobj, E *dobj,
    E *pres, E *dres)
{
  const INT m3 = 3 * M;
  E *Mn = work, *L = Mn + n * n, *Gt = L + n * n, *W = Gt + m3 * n;
  E *Winv = W + 9 * M, *s = Winv + 9 * M, *rz = s + m3, *Wrz = rz + m3;
  E *lam = Wrz + m3, *ds = lam + m3, *dy = ds + m3, *t = dy + m3;
  E *rx = t + m3, *dx = rx + n;
  const E hn = E_SQRT(dot(h, h, m3)), qn = E_SQRT(dot(q, q, n));
  E marg, dres_min = EK(0.0);
  int status;
  INT i, j, k, l;

  if (!warm)
    for (k = 0; k < n; k++)
      x[k] = EK(0.0);
  for (j = 0; j < m3; j++)
    s[j] = h[j] - dot(G + j * n, x, n);
  marg = cone_margin(s, M);
  if (marg <= EK(0.0))
    for (i = 0; i < M; i++)
      s[3 * i] += EK(1.0) - marg;
  if (warm)
  {
    E ymax = EK(0.0), shift;
    for (i = 0; i < M; i++)
      if (e_abs(y[3 * i]) > ymax)
        ymax = e_abs(y[3 * i]);
    shift = EK(0.1) * ymax - cone_margin(y, M);
    shift = IF(shift > EK(0.0), shift, EK(0.0)) + EK(1E-3) * ymax;
    for (i = 0; i < M; i++)
      y[3 * i] += shift;
  }
  else
    for (i = 0; i < M; i++)
    {
      y[3 * i] = EK(1.0);
      y[3 * i + 1] = y[3 * i + 2] = EK(0.0);
    }

  for (k = 0;; k++)
  {
    E gap, mu, mu_aff, sigma, a, ap, ad, dmax, reg, slack;

    for (l = 0; l < n; l++)
      rx[l] = q[l];
    for (j = 0; j < m3; j++)
    {
      for (l = 0; l < n; l++)
        rx[l] += G[j * n + l] * y[j];
      rz[j] = dot(G + j * n, x, n) + s[j] - h[j];
    }
    gap = dot(s, y, m3);
    mu = gap / (E)M;
    *pobj = dot(q, x, n);
    *dobj = -dot(h, y, m3);
    *pres = E_SQRT(dot(rz, rz, m3)) / (EK(1.0) + hn);
    *dres = E_SQRT(dot(rx, rx, n)) / (EK(1.0) + qn);
    slack = tol * (EK(1.0) + e_abs(*pobj));

    if (*pres < tol && *dres < tol && gap < slack)
    {
      status = SOCP3_OPTIMAL;
      break;
    }
    if (k == 0)
      dres_min = *dres;
    if (*dres > EK(100.0) * dres_min && dres_min < EK(10.0) * tol
        && gap < slack)
    {
      status = SOCP3_STALLED;
      break;
    }
    if (*dres < dres_min)
      dres_min = *dres;
    if (early && *dres < tol && *dobj > EK(10.0) * tol)
    {
      status = SOCP3_INFEASIBLE;
      break;
    }
    if (early && *pres < tol && *pobj < EK(-10.0) * tol)
    {
      status = SOCP3_FEASIBLE;
      break;
    }
    if (k == maxit)
    {
      status = SOCP3_MAXIT;
      break;
    }

    for (j = 0; j < n * n; j++)
      Mn[j] = EK(0.0);
    for (i = 0; i < M; i++)
    {
      E *Wi = W + 9 * i, *Wvi = Winv + 9 * i;
      nt_scaling(Wi, Wvi, s + 3 * i, y + 3 * i);
      mv3(lam + 3 * i, Wi, y + 3 * i);
      mv3(Wrz + 3 * i, Wvi, rz + 3 * i);
      for (l = 0; l < n; l++)
      {
        const E g[3] = {G[(3 * i) * n + l], G[(3 * i + 1) * n + l],
            G[(3 * i + 2) * n + l]};
        E u[3];
        mv3(u, Wvi, g);
        Gt[(3 * i) * n + l] = u[0];
        Gt[(3 * i + 1) * n + l] = u[1];
        Gt[(3 * i + 2) * n + l] = u[2];
      }
    }
    for (j = 0; j < m3; j++)
    {
      const E *g = Gt + j * n;
      for (l = 0; l < n; l++)
        for (i = 0; i <= l; i++)
          Mn[l * n + i] += g[l] * g[i];
    }
    if (!all_finite(Mn, n * n))
    {
      status = SOCP3_NUMERICAL;
      break;
    }
    dmax = EK(0.0);
    for (l = 0; l < n; l++)
      if (Mn[l * n + l] > dmax)
        dmax = Mn[l * n + l];
    reg = EK(100.0) * eps * dmax;
    for (l = 0; l < 8 && !cholesky(L, Mn, reg, n); l++)
      reg *= EK(100.0);
    if (l == 8)
    {
      status = SOCP3_NUMERICAL;
      break;
    }

    for (j = 0; j < m3; j++)
      t[j] = -lam[j];
    newton(dx, dy, ds, t, M, n, G, Gt, Winv, L, rx, rz, Wrz);
    if (!(all_finite(dx, n) && all_finite(dy, m3) && all_finite(ds, m3)))
    {
      status = SOCP3_NUMERICAL;
      break;
    }
    ap = step_to_boundary(s, ds, M);
    ad = step_to_boundary(y, dy, M);
    a = IF(ap < ad, ap, ad);
    a = IF(a < EK(1.0), a, EK(1.0));
    mu_aff = EK(0.0);
    for (j = 0; j < m3; j++)
      mu_aff += (s[j] + a * ds[j]) * (y[j] + a * dy[j]);
    mu_aff /= (E)M;
    sigma = mu_aff / mu;
    sigma = IF(sigma < EK(0.0), EK(0.0), IF(sigma > EK(1.0), EK(1.0), sigma));
    sigma = sigma * sigma * sigma;

    for (i = 0; i < M; i++)
    {
      const E *li = lam + 3 * i;
      if (!(li[0] > EK(0.0)) || li[0] * li[0] - li[1] * li[1] - li[2] * li[2]
          <= eps * li[0] * li[0])
        break;
    }
    if (i < M)
    {
      status = SOCP3_NUMERICAL;
      break;
    }
    for (i = 0; i < M; i++)
    {
      E ws[3], wy[3], w[3];
      mv3(ws, Winv + 9 * i, ds + 3 * i);
      mv3(wy, W + 9 * i, dy + 3 * i);
      w[0] = sigma * mu - dot(ws, wy, 3);
      w[1] = -(ws[0] * wy[1] + wy[0] * ws[1]);
      w[2] = -(ws[0] * wy[2] + wy[0] * ws[2]);
      jsolve(t + 3 * i, lam + 3 * i, w);
    }
    if (!all_finite(t, m3))
      for (i = 0; i < M; i++)
      {
        const E w[3] = {sigma * mu, EK(0.0), EK(0.0)};
        jsolve(t + 3 * i, lam + 3 * i, w);
      }
    if (!all_finite(t, m3))
    {
      status = SOCP3_NUMERICAL;
      break;
    }
    for (j = 0; j < m3; j++)
      t[j] -= lam[j];
    newton(dx, dy, ds, t, M, n, G, Gt, Winv, L, rx, rz, Wrz);
    if (!(all_finite(dx, n) && all_finite(dy, m3) && all_finite(ds, m3)))
    {
      status = SOCP3_NUMERICAL;
      break;
    }
    ap = step_to_boundary(s, ds, M);
    ad = step_to_boundary(y, dy, M);
    a = EK(0.99) * IF(ap < ad, ap, ad);
    a = IF(a < EK(1.0), a, EK(1.0));
    for (l = 0; l < 40 && !(interior(s, ds, a, M) && interior(y, dy, a, M));
        l++)
      a *= EK(0.5);
    for (l = 0; l < n; l++)
      x[l] += a * dx[l];
    for (j = 0; j < m3; j++)
    {
      s[j] += a * ds[j];
      y[j] += a * dy[j];
    }
  }

  *it = k;
  return status;
}
