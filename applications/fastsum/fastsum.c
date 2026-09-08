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

/*! \file fastsum.c
 *  \brief Fast NFFT-based summation algorithm.
 *
 *  \author Markus Fenn
 *  \date 2003-2006
 */
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "nfft3.h"
#include "fastsum.h"

// Required for test if (ths->k == one_over_x)
#include "kernels.h"

/**
 * \addtogroup applications_fastsum
 * \{
 */

/** max */
static int max_i(int a, int b)
{
  return a >= b ? a : b;
}

/** factorial */
static NFFT_R fak(int n)
{
  if (n <= 1)
    return NFFT_K(1.0);
  else
    return (NFFT_R)(n) * fak(n - 1);
}

/** binomial coefficient */
static NFFT_R binom(int n, int m)
{
  return fak(n) / fak(m) / fak(n - m);
}

/** basis polynomial for regularized kernel */
static NFFT_R BasisPoly(int m, int r, NFFT_R xx)
{
  int k;
  NFFT_R sum = NFFT_K(0.0);

  for (k = 0; k <= m - r; k++)
  {
    sum += binom(m + k, k) * NFFT_POW((xx + NFFT_K(1.0)) / NFFT_K(2.0), (NFFT_R) k);
  }
  return sum * NFFT_POW((xx + NFFT_K(1.0)), (NFFT_R) r) * NFFT_POW(NFFT_K(1.0) - xx, (NFFT_R) (m + 1))
      / (NFFT_R)(1 << (m + 1)) / fak(r); /* 1<<(m+1) = 2^(m+1) */
}

/** regularized kernel with K_I arbitrary and K_B smooth to zero */
NFFT_C regkern(kernel k, NFFT_R xx, int p, const NFFT_R *param, NFFT_R a, NFFT_R b)
{
  int r;
  NFFT_C sum = NFFT_K(0.0);

  if (xx < -NFFT_K(0.5))
    xx = -NFFT_K(0.5);
  if (xx > NFFT_K(0.5))
    xx = NFFT_K(0.5);
  if ((xx >= -NFFT_K(0.5) + b && xx <= -a) || (xx >= a && xx <= NFFT_K(0.5) - b))
  {
    return k(xx, 0, param);
  }
  else if (xx < -NFFT_K(0.5) + b)
  {
    sum = (k(-NFFT_K(0.5), 0, param) + k(NFFT_K(0.5), 0, param)) / NFFT_K(2.0)
        * BasisPoly(p - 1, 0, NFFT_K(2.0) * xx / b + (NFFT_K(1.0) - b) / b);
    for (r = 0; r < p; r++)
    {
      sum += NFFT_POW(-b / NFFT_K(2.0), (NFFT_R) r) * k(-NFFT_K(0.5) + b, r, param)
          * BasisPoly(p - 1, r, -NFFT_K(2.0) * xx / b + (b - NFFT_K(1.0)) / b);
    }
    return sum;
  }
  else if ((xx > -a) && (xx < a))
  {
    for (r = 0; r < p; r++)
    {
      sum +=
          NFFT_POW(a, (NFFT_R) r)
              * (k(-a, r, param) * BasisPoly(p - 1, r, xx / a)
                  + k(a, r, param) * BasisPoly(p - 1, r, -xx / a)
                      * (r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx > NFFT_K(0.5) - b)
  {
    sum = (k(-NFFT_K(0.5), 0, param) + k(NFFT_K(0.5), 0, param)) / NFFT_K(2.0)
        * BasisPoly(p - 1, 0, -NFFT_K(2.0) * xx / b + (NFFT_K(1.0) - b) / b);
    for (r = 0; r < p; r++)
    {
      sum += NFFT_POW(b / NFFT_K(2.0), (NFFT_R) r) * k(NFFT_K(0.5) - b, r, param)
          * BasisPoly(p - 1, r, NFFT_K(2.0) * xx / b - (NFFT_K(1.0) - b) / b);
    }
    return sum;
  }
  return k(xx, 0, param);
}

/** regularized kernel with K_I arbitrary and K_B periodized
 *  (used in 1D)
 */
static NFFT_C regkern1(kernel k, NFFT_R xx, int p, const NFFT_R *param, NFFT_R a, NFFT_R b)
{
  int r;
  NFFT_C sum = NFFT_K(0.0);

  if (xx < -NFFT_K(0.5))
    xx = -NFFT_K(0.5);
  if (xx > NFFT_K(0.5))
    xx = NFFT_K(0.5);
  if ((xx >= -NFFT_K(0.5) + b && xx <= -a) || (xx >= a && xx <= NFFT_K(0.5) - b))
  {
    return k(xx, 0, param);
  }
  else if ((xx > -a) && (xx < a))
  {
    for (r = 0; r < p; r++)
    {
      sum +=
          NFFT_POW(a, (NFFT_R) r)
              * (k(-a, r, param) * BasisPoly(p - 1, r, xx / a)
                  + k(a, r, param) * BasisPoly(p - 1, r, -xx / a)
                      * (r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx < -NFFT_K(0.5) + b)
  {
    for (r = 0; r < p; r++)
    {
      sum += NFFT_POW(b, (NFFT_R) r)
          * (k(NFFT_K(0.5) - b, r, param) * BasisPoly(p - 1, r, (xx + NFFT_K(0.5)) / b)
              + k(-NFFT_K(0.5) + b, r, param) * BasisPoly(p - 1, r, -(xx + NFFT_K(0.5)) / b)
                  * (r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx > NFFT_K(0.5) - b)
  {
    for (r = 0; r < p; r++)
    {
      sum += NFFT_POW(b, (NFFT_R) r)
          * (k(NFFT_K(0.5) - b, r, param) * BasisPoly(p - 1, r, (xx - NFFT_K(0.5)) / b)
              + k(-NFFT_K(0.5) + b, r, param) * BasisPoly(p - 1, r, -(xx - NFFT_K(0.5)) / b)
                  * (r & 1 ? -1 : 1));
    }
    return sum;
  }
  return k(xx, 0, param);
}

/** regularized kernel for even kernels with K_I even and K_B mirrored */

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth to K(1/2) (used in dD, d>1)
 */
static NFFT_C regkern3(kernel k, NFFT_R xx, int p, const NFFT_R *param, NFFT_R a, NFFT_R b)
{
  int r;
  NFFT_C sum = NFFT_K(0.0);

  xx = NFFT_FABS(xx);

  if (xx >= NFFT_K(0.5))
  {
    xx = NFFT_K(0.5);
  }
  /* else */
  if ((a <= xx) && (xx <= NFFT_K(0.5) - b))
  {
    return k(xx, 0, param);
  }
  else if (xx < a)
  {
    for (r = 0; r < p; r++)
    {
      sum += NFFT_POW(-a, (NFFT_R) r) * k(a, r, param)
          * (BasisPoly(p - 1, r, xx / a) + BasisPoly(p - 1, r, -xx / a));
    }
    return sum;
  }
  else if ((NFFT_K(0.5) - b < xx) && (xx <= NFFT_K(0.5)))
  {
    sum = k(NFFT_K(0.5), 0, param) * BasisPoly(p - 1, 0, -NFFT_K(2.0) * xx / b + (NFFT_K(1.0) - b) / b);
    for (r = 0; r < p; r++)
    {
      sum += NFFT_POW(b / NFFT_K(2.0), (NFFT_R) r) * k(NFFT_K(0.5) - b, r, param)
          * BasisPoly(p - 1, r, NFFT_K(2.0) * xx / b - (NFFT_K(1.0) - b) / b);
    }
    return sum;
  }
  return NFFT_K(0.0);
}

/** linear spline interpolation in near field with even kernels */

/** cubic spline interpolation in near field with even kernels */
NFFT_C kubintkern(const NFFT_R x, const NFFT_C *Add, const int Ad, const NFFT_R a)
{
  NFFT_R c, c1, c2, c3, c4;
  int r;
  NFFT_C f0, f1, f2, f3;
  c = x * (NFFT_R)(Ad) / a;
  r = (int)(NFFT_LRINT(c));
  r = abs(r);
  if (r == 0)
  {
    f0 = Add[r + 1];
    f1 = Add[r];
    f2 = Add[r + 1];
    f3 = Add[r + 2];
  }
  else
  {
    f0 = Add[r - 1];
    f1 = Add[r];
    f2 = Add[r + 1];
    f3 = Add[r + 2];
  }
  c = NFFT_FABS(c);
  c1 = c - (NFFT_R)(r);
  c2 = c1 + NFFT_K(1.0);
  c3 = c1 - NFFT_K(1.0);
  c4 = c1 - NFFT_K(2.0);
  return (-f0 * c1 * c3 * c4 / NFFT_K(6.0) + f1 * c2 * c3 * c4 / NFFT_K(2.0)
      - f2 * c2 * c1 * c4 / NFFT_K(2.0) + f3 * c2 * c1 * c3 / NFFT_K(6.0));
}

/** cubic spline interpolation in near field with arbitrary kernels */
static NFFT_C kubintkern1(const NFFT_R x, const NFFT_C *Add, const int Ad, const NFFT_R a)
{
  NFFT_R c, c1, c2, c3, c4;
  int r;
  NFFT_C f0, f1, f2, f3;
  Add += 2;
  c = (x + a) * (NFFT_R)(Ad) / NFFT_K(2.0) / a;
  r = (int)(NFFT_LRINT(c));
  r = abs(r);
  /*if (r==0) {f0=Add[r];f1=Add[r];f2=Add[r+1];f3=Add[r+2];}
   else */
  {
    f0 = Add[r - 1];
    f1 = Add[r];
    f2 = Add[r + 1];
    f3 = Add[r + 2];
  }
  c = NFFT_FABS(c);
  c1 = c - (NFFT_R)(r);
  c2 = c1 + NFFT_K(1.0);
  c3 = c1 - NFFT_K(1.0);
  c4 = c1 - NFFT_K(2.0);
  return (-f0 * c1 * c3 * c4 / NFFT_K(6.0) + f1 * c2 * c3 * c4 / NFFT_K(2.0)
      - f2 * c2 * c1 * c4 / NFFT_K(2.0) + f3 * c2 * c1 * c3 / NFFT_K(6.0));
}

/** quicksort algorithm for source knots and associated coefficients */
static void quicksort(int d, int t, NFFT_R *x, NFFT_C *alpha, int *permutation_x_alpha, int N)
{
  int lpos = 0;
  int rpos = N - 1;
  NFFT_R pivot = x[(N / 2) * d + t];

  int k;
  NFFT_R temp1;
  NFFT_C temp2;
  int temp_int;

  while (lpos <= rpos)
  {
    while (x[lpos * d + t] < pivot)
      lpos++;
    while (x[rpos * d + t] > pivot)
      rpos--;
    if (lpos <= rpos)
    {
      for (k = 0; k < d; k++)
      {
        temp1 = x[lpos * d + k];
        x[lpos * d + k] = x[rpos * d + k];
        x[rpos * d + k] = temp1;
      }
      temp2 = alpha[lpos];
      alpha[lpos] = alpha[rpos];
      alpha[rpos] = temp2;
      
      if (permutation_x_alpha)   /** store the permutation of x */
      {
        temp_int = permutation_x_alpha[lpos];
        permutation_x_alpha[lpos] = permutation_x_alpha[rpos];
        permutation_x_alpha[rpos] = temp_int;      
      }
      
      lpos++;
      rpos--;
    }
  }
  if (0 < rpos)
    quicksort(d, t, x, alpha, permutation_x_alpha, rpos + 1);
  if (lpos < N - 1)
    quicksort(d, t, x + lpos * d, alpha + lpos, permutation_x_alpha ? permutation_x_alpha + lpos : NULL, N - lpos);
}

/** initialize box-based search data structures */
static void BuildBox(fastsum_plan *ths)
{
  int t, l;
  int *box_index;
  NFFT_R val[ths->d];

  box_index = (int *) NFFT(malloc)((size_t)(ths->box_count) * sizeof(int));
  for (t = 0; t < ths->box_count; t++)
    box_index[t] = 0;

  for (l = 0; l < ths->N_total; l++)
  {
    int ind = 0;
    for (t = 0; t < ths->d; t++)
    {
      val[t] = ths->x[ths->d * l + t] + NFFT_K(0.25) - ths->eps_B / NFFT_K(2.0);
      ind *= ths->box_count_per_dim;
      ind += (int) (val[t] / ths->eps_I);
    }
    box_index[ind]++;
  }

  ths->box_offset[0] = 0;
  for (t = 1; t <= ths->box_count; t++)
  {
    ths->box_offset[t] = ths->box_offset[t - 1] + box_index[t - 1];
    box_index[t - 1] = ths->box_offset[t - 1];
  }

  for (l = 0; l < ths->N_total; l++)
  {
    int ind = 0;
    for (t = 0; t < ths->d; t++)
    {
      val[t] = ths->x[ths->d * l + t] + NFFT_K(0.25) - ths->eps_B / NFFT_K(2.0);
      ind *= ths->box_count_per_dim;
      ind += (int) (val[t] / ths->eps_I);
    }

    ths->box_alpha[box_index[ind]] = ths->alpha[l];

    for (t = 0; t < ths->d; t++)
    {
      ths->box_x[ths->d * box_index[ind] + t] = ths->x[ths->d * l + t];
    }
    box_index[ind]++;
  }
  NFFT(free)(box_index);
}

/** inner computation function for box-based near field correction */
static inline NFFT_C calc_SearchBox(int d, NFFT_R *y, NFFT_R *x, NFFT_C *alpha, int start,
    int end_lt, const NFFT_C *Add, const int Ad, int p, NFFT_R a, const kernel k,
    const NFFT_R *param, const unsigned flags)
{
  NFFT_C result = NFFT_K(0.0);

  int m, l;
  NFFT_R r;

  for (m = start; m < end_lt; m++)
  {
    if (d == 1)
    {
      r = y[0] - x[m];
    }
    else
    {
      r = NFFT_K(0.0);
      for (l = 0; l < d; l++)
        r += (y[l] - x[m * d + l]) * (y[l] - x[m * d + l]);
      r = NFFT_SQRT(r);
    }
    if (NFFT_FABS(r) < a)
    {
      result += alpha[m] * k(r, 0, param); /* alpha*(kern-regkern) */
      if (d == 1)
      {
        if (flags & EXACT_NEARFIELD)
          result -= alpha[m] * regkern1(k, r, p, param, a, NFFT_K(1.0) / NFFT_K(16.0)); /* exact value (in 1D)  */
        else
          result -= alpha[m] * kubintkern1(r, Add, Ad, a); /* spline approximation */
      }
      else
      {
        if (flags & EXACT_NEARFIELD)
          result -= alpha[m] * regkern(k, r, p, param, a, NFFT_K(1.0) / NFFT_K(16.0)); /* exact value (in dD)  */
        else
#if defined(NF_KUB)
          result -= alpha[m] * kubintkern(r, Add, Ad, a); /* spline approximation */
#elif defined(NF_QUADR)
        result -= alpha[m]*quadrintkern(r,Add,Ad,a); /* spline approximation */
#elif defined(NF_LIN)
        result -= alpha[m]*linintkern(r,Add,Ad,a); /* spline approximation */
#else
#error define interpolation method
#endif
      }
    }
  }
  return result;
}

/** box-based near field correction */
static NFFT_C SearchBox(NFFT_R *y, fastsum_plan *ths)
{
  NFFT_C val = NFFT_K(0.0);
  int t;
  int y_multiind[ths->d];
  int multiindex[ths->d];
  int y_ind;

  for (t = 0; t < ths->d; t++)
  {
    y_multiind[t] = (int)(NFFT_LRINT((y[t] + NFFT_K(0.25) - ths->eps_B / NFFT_K(2.0)) / ths->eps_I));
  }

  if (ths->d == 1)
  {
    for (y_ind = max_i(0, y_multiind[0] - 1);
        y_ind < ths->box_count_per_dim && y_ind <= y_multiind[0] + 1; y_ind++)
    {
      val += calc_SearchBox(ths->d, y, ths->box_x, ths->box_alpha,
          ths->box_offset[y_ind], ths->box_offset[y_ind + 1], ths->Add, ths->Ad,
          ths->p, ths->eps_I, ths->k, ths->kernel_param, ths->flags);
    }
  }
  else if (ths->d == 2)
  {
    for (multiindex[0] = max_i(0, y_multiind[0] - 1);
        multiindex[0] < ths->box_count_per_dim
            && multiindex[0] <= y_multiind[0] + 1; multiindex[0]++)
      for (multiindex[1] = max_i(0, y_multiind[1] - 1);
          multiindex[1] < ths->box_count_per_dim
              && multiindex[1] <= y_multiind[1] + 1; multiindex[1]++)
      {
        y_ind = (ths->box_count_per_dim * multiindex[0]) + multiindex[1];
        val += calc_SearchBox(ths->d, y, ths->box_x, ths->box_alpha,
            ths->box_offset[y_ind], ths->box_offset[y_ind + 1], ths->Add,
            ths->Ad, ths->p, ths->eps_I, ths->k, ths->kernel_param, ths->flags);
      }
  }
  else if (ths->d == 3)
  {
    for (multiindex[0] = max_i(0, y_multiind[0] - 1);
        multiindex[0] < ths->box_count_per_dim
            && multiindex[0] <= y_multiind[0] + 1; multiindex[0]++)
      for (multiindex[1] = max_i(0, y_multiind[1] - 1);
          multiindex[1] < ths->box_count_per_dim
              && multiindex[1] <= y_multiind[1] + 1; multiindex[1]++)
        for (multiindex[2] = max_i(0, y_multiind[2] - 1);
            multiindex[2] < ths->box_count_per_dim
                && multiindex[2] <= y_multiind[2] + 1; multiindex[2]++)
        {
          y_ind = ((ths->box_count_per_dim * multiindex[0]) + multiindex[1])
              * ths->box_count_per_dim + multiindex[2];
          val += calc_SearchBox(ths->d, y, ths->box_x, ths->box_alpha,
              ths->box_offset[y_ind], ths->box_offset[y_ind + 1], ths->Add,
              ths->Ad, ths->p, ths->eps_I, ths->k, ths->kernel_param,
              ths->flags);
        }
  }
  else
  {
    return 0.0/0.0; //exit(EXIT_FAILURE);
  }
  return val;
}

/** recursive sort of source knots dimension by dimension to get tree structure */
static void BuildTree(int d, int t, NFFT_R *x, NFFT_C *alpha, int *permutation_x_alpha, int N)
{
  if (N > 1)
  {
    int m = N / 2;

    quicksort(d, t, x, alpha, permutation_x_alpha, N);

    BuildTree(d, (t + 1) % d, x, alpha, permutation_x_alpha, m);
    BuildTree(d, (t + 1) % d, x + (m + 1) * d, alpha + (m + 1), permutation_x_alpha ? permutation_x_alpha + (m + 1) : NULL, N - m - 1);
  }
}

/** fast search in tree of source knots for near field computation*/
static NFFT_C SearchTree(const int d, const int t, const NFFT_R *x, const NFFT_C *alpha,
    const NFFT_R *xmin, const NFFT_R *xmax, const int N, const kernel k, const NFFT_R *param,
    const int Ad, const NFFT_C *Add, const int p, const unsigned flags)
{
  if (N == 0)
  {
      return NFFT_K(0.0);
  }
  else
  {
      int m = N / 2;
      NFFT_R Min = xmin[t];
      NFFT_R Max = xmax[t];
      NFFT_R Median = x[m * d + t];
      NFFT_R a = NFFT_FABS(Max - Min) / 2;
      int l;
      int E = 0;
      NFFT_R r;

      if (Min > Median)
          return SearchTree(d, (t + 1) % d, x + (m + 1) * d, alpha + (m + 1), xmin,
                  xmax, N - m - 1, k, param, Ad, Add, p, flags);
      else if (Max < Median)
          return SearchTree(d, (t + 1) % d, x, alpha, xmin, xmax, m, k, param, Ad,
                  Add, p, flags);
      else
      {
          NFFT_C result = NFFT_K(0.0);
          E = 0;

          for (l = 0; l < d; l++)
          {
              if (x[m * d + l] > xmin[l] && x[m * d + l] < xmax[l])
                  E++;
          }

          if (E == d)
          {
              if (d == 1)
              {
                  r = xmin[0] + a - x[m]; /* remember: xmin+a = y */
              }
              else
              {
                  r = NFFT_K(0.0);
                  for (l = 0; l < d; l++)
                      r += (xmin[l] + a - x[m * d + l]) * (xmin[l] + a - x[m * d + l]); /* remember: xmin+a = y */
                  r = NFFT_SQRT(r);
              }
              if (NFFT_FABS(r) < a)
              {
                  result += alpha[m] * k(r, 0, param); /* alpha*(kern-regkern) */
                  if (d == 1)
                  {
                      if (flags & EXACT_NEARFIELD)
                          result -= alpha[m] * regkern1(k, r, p, param, a, NFFT_K(1.0) / NFFT_K(16.0)); /* exact value (in 1D)  */
                      else
                          result -= alpha[m] * kubintkern1(r, Add, Ad, a); /* spline approximation */
                  }
                  else
                  {
                      if (flags & EXACT_NEARFIELD)
                          result -= alpha[m] * regkern(k, r, p, param, a, NFFT_K(1.0) / NFFT_K(16.0)); /* exact value (in dD)  */
                      else
#if defined(NF_KUB)
                          result -= alpha[m] * kubintkern(r, Add, Ad, a); /* spline approximation */
#elif defined(NF_QUADR)
                      result -= alpha[m]*quadrintkern(r,Add,Ad,a); /* spline approximation */
#elif defined(NF_LIN)
                      result -= alpha[m]*linintkern(r,Add,Ad,a); /* spline approximation */
#else
#error define interpolation method
#endif
                  }
              }
          }
          result += SearchTree(d, (t + 1) % d, x + (m + 1) * d, alpha + (m + 1), xmin,
                  xmax, N - m - 1, k, param, Ad, Add, p, flags);
          result += SearchTree(d, (t + 1) % d, x, alpha, xmin, xmax, m, k, param, Ad, Add,
                        p, flags);
          return result;
      }
  }
}

static void fastsum_precompute_kernel(fastsum_plan *ths)
{
  int j, k, t;
  NFFT_INT N[ths->d];
  int n_total;
#ifdef MEASURE_TIME
  double t0, t1;
#endif

  ths->MEASURE_TIME_t[0] = NFFT_K(0.0);

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif
  /** precompute spline values for near field */
  if (ths->eps_I > 0.0 && !(ths->flags & EXACT_NEARFIELD))
  {
    if (ths->d == 1)
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k)
#endif
      for (k = -ths->Ad / 2 - 2; k <= ths->Ad / 2 + 2; k++)
        ths->Add[k + ths->Ad / 2 + 2] = regkern1(ths->k,
            ths->eps_I * (NFFT_R) k / (NFFT_R)(ths->Ad) * NFFT_K(2.0), ths->p, ths->kernel_param,
            ths->eps_I, ths->eps_B);
    else
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k)
#endif
      for (k = 0; k <= ths->Ad + 2; k++)
        ths->Add[k] = regkern3(ths->k, ths->eps_I * (NFFT_R) k / (NFFT_R)(ths->Ad), ths->p,
            ths->kernel_param, ths->eps_I, ths->eps_B);
  }
#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[0] += t1 - t0;
#endif

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif
  /** precompute Fourier coefficients of regularised kernel*/
  n_total = 1;
  for (t = 0; t < ths->d; t++)
    n_total *= ths->n;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(j,k,t)
#endif
  for (j = 0; j < n_total; j++)
  {
    if (ths->d == 1)
      ths->b[j] = regkern1(ths->k, (NFFT_R) - (j / (NFFT_R)(ths->n) - NFFT_K(0.5)), ths->p,
          ths->kernel_param, ths->eps_I, ths->eps_B) / (NFFT_R)(n_total);
    else
    {
      k = j;
      ths->b[j] = NFFT_K(0.0);
      for (t = 0; t < ths->d; t++)
      {
        ths->b[j] += ((NFFT_R) (k % (ths->n)) / (NFFT_R)(ths->n) - NFFT_K(0.5))
            * ((NFFT_R) (k % (ths->n)) / (NFFT_R)(ths->n) - NFFT_K(0.5));
        k = k / (ths->n);
      }
      ths->b[j] = regkern3(ths->k, NFFT_SQRT(NFFT_CREAL(ths->b[j])), ths->p, ths->kernel_param,
          ths->eps_I, ths->eps_B) / (NFFT_R)(n_total);
    }
  }

  for (t = 0; t < ths->d; t++)
    N[t] = ths->n;

  NFFT(fftshift_complex)(ths->b, (int)(ths->d), N);
  FFTW(execute)(ths->fft_plan);
  NFFT(fftshift_complex)(ths->b, (int)(ths->d), N);
#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[0] += t1 - t0;
#endif
}

void fastsum_init_guru_kernel(fastsum_plan *ths, int d, kernel k, NFFT_R *param,
    unsigned flags, int nn, int p, NFFT_R eps_I, NFFT_R eps_B)
{
  int t;
  int N[d];
  int n_total;
#ifdef _OPENMP
  int nthreads = NFFT(get_num_threads)();
#endif

  ths->d = d;

  ths->k = k;
  ths->kernel_param = param;

  ths->flags = flags;

  ths->p = p;
  ths->eps_I = eps_I; /* =(NFFT_R)ths->p/(NFFT_R)nn; *//** inner boundary */
  ths->eps_B = eps_B; /* =NFFT_K(1.0)/NFFT_K(16.0); *//** outer boundary */

  /** init spline for near field computation */
  if (ths->eps_I > 0.0 && !(ths->flags & EXACT_NEARFIELD))
  {
    if (ths->d == 1)
    {
      ths->Ad = 4 * (ths->p) * (ths->p);
      ths->Add = (NFFT_C *) NFFT(malloc)((size_t)(ths->Ad + 5) * (sizeof(NFFT_C)));
    }
    else
    {
      if (ths->k == one_over_x)
      {
        NFFT_R delta = NFFT_K(1e-8);
        switch (p)
        {
          case 2:
            delta = NFFT_K(1e-3);
            break;
          case 3:
            delta = NFFT_K(1e-4);
            break;
          case 4:
            delta = NFFT_K(1e-5);
            break;
          case 5:
            delta = NFFT_K(1e-6);
            break;
          case 6:
            delta = NFFT_K(1e-6);
            break;
          case 7:
            delta = NFFT_K(1e-7);
            break;
          default:
            delta = NFFT_K(1e-8);
        }

#if defined(NF_KUB)
        ths->Ad = max_i(10, (int)(NFFT_LRINT(NFFT_CEIL(NFFT_K(1.4) / NFFT_POW(delta, NFFT_K(1.0) / NFFT_K(4.0))))));
        ths->Add = (NFFT_C *) NFFT(malloc)((size_t)(ths->Ad + 3) * (sizeof(NFFT_C)));
#elif defined(NF_QUADR)
        ths->Ad = (int)(NFFT_LRINT(NFFT_CEIL(NFFT_K(2.2)/NFFT_POW(delta,NFFT_K(1.0)/NFFT_K(3.0)))));
        ths->Add = (NFFT_C *)NFFT(malloc)((size_t)(ths->Ad+3)*(sizeof(NFFT_C)));
#elif defined(NF_LIN)
        ths->Ad = (int)(NFFT_LRINT(NFFT_CEIL(NFFT_K(1.7)/pow(delta,NFFT_K(1.0)/NFFT_K(2.0)))));
        ths->Add = (NFFT_C *)NFFT(malloc)((size_t)(ths->Ad+3)*(sizeof(NFFT_C)));
#else
#error define NF_LIN or NF_QUADR or NF_KUB
#endif
      }
      else
      {
        ths->Ad = 2 * (ths->p) * (ths->p);
        ths->Add = (NFFT_C *) NFFT(malloc)((size_t)(ths->Ad + 3) * (sizeof(NFFT_C)));
      }
    } /* multi-dimensional case */
  } /* !EXACT_NEARFIELD == spline approximation in near field AND eps_I > 0 */

  ths->n = nn;
  for (t = 0; t < d; t++)
  {
    N[t] = nn;
  }

  /** init d-dimensional FFTW plan */
  n_total = 1;
  for (t = 0; t < d; t++)
    n_total *= nn;

  ths->b = (NFFT_C*) NFFT(malloc)((size_t)(n_total) * sizeof(NFFT_C));
  ths->f_hat = (NFFT_C*) NFFT(malloc)((size_t)(n_total) * sizeof(NFFT_C));
#if defined(_OPENMP) && defined(HAVE_FFTW_THREADS)
  #pragma omp critical (nfft_omp_critical_fftw_plan)
  {
    FFTW(plan_with_nthreads)(nthreads);
#endif

  ths->fft_plan = FFTW(plan_dft)(d, N, ths->b, ths->b, FFTW_FORWARD,
      FFTW_ESTIMATE);

#if defined(_OPENMP) && defined(HAVE_FFTW_THREADS)
}
#endif

  fastsum_precompute_kernel(ths);
}

void fastsum_init_guru_source_nodes(fastsum_plan *ths, int N_total, int nn_oversampled, int m)
{
  int t;
  int N[ths->d], n[ths->d];
  unsigned sort_flags_adjoint = 0U;

  if (ths->d > 1)
  {
#ifdef _OPENMP
    sort_flags_adjoint = NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT;
#else
    sort_flags_adjoint = NFFT_SORT_NODES;
#endif
  }

  ths->N_total = N_total;

  ths->x = (NFFT_R *) NFFT(malloc)((size_t)(ths->d * N_total) * (sizeof(NFFT_R)));
  ths->alpha = (NFFT_C *) NFFT(malloc)((size_t)(N_total) * (sizeof(NFFT_C)));

  /** init d-dimensional NFFT plan */
  for (t = 0; t < ths->d; t++)
  {
    N[t] = ths->n;
    n[t] = nn_oversampled;
  }

  NFFT(init_guru)(&(ths->mv1), ths->d, N, N_total, n, m,
      sort_flags_adjoint |
      PRE_PHI_HUT | PRE_PSI | /*MALLOC_X | MALLOC_F_HAT | MALLOC_F |*/ FFTW_INIT
          | ((ths->d == 1) ? FFT_OUT_OF_PLACE : 0U),
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  ths->mv1.x = ths->x;
  ths->mv1.f = ths->alpha;
  ths->mv1.f_hat = ths->f_hat;
  
  ths->box_offset = NULL;
  ths->box_alpha = NULL;
  ths->box_x = NULL;
  ths->permutation_x_alpha = NULL;

  if (ths->flags & NEARFIELD_BOXES)
  {
    if (ths->eps_I > 0.0)
    {
      ths->box_count_per_dim = (int)(NFFT_LRINT(NFFT_FLOOR((NFFT_K(0.5) - ths->eps_B) / ths->eps_I))) + 1;
      ths->box_count = 1;
      for (t = 0; t < ths->d; t++)
        ths->box_count *= ths->box_count_per_dim;

      ths->box_offset = (int *) NFFT(malloc)((size_t)(ths->box_count + 1) * sizeof(int));

      ths->box_alpha = (NFFT_C *) NFFT(malloc)((size_t)(ths->N_total) * (sizeof(NFFT_C)));

      ths->box_x = (NFFT_R *) NFFT(malloc)((size_t)(ths->d * ths->N_total) * sizeof(NFFT_R));
    } /* eps_I > 0 */
  } /* NEARFIELD_BOXES */
  else
  {
    if ((ths->flags & STORE_PERMUTATION_X_ALPHA) && (ths->eps_I > 0.0))
    {
      ths->permutation_x_alpha = (int *) NFFT(malloc)((size_t)(ths->N_total) * (sizeof(int)));
      for (int i=0; i<ths->N_total; i++)
        ths->permutation_x_alpha[i] = i;
    }
  } /* search tree */
}

void fastsum_init_guru_target_nodes(fastsum_plan *ths, int M_total, int nn_oversampled, int m)
{
  int t;
  int N[ths->d], n[ths->d];
  unsigned sort_flags_trafo = 0U;

  if (ths->d > 1)
    sort_flags_trafo = NFFT_SORT_NODES;

  ths->M_total = M_total;

  ths->y = (NFFT_R *) NFFT(malloc)((size_t)(ths->d * M_total) * (sizeof(NFFT_R)));
  ths->f = (NFFT_C *) NFFT(malloc)((size_t)(M_total) * (sizeof(NFFT_C)));

  /** init d-dimensional NFFT plan */
  for (t = 0; t < ths->d; t++)
  {
    N[t] = ths->n;
    n[t] = nn_oversampled;
  }

  NFFT(init_guru)(&(ths->mv2), ths->d, N, M_total, n, m,
      sort_flags_trafo |
      PRE_PHI_HUT | PRE_PSI | /*MALLOC_X | MALLOC_F_HAT | MALLOC_F |*/ FFTW_INIT
          | ((ths->d == 1) ? FFT_OUT_OF_PLACE : 0U),
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  ths->mv2.x = ths->y;
  ths->mv2.f = ths->f;
  ths->mv2.f_hat = ths->f_hat;
}

/** initialization of fastsum plan */
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total,
    kernel k, NFFT_R *param, unsigned flags, int nn, int m, int p, NFFT_R eps_I, NFFT_R eps_B)
{
  fastsum_init_guru_kernel(ths, d, k, param, flags, nn, p, eps_I, eps_B);
  fastsum_init_guru_source_nodes(ths, N_total, 2*nn, m);
  fastsum_init_guru_target_nodes(ths, M_total, 2*nn, m);
}

/** finalization of fastsum plan */
void fastsum_finalize_source_nodes(fastsum_plan *ths)
{
  NFFT(free)(ths->x);
  NFFT(free)(ths->alpha);

  NFFT(finalize)(&(ths->mv1));

  if (ths->flags & NEARFIELD_BOXES)
  {
    if (ths->eps_I > 0.0)
    {
      NFFT(free)(ths->box_offset);
      NFFT(free)(ths->box_alpha);
      NFFT(free)(ths->box_x);
    }
  } /* NEARFIELD_BOXES */
  else
  {
    if (ths->permutation_x_alpha)
      NFFT(free)(ths->permutation_x_alpha);
  } /* search tree */
}

/** finalization of fastsum plan */
void fastsum_finalize_target_nodes(fastsum_plan *ths)
{
  NFFT(free)(ths->y);
  NFFT(free)(ths->f);

  NFFT(finalize)(&(ths->mv2));
}

/** finalization of fastsum plan */
void fastsum_finalize_kernel(fastsum_plan *ths)
{
  if (ths->eps_I > 0.0 && !(ths->flags & EXACT_NEARFIELD))
    NFFT(free)(ths->Add);

#ifdef _OPENMP
  #pragma omp critical (nfft_omp_critical_fftw_plan)
  {
#endif
  FFTW(destroy_plan)(ths->fft_plan);
#ifdef _OPENMP
}
#endif

  NFFT(free)(ths->b);
  NFFT(free)(ths->f_hat);
}

/** finalization of fastsum plan */
void fastsum_finalize(fastsum_plan *ths)
{
  fastsum_finalize_target_nodes(ths);
  fastsum_finalize_source_nodes(ths);
  fastsum_finalize_kernel(ths);
}

/** direct computation of sums */
void fastsum_exact(fastsum_plan *ths)
{
  int j, k;
  int t;
  NFFT_R r;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(j,k,t,r)
#endif
  for (j = 0; j < ths->M_total; j++)
  {
    ths->f[j] = NFFT_K(0.0);
    for (k = 0; k < ths->N_total; k++)
    {
      if (ths->d == 1)
        r = ths->y[j] - ths->x[k];
      else
      {
        r = NFFT_K(0.0);
        for (t = 0; t < ths->d; t++)
          r += (ths->y[j * ths->d + t] - ths->x[k * ths->d + t])
              * (ths->y[j * ths->d + t] - ths->x[k * ths->d + t]);
        r = NFFT_SQRT(r);
      }
      ths->f[j] += ths->alpha[k] * ths->k(r, 0, ths->kernel_param);
    }
  }
}

/** precomputation for fastsum */
void fastsum_precompute_source_nodes(fastsum_plan *ths)
{
#ifdef MEASURE_TIME
  double t0, t1;
#endif

  ths->MEASURE_TIME_t[1] = NFFT_K(0.0);
  ths->MEASURE_TIME_t[3] = NFFT_K(0.0);

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif

  if (ths->eps_I > 0.0)
  {
    if (ths->flags & NEARFIELD_BOXES)
      BuildBox(ths);
    else
    /** sort source knots for search tree */
      BuildTree(ths->d, 0, ths->x, ths->alpha, ths->permutation_x_alpha, ths->N_total);
  } /* eps_I > 0 */

#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[3] += t1 - t0;
#endif

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif
  /** init NFFT plan for transposed transform in first step*/

  /** precompute psi, the entries of the matrix B */
  if (ths->mv1.flags & PRE_LIN_PSI)
    NFFT(precompute_lin_psi)(&(ths->mv1));

  if (ths->mv1.flags & PRE_PSI)
    NFFT(precompute_psi)(&(ths->mv1));

  if (ths->mv1.flags & PRE_FULL_PSI)
    NFFT(precompute_full_psi)(&(ths->mv1));
#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[1] += t1 - t0;
#endif

}

/** precomputation for fastsum */
void fastsum_precompute_target_nodes(fastsum_plan *ths)
{
#ifdef MEASURE_TIME
  double t0, t1;
#endif

  ths->MEASURE_TIME_t[2] = NFFT_K(0.0);

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif
  /** init NFFT plan for transform in third step*/

  /** precompute psi, the entries of the matrix B */
  if (ths->mv2.flags & PRE_LIN_PSI)
    NFFT(precompute_lin_psi)(&(ths->mv2));

  if (ths->mv2.flags & PRE_PSI)
    NFFT(precompute_psi)(&(ths->mv2));

  if (ths->mv2.flags & PRE_FULL_PSI)
    NFFT(precompute_full_psi)(&(ths->mv2));
#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[2] += t1 - t0;
#endif
}

/** precomputation for fastsum */
void fastsum_precompute(fastsum_plan *ths)
{
  fastsum_precompute_source_nodes(ths);
  fastsum_precompute_target_nodes(ths);
}

/** fast NFFT-based summation */
void fastsum_trafo(fastsum_plan *ths)
{
  int j, k, t;
#ifdef MEASURE_TIME
  double t0, t1;
#endif

  ths->MEASURE_TIME_t[4] = NFFT_K(0.0);
  ths->MEASURE_TIME_t[5] = NFFT_K(0.0);
  ths->MEASURE_TIME_t[6] = NFFT_K(0.0);
  ths->MEASURE_TIME_t[7] = NFFT_K(0.0);

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif
  /** first step of algorithm */
  NFFT(adjoint)(&(ths->mv1));
#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[4] += t1 - t0;
#endif

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif
  /** second step of algorithm */
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < ths->mv2.N_total; k++)
    ths->mv2.f_hat[k] = ths->b[k] * ths->mv1.f_hat[k];
#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[5] += t1 - t0;
#endif

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif
  /** third step of algorithm */
  NFFT(trafo)(&(ths->mv2));
#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[6] += t1 - t0;
#endif

#ifdef MEASURE_TIME
  t0 = NFFT(clock_gettime_seconds)();
#endif

  /** write far field to output */
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(j)
#endif
  for (j = 0; j < ths->M_total; j++)
    ths->f[j] = ths->mv2.f[j];

  if (ths->eps_I > 0.0)
  {
    /** add near field */
  #ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j,k,t)
  #endif
    for (j = 0; j < ths->M_total; j++)
    {
      NFFT_R ymin[ths->d], ymax[ths->d]; /** limits for d-dimensional near field box */

      if (ths->flags & NEARFIELD_BOXES)
        ths->f[j] += SearchBox(ths->y + ths->d * j, ths);
      else
      {
        for (t = 0; t < ths->d; t++)
        {
          ymin[t] = ths->y[ths->d * j + t] - ths->eps_I;
          ymax[t] = ths->y[ths->d * j + t] + ths->eps_I;
        }
        ths->f[j]
         += SearchTree(ths->d, 0, ths->x, ths->alpha, ymin, ymax, ths->N_total,
             ths->k, ths->kernel_param, ths->Ad, ths->Add, ths->p, ths->flags);
      }
    }
  }

#ifdef MEASURE_TIME
  t1 = NFFT(clock_gettime_seconds)();
  ths->MEASURE_TIME_t[7] += t1 - t0;
#endif
}
/* \} */

/* fastsum.c */
