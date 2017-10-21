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

#include <math.h>
#include <stdio.h>
#include "infft.h"
#include "legendre.h"
#include "infft.h"

/* One over sqrt(pi) */
DK(KSQRTPII,0.56418958354775628694807945156077258584405062932900);

static inline R alpha_al(const int k, const int n)
{
  if (k > 0)
  {
    if (k < n)
      return IF(k%2,K(1.0),K(-1.0));
    else
      return SQRT(((R)(2*k+1))/((R)(k-n+1)))*SQRT((((R)(2*k+1))/((R)(k+n+1))));
  }
  else if (k == 0)
  {
    if (n == 0)
      return K(1.0);
    else
      return IF(n%2,K(0.0),K(-1.0));
  }
  return K(0.0);
}

static inline R beta_al(const int k, const int n)
{
  if (0 <= k && k < n)
    return K(1.0);
  else
    return K(0.0);
}

static inline R gamma_al(const int k, const int n)
{
  if (k == -1)
    return SQRT(KSQRTPII*nfft_lambda((R)(n),K(0.5)));
  else if (k <= n)
    return K(0.0);
  else
    return -SQRT(((R)(k-n))/((R)(k-n+1))*((R)(k+n))/((R)(k+n+1)));
}

void alpha_al_row(R *alpha, const int N, const int n)
{
  int j;
  R *p = alpha;
  for (j = -1; j <= N; j++)
    *p++ = alpha_al(j,n);
}

void beta_al_row(R *beta, const int N, const int n)
{
  int j;
  R *p = beta;
  for (j = -1; j <= N; j++)
    *p++ = beta_al(j,n);
}

void gamma_al_row(R *gamma, const int N, const int n)
{
  int j;
  R *p = gamma;
  for (j = -1; j <= N; j++)
    *p++ = gamma_al(j,n);
}

inline void alpha_al_all(R *alpha, const int N)
{
  int i,j;
  R *p = alpha;
  for (i = 0; i <= N; i++)
    for (j = -1; j <= N; j++)
      *p++ = alpha_al(j,i);
}

inline void beta_al_all(R *alpha, const int N)
{
  int i,j;
  R *p = alpha;
  for (i = 0; i <= N; i++)
    for (j = -1; j <= N; j++)
      *p++ = beta_al(j,i);
}

inline void gamma_al_all(R *alpha, const int N)
{
  int i,j;
  R *p = alpha;
  for (i = 0; i <= N; i++)
    for (j = -1; j <= N; j++)
      *p++ = gamma_al(j,i);
}

void eval_al(R *x, R *y, const int size, const int k, R *alpha,
  R *beta, R *gamma)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  R a,b,x_val_act,a_old;
  R *x_act, *y_act;
  R *alpha_act, *beta_act, *gamma_act;

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;

    if (k == 0)
    {
      *y_act = 1.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gamma[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
	       b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
    }
    x_act++;
    y_act++;
  }
}

int eval_al_thresh(R *x, R *y, const int size, const int k, R *alpha,
  R *beta, R *gamma, R threshold)
{
  /* Evaluate the associated Legendre polynomial P_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i,j;
  R a,b,x_val_act,a_old;
  R *x_act, *y_act;
  R *alpha_act, *beta_act, *gamma_act;

  /* Traverse all nodes. */
  x_act = x;
  y_act = y;
  for (i = 0; i < size; i++)
  {
    a = 1.0;
    b = 0.0;
    x_val_act = *x_act;

    if (k == 0)
    {
     *y_act = 1.0;
    }
    else
    {
      alpha_act = &(alpha[k]);
      beta_act = &(beta[k]);
      gamma_act = &(gamma[k]);
      for (j = k; j > 1; j--)
      {
        a_old = a;
        a = b + a_old*((*alpha_act)*x_val_act+(*beta_act));
	       b = a_old*(*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a*((*alpha_act)*x_val_act+(*beta_act))+b);
      if (fabs(*y_act) > threshold)
      {
        return 1;
      }
    }
    x_act++;
    y_act++;
  }
  return 0;
}
