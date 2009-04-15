/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id$ */

#include <math.h>
#include <stdio.h>
#include "infft.h"
#include "wigner.h"
#include "util.h"

double SO3_gamma(int m1, int m2, int j)
{
  double dj, dm1, dm2, M, mini;

  dj = (double) j;
  dm1 = (double) m1;
  dm2 = (double) m2;
  M = (double) (ABS(m1) > ABS(m2) ? ABS(m1) : ABS(m2));
  mini = (double) (ABS(m1) < ABS(m2) ? ABS(m1) : ABS(m2));

  if (j == -1)
  {
    /* Constant is ((2n)!)^(1/2) / (2^n n!). */
    const int mu = ABS(m2 - m1), nu = ABS(m2 + m1), La = (mu + nu) / 2;
    const R eps = IF(m1 > m2 || !((m1 + m2) % 2), K(1.0), K(-1.0));
    return eps * POW(K(0.5),La) * nfft_lambda2((R)mu,(R)nu);
  }
  else if (j <= M)
  {
    return (0.0);
  }
  else
    return -(((R)(j+1))/((R)j)) *
    SQRT(
       (((R)(j-m1))/((R)(j+1-m1)))
      *(((R)(j+m1))/((R)(j+1+m1)))
      *(((R)(j-m2))/((R)(j+1-m2)))
      *(((R)(j+m2))/((R)(j+1+m2))));
}

inline double SO3_alpha(int m1, int m2, int j)
{
  double dj, dm1, dm2, M, mini, neg;
  dj = (double) j;
  dm1 = (double) m1;
  dm2 = (double) m2;
  M = (double) (ABS(m1) > ABS(m2) ? ABS(m1) : ABS(m2));
  mini = (double) (ABS(m1) < ABS(m2) ? ABS(m1) : ABS(m2));

  /**(-1) - Faktor der für k xor m negativ gebraucht wird */

  if ((dm1 < 0 && dm2 >= 0) || (dm2 < 0 && dm1 >= 0))
  {
    neg = -1.0;
  }
  else
  {
    neg = 1.0;
  }

  if (j == -1)
  {
    return (0.0);
  }
  else if (j == 0)
  {
    if (dm1 == dm2)
    {
      return 1.0;
    }
    else
    {
      return (int) (dm1 + dm2) % 2 == 0 ? -1.0 : 0.0;
    }
  }
  else if (j < M - mini)
  {
    return j % 2 == 0 ? -1.0 : 1.0;
  }
  else if (j < M)
  {
    return 1.0 * neg;
  }
  else
    return
      SQRT(
        (((R)(j+1))/((R)(j+1-m1)))
        *(((R)(j+1))/((R)(j+1+m1)))
        *(((R)(2*j+1))/((R)(j+1-m2)))
        *(((R)(2*j+1))/((R)(j+1+m2))));
}

double SO3_beta(int m1, int m2, int j)
{
  double dj, dm1, dm2, M, mini;

  dj = (double) j;
  dm1 = (double) m1;
  dm2 = (double) m2;
  M = (double) (ABS(m1) > ABS(m2) ? ABS(m1) : ABS(m2));
  mini = (double) (ABS(m1) < ABS(m2) ? ABS(m1) : ABS(m2));

  if (0 <= j && j < M)
    return (1.0);
  else if (dm1 == 0. || dm2 == 0.)
    return (0.0);
  else if (j < 0)
    return K(0.0);
  else
    return -COPYSIGN(
      SQRT(
        (((R)(m1*m2))/((R)(j+1-m1)))
        *(((R)(m1*m2))/((R)(j+1+m1)))
        *(((R)(2*j+1))/((R)(j+1-m2)))
        *(((R)(2*j+1))/((R)(j+1+m2))))/((R)j),
      (R)(m1*m2));
}

/*compute the coefficients for all degrees*/

inline void SO3_alpha_row(double *alpha, int N, int k, int m)
{
  int j;
  double *alpha_act = alpha;
  for (j = -1; j <= N; j++)
  {
    *alpha_act = SO3_alpha(k, m, j);
    alpha_act++;
  }
}

inline void SO3_beta_row(double *beta, int N, int k, int m)
{
  int j;
  double *beta_act = beta;
  for (j = -1; j <= N; j++)
  {
    *beta_act = SO3_beta(k, m, j);
    beta_act++;
  }
}

inline void SO3_gamma_row(double *gamma, int N, int k, int m)
{
  int j;
  double *gamma_act = gamma;
  for (j = -1; j <= N; j++)
  {
    *gamma_act = SO3_gamma(k, m, j);
    gamma_act++;
  }
}

/*compute for all degrees l and orders k*/

inline void SO3_alpha_matrix(double *alpha, int N, int m)
{
  int i, j;
  double *alpha_act = alpha;
  for (i = -N; i <= N; i++)
  {
    for (j = -1; j <= N; j++)
    {
      *alpha_act = SO3_alpha(i, m, j);
      alpha_act++;
    }
  }
}

inline void SO3_beta_matrix(double *alpha, int N, int m)
{
  int i, j;
  double *alpha_act = alpha;
  for (i = -N; i <= N; i++)
  {
    for (j = -1; j <= N; j++)
    {
      *alpha_act = SO3_beta(i, m, j);
      alpha_act++;
    }
  }
}

inline void SO3_gamma_matrix(double *alpha, int N, int m)
{
  int i, j;
  double *alpha_act = alpha;
  for (i = -N; i <= N; i++)
  {
    for (j = -1; j <= N; j++)
    {
      *alpha_act = SO3_gamma(i, m, j);
      alpha_act++;
    }
  }
}

/*compute all 3termrecurrence coeffs*/

inline void SO3_alpha_all(double *alpha, int N)
{
  int q;
  int i, j, m;
  double *alpha_act = alpha;
  q = 0;
  for (m = -N; m <= N; m++)
  {
    for (i = -N; i <= N; i++)
    {
      for (j = -1; j <= N; j++)
      {
        *alpha_act = SO3_alpha(i, m, j);
        fprintf(stdout, "alpha_all_%d^[%d,%d]=%f\n", j, i, m,
            SO3_alpha(i, m, j));
        alpha_act++;
        q = q + 1;

      }
    }
  }
}

inline void SO3_beta_all(double *alpha, int N)
{
  int i, j, m;
  double *alpha_act = alpha;
  for (m = -N; m <= N; m++)
  {
    for (i = -N; i <= N; i++)
    {
      for (j = -1; j <= N; j++)
      {
        *alpha_act = SO3_beta(i, m, j);
        alpha_act++;
      }
    }
  }
}

inline void SO3_gamma_all(double *alpha, int N)
{
  int i, j, m;
  double *alpha_act = alpha;
  for (m = -N; m <= N; m++)
  {
    for (i = -N; i <= N; i++)
    {
      for (j = -1; j <= N; j++)
      {
        *alpha_act = SO3_gamma(i, m, j);
        alpha_act++;
      }
    }
  }
}

inline void eval_wigner(double *x, double *y, int size, int k, double *alpha,
    double *beta, double *gamma)
{
  /* Evaluate the wigner function d_{k,nleg} (l,x) for the vector
   * of knots  x[0], ..., x[size-1] by the Clenshaw algorithm
   */
  int i, j;
  double a, b, x_val_act, a_old;
  double *x_act, *y_act;
  double *alpha_act, *beta_act, *gamma_act;

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
        a = b + a_old * ((*alpha_act) * x_val_act + (*beta_act));
        b = a_old * (*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a * ((*alpha_act) * x_val_act + (*beta_act)) + b);
    }
    x_act++;
    y_act++;
  }
}

inline int eval_wigner_thresh(double *x, double *y, int size, int k,
    double *alpha, double *beta, double *gamma, double threshold)
{

  int i, j;
  double a, b, x_val_act, a_old;
  double *x_act, *y_act;
  double *alpha_act, *beta_act, *gamma_act;

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
        a = b + a_old * ((*alpha_act) * x_val_act + (*beta_act));
        b = a_old * (*gamma_act);
        alpha_act--;
        beta_act--;
        gamma_act--;
      }
      *y_act = (a * ((*alpha_act) * x_val_act + (*beta_act)) + b);
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

/************************************************************************/
/* L2 normed wigner little d, WHERE THE DEGREE OF THE FUNCTION IS EQUAL
 TO ONE OF ITS ORDERS. This is the function to use when starting the
 three-term recurrence at orders (m1,m2)

 Note that, by definition, since I am starting the recurrence with this
 function, that the degree j of the function is equal to max(abs(m1), abs(m2) ).
 */

double wigner_start(int m1, int m2, double theta)
{

  int i, l, delta;
  int cosPower, sinPower;
  int absM1, absM2;
  double dl, dm1, dm2, normFactor, sinSign;
  double dCP, dSP;
  double max;
  double min;

  max = (double) (ABS(m1) > ABS(m2) ? ABS(m1) : ABS(m2));
  min = (double) (ABS(m1) < ABS(m2) ? ABS(m1) : ABS(m2));

  l = max;
  delta = l - min;

  absM1 = ABS(m1);
  absM2 = ABS(m2);
  dl = (double) l;
  dm1 = (double) m1;
  dm2 = (double) m2;
  sinSign = 1.;
  normFactor = 1.;

  for (i = 0; i < delta; i++)
    normFactor *= SQRT((2. * dl - ((double) i)) / (((double) i) + 1.));

  /* need to adjust to make the L2-norm equal to 1 */

  normFactor *= SQRT((2. * dl + 1.) / 2.);

  if (l == absM1)
    if (m1 >= 0)
    {
      cosPower = l + m2;
      sinPower = l - m2;
      if ((l - m2) % 2)
        sinSign = -1.;
    }
    else
    {
      cosPower = l - m2;
      sinPower = l + m2;
    }
  else if (m2 >= 0)
  {
    cosPower = l + m1;
    sinPower = l - m1;
  }
  else
  {
    cosPower = l - m1;
    sinPower = l + m1;
    if ((l + m1) % 2)
      sinSign = -1.;
  }

  dCP = (double) cosPower;
  dSP = (double) sinPower;

  return normFactor * sinSign * POW(sin(theta / 2), dSP) * POW(cos(theta / 2),
      dCP);
}
