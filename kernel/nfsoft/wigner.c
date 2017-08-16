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
#include "wigner.h"
#include "infft.h"

double SO3_alpha(const int m1, const int m2, const int j)
{
  const int M = MAX(ABS(m1),ABS(m2)), mini = MIN(ABS(m1),ABS(m2));

  if (j < 0)
    return K(0.0);
  else if (j == 0)
  {
    if (m1 == 0 && m2 == 0)
      return K(1.0);
    if (m1 == m2)
      return K(0.5);
    else
      return IF((m1+m2)%2,K(0.0),K(-0.5));
  }
  else if (j < M - mini)
    return IF(j%2,K(0.5),K(-0.5));
  else if (j < M)
    return K(0.5) * SIGNF((R)m1)*SIGNF((R)m2);
  else
    return
      SQRT(((R)(j+1))/((R)(j+1-m1)))
      * SQRT(((R)(2*j+1))/((R)(j+1+m1)))
      * SQRT(((R)(j+1))/((R)(j+1-m2)))
      * SQRT(((R)(2*j+1))/((R)(j+1+m2)));
}

double SO3_beta(const int m1, const int m2, const int j)
{
  if (j < 0)
    return K(0.0);
  else if (j < MAX(ABS(m1),ABS(m2)))
    return K(0.5);
  else if (m1 == 0 || m2 == 0)
    return K(0.0);
  else
  {
    const R m1a = FABS((R)m1), m2a = FABS((R)m2);
    return -COPYSIGN(
      ((SQRT(m1a)*SQRT(m2a))/((R)j))
      * SQRT(m1a/((R)(j+1-m1)))
      * SQRT(((R)(2*j+1))/((R)(j+1+m1)))
      * SQRT(m2a/((R)(j+1-m2)))
      * SQRT(((R)(2*j+1))/((R)(j+1+m2))),
      SIGNF((R)m1)*SIGNF((R)m2));
  }
}

double SO3_gamma(const int m1, const int m2, const int j)
{
  if (MAX(ABS(m1),ABS(m2)) < j)
    return -(((R)(j+1))/((R)j)) * SQRT((((R)(j-m1))/((R)(j+1-m1)))
        *(((R)(j+m1))/((R)(j+1+m1)))*(((R)(j-m2))/((R)(j+1-m2)))
        *(((R)(j+m2))/((R)(j+1+m2))));
  else if (j == -1)
    return IF(m1 > m2 || !((m1 + m2) % 2), K(1.0), K(-1.0))
      * nfft_lambda2((R)ABS(m2 - m1),(R)ABS(m2 + m1));
  else
    return K(0.0);
}

/*compute the coefficients for all degrees*/

inline void SO3_alpha_row(double *alpha, int N, int k, int m)
{
  int j;
  double *alpha_act = alpha;
  for (j = -1; j <= N; j++)
    *alpha_act++ = SO3_alpha(k, m, j);
}

inline void SO3_beta_row(double *beta, int N, int k, int m)
{
  int j;
  double *beta_act = beta;
  for (j = -1; j <= N; j++)
    *beta_act++ = SO3_beta(k, m, j);
}

inline void SO3_gamma_row(double *gamma, int N, int k, int m)
{
  int j;
  double *gamma_act = gamma;
  for (j = -1; j <= N; j++)
    *gamma_act++ = SO3_gamma(k, m, j);
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
