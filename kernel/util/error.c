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

#include "infft.h"

static R cnrm1(const C *x, const INT n)
{
  INT k;
  R nrm = K(0.0);

  for (k = 0; k < n; k++)
    nrm += CABS(x[k]);

  return nrm;
}

static R nrm1(const R *x, const INT n)
{
  INT k;
  R nrm = K(0.0);

  for (k = 0; k < n; k++)
    nrm += FABS(x[k]);

  return nrm;
}

static R cerr2(const C *x, const C *y, const INT n)
{
  INT k;
  R err = K(0.0);

  if (!y)
  {
    for (k = 0; k < n; k++)
      err += CONJ(x[k]) * x[k];
  }
  else
  {
    for (k = 0; k < n; k++)
      err += CONJ(x[k]-y[k]) * (x[k]-y[k]);
  }

  return SQRT(err);
}

static R err2(const R *x, const R *y, const INT n)
{
  INT k;
  R err = K(0.0);

  if (!y)
  {
    for (k = 0; k < n; k++)
      err += x[k]*x[k];
  }
  else
  {
    for (k = 0; k < n; k++)
      err += (x[k]-y[k]) * (x[k]-y[k]);
  }

  return SQRT(err);
}

static R cerri(const C *x, const C *y, const INT n)
{
  INT k;
  R err = K(0.0), t;

  if (!y)
  {
    for (k = 0; k < n; k++)
    {
      t = CABS(x[k]);
      err = MAX(err, t);
    }
  }
  else
  {
    for (k = 0; k < n; k++)
    {
      t = CABS(x[k] - y[k]);
      err = MAX(err, t);
    }
  }

  return err;
}

static R erri(const R *x, const R *y, const INT n)
{
  INT k;
  R err = K(0.0), t;

  if (!y)
  {
    for (k = 0; k < n; k++)
    {
      t = FABS(x[k]);
      err = MAX(err, t);
    }
  }
  else
  {
    for (k = 0; k < n; k++)
    {
      t = FABS(x[k] - y[k]);
      err = MAX(err, t);
    }
  }

  return err;
}

/** computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
R Y(error_l_infty_complex)(const C *x, const C *y, const INT n)
{
  return (cerri(x, y, n)/cerri(x, 0, n));
}

/** computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
R Y(error_l_infty_double)(const R *x, const R *y, const INT n)
{
  return (erri(x, y, n)/erri(x, 0, n));
}

/** computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
R Y(error_l_infty_1_complex)(const C *x, const C *y, const INT n,
  const C *z, const INT m)
{
  return (cerri(x, y, n)/cnrm1(z, m));
}

/** computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
R Y(error_l_infty_1_double)(const R *x, const R *y, const INT n, const R *z,
  const INT m)
{
  return (erri(x, y, n)/nrm1(z, m));
}

/** computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
R Y(error_l_2_complex)(const C *x, const C *y, const INT n)
{
  return (cerr2(x, y, n)/cerr2(x, 0, n));
}

/** computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
R Y(error_l_2_double)(const R *x, const R *y, const INT n)
{
  return (err2(x, y, n)/err2(x, NULL, n));
}
