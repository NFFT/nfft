/*
 * Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
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

INT Y(exp2i)(const INT a)
{
  return (1U << a);
}

INT Y(log2i)(const INT m)
{
  INT l = 0;
  INT mm = m;

  while (mm > 0)
  {
    mm = (mm >> 1);
    l++;
  }
  return (l-1);
}

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
INT Y(next_power_of_2)(const INT N)
{
  INT n,i,logn;
  INT N_is_not_power_of_2=0;

  if (N == 0)
    return 1;
  else if (N == 1)
    return 2;
  else
  {
    n = N;
    logn = 0;
    while (n != 1)
    {
      if (n%2 == 1)
        N_is_not_power_of_2=1;
      n = n/2;
      logn++;
    }

    if (!N_is_not_power_of_2)
      logn--;

    for (i = 0; i <= logn; i++)
      n = n*2;

    return n;
  }
}

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
void Y(next_power_of_2_exp)(const INT N, INT *N2, INT *t)
{
  INT n,i,logn;
  INT N_is_not_power_of_2=0;

  if (N == 0)
  {
    *N2 = 1;
    *t = 0;
  }
  else
  {
    n = N;
    logn = 0;
    while (n != 1)
    {
      if (n%2 == 1)
      {
          N_is_not_power_of_2=1;
      }
      n = n/2;
      logn++;
    }

    if (!N_is_not_power_of_2)
    {
      logn--;
    }

    for (i = 0; i <= logn; i++)
    {
      n = n*2;
    }

    *N2 = n;
    *t = logn+1;
  }
}
