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

INT Y(exp2i)(const INT a)
{
  return (1U << a);
}

/**
 * Return floor(log2(x)) for an integer input x. If x is zero or negative this
 * method returns -1.
 *
 * see https://graphics.stanford.edu/~seander/bithacks.html#IntegerLogDeBruijn
 */
INT Y(log2i)(const INT m)
{
  /* Special case, zero or negative input. */
  if (m <= 0)
    return -1;

#if SIZEOF_PTRDIFF_T == 8
  /* Hash map with return values based on De Bruijn sequence. */
  static INT debruijn[64] =
  {
    0, 58, 1, 59, 47, 53, 2, 60, 39, 48, 27, 54, 33, 42, 3, 61, 51, 37, 40, 49,
    18, 28, 20, 55, 30, 34, 11, 43, 14, 22, 4, 62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21, 56, 45, 25, 31, 35, 16, 9, 12, 44, 24, 15,
    8, 23, 7, 6, 5, 63
  };

  register uint64_t v = (uint64_t)(m);

  /* Round down to one less than a power of 2. */
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;

  /* 0x03f6eaf2cd271461 is a hexadecimal representation of a De Bruijn
   * sequence for binary words of length 6. The binary representation
   * starts with 000000111111. This is required to make it work with one less
   * than a power of 2 instead of an actual power of 2.
   */
  return debruijn[(uint64_t)(v * 0x03f6eaf2cd271461LU) >> 58];
#elif SIZEOF_PTRDIFF_T == 4
  /* Hash map with return values based on De Bruijn sequence. */
  static INT debruijn[32] =
  {
    0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30, 8, 12, 20, 28,
    15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
  };

  register uint32_t v = (uint32_t)(m);

  /* Round down to one less than a power of 2. */
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;

  /* 0x07C4ACDD is a hexadecimal representation of a De Bruijn sequence for
   * binary words of length 5. The binary representation starts with
   * 0000011111. This is required to make it work with one less than a power of
   * 2 instead of an actual power of 2.
   */
  return debruijn[(uint32_t)(v * 0x07C4ACDDU) >> 27];
#else
#error Incompatible size of ptrdiff_t
#endif
}

/**
 * Calculate next power of two larger or equal to a given integer value.
 *
 * If the input is negative, this method returns -1. As a special case, 2 is
 * returned when the input is 1. In all other cases, the smallest power of 2
 * larger or equal to the input is returned.
 *
 * see https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
 */
INT Y(next_power_of_2)(const INT x)
{
    if (x < 0)
        return -1;
    else if (x < 2)
        return x + 1;
    else
    {
        uint64_t v = (uint64_t)x;

        /* Subtract one, so we return the input if it is a power of two. */
        v--;

        /* Round down to one less than a power of two. */
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
#if SIZEOF_PTRDIFF_T >= 2
        v |= v >> 8;
#endif
#if SIZEOF_PTRDIFF_T >= 4
        v |= v >> 16;
#endif
#if SIZEOF_PTRDIFF_T >= 8
        v |= v >> 32;
#endif
        /* Add one to get power of two. */
        v++;

        return v;
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

void Y(next_power_of_2_exp_int)(const int N, int *N2, int *t)
{
  int n,i,logn;
  int N_is_not_power_of_2=0;

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
