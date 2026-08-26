/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* Table-sizing helpers for the wisdom store. Sizes are small, so trial
 * division is plenty. */

int Y(is_prime)(INT n) {
  INT d;
  if (n < 2)
    return 0;
  if (n < 4)
    return 1; /* 2 and 3 */
  if (n % 2 == 0)
    return 0;
  for (d = 3; d * d <= n; d += 2)
    if (n % d == 0)
      return 0;
  return 1;
}

/* Smallest prime >= n. */
INT Y(next_prime)(INT n) {
  INT c;
  if (n <= 2)
    return 2;
  c = (n % 2 == 0) ? n + 1 : n; /* start at an odd candidate >= n */
  while (!Y(is_prime)(c))
    c += 2;
  return c;
}
