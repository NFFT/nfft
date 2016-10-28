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

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"
#include "version.h"

static INT _log2i(const INT m)
{
  INT l = 0;
  INT mm = m;

  if (m <= 0)
    return -1;

  while (mm > (INT)(0))
  {
    mm = (mm >> 1);
    l++;
  }

  return (l-1);
}

void X(check_log2i)(void)
{
    INT i;
    INT j;

    {
        INT r = Y(log2i)(0);
        int ok = r == -1;
        printf("log2i(%td) = %td -> %s\n", (INT)(0), r, ok ? "OK" : "FAIL");
        CU_ASSERT(ok)
    }

    {
        INT r = Y(log2i)(-1);
        int ok = r == -1;
        printf("log2i(%td) = %td -> %s\n", (INT)(-1), r, ok ? "OK" : "FAIL");
        CU_ASSERT(ok)
    }

    for (i = 0, j = 1; i < 8 * SIZEOF_PTRDIFF_T - 1; i++, j <<= 1)
    {
        {
            INT r = Y(log2i)(j);
            INT r2 = _log2i(j);
            int ok = r == r2;
            printf("log2i(%td) = %td -> %s\n", j, r, ok ? "OK" : "FAIL");
            CU_ASSERT(ok)
        }
        {
            INT r = Y(log2i)(j - 1);
            INT r2 = _log2i(j - 1);
            int ok = r == r2;
            printf("log2i(%td) = %td -> %s\n", j - 1, r, ok ? "OK" : "FAIL");
            CU_ASSERT(ok)
        }
    }
}

