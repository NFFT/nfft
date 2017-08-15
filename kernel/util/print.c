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

#include <stdio.h>
#include "infft.h"

/** Print real vector to standard output. */
void Y(vpr_double)(R *x, const INT n, const char *text)
{
  INT k;

  if (x == NULL)
  {
    printf("null pointer\n");
    fflush(stdout);
    exit(-1);
  }

  if (text != NULL)
  {
    printf ("\n %s, adr=%p\n", text, (void*)x);

    for (k = 0; k < n; k++)
    {
      if (k%8 == 0)
        printf("%6td.\t", k);

      printf("%+.1" __FES__ ",", x[k]);

      if (k%8 == 7)
        printf("\n");
    }

    if (n%8 != 0)
      printf("\n");
  }
  else
    for (k = 0; k < n; k++)
      printf("%+" __FES__ ",\n", x[k]);

  fflush(stdout);
}

/** Print complex vector to standard output. */
void Y(vpr_complex)(C *x, const INT n, const char *text)
{
  INT k;

  if(text != NULL)
  {
    printf("\n %s, adr=%p\n", text, (void*)x);
    for (k = 0; k < n; k++)
    {
      if (k%4 == 0)
        printf("%6td.\t", k);

      printf("%+.1" __FES__ "%+.1" __FES__ "i,", CREAL(x[k]), CIMAG(x[k]));

      if (k%4==3)
        printf("\n");
    }
    if (n%4!=0)
      printf("\n");
  }
  else
    for (k = 0; k < n; k++)
      printf("%+" __FES__ "%+" __FES__ "i,\n", CREAL(x[k]), CIMAG(x[k]));

  fflush(stdout);
}
