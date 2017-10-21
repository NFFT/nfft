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

#include<stdlib.h>

#include "api.h"

Y(malloc_type_function) Y(malloc_hook) = 0;
Y(free_type_function) Y(free_hook) = 0;
Y(die_type_function) Y(die_hook) = 0;

void *Y(malloc)(size_t n)
{
  void *p;

  if (Y(malloc_hook))
    return Y(malloc_hook)(n);

  if (n == 0)
    n = 1;

  p = FFTW(malloc)(n);

  if (!p)
    Y(die)(STRINGIZE(Y(malloc)) ": out of memory\n");

  return p;
}

void Y(free)(void *p)
{
  if (p)
  {
    if (Y(free_hook))
    {
      Y(free_hook)(p);
      return;
    }
    FFTW(free)(p);
  }
}

void Y(die)(const char *s)
{
  if (Y(die_hook))
    Y(die_hook)(s);

  exit(EXIT_FAILURE);
}
