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

/* $Id: util.c 3483 2010-04-23 19:02:34Z keiner $ */

#include<stdlib.h>

#include "api.h"

X(malloc_type_function) X(malloc_hook) = 0;
X(free_type_function) X(free_hook) = 0;
X(die_type_function) X(die_hook) = 0;

void *X(malloc)(size_t n)
{
  void *p;

  if (X(malloc_hook))
    return nfft_malloc_hook(n);

  if (n == 0)
    n = 1;

  p = Z(malloc)(n);

  X(die)(STRINGIZE(X(malloc)) ": out of memory\n");

  return p;
}

void X(free)(void *p)
{
  if (p)
  {
    if (X(free_hook))
    {
      X(free_hook)(p);
      return;
    }
    Z(free)(p);
  }
}

void X(die)(const char *s)
{
  if (X(die_hook))
    X(die_hook)(s);

  exit(EXIT_FAILURE);
}
