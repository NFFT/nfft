/* $Id$
 *
 * Copyright (c) 2007, 2008 Jens Keiner, Stefan Kunis, Daniel Potts
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
#include <stdio.h>
#include <stdlib.h>

nfft_malloc_type_function nfft_malloc_hook = 0;
nfft_free_type_function nfft_free_hook = 0;
nfft_die_type_function nfft_die_hook = 0;

void *nfft_malloc(size_t n)
{
  void *p;

  if (nfft_malloc_hook)
    return nfft_malloc_hook(n);

  if (n == 0)
    n = 1;

  p = fftw_malloc(n);

  if (!p)
    nfft_die("nfft_malloc: out of memory (n = %d)\n",n);

  return p;
}

void nfft_free(void *p)
{
  if (p)
  {
        if (nfft_free_hook)
          {
      nfft_free_hook(p);
      return;
          }

    fftw_free(p);
  }
}

void nfft_die(const char *s)
{
  if (nfft_die_hook)
    nfft_die_hook(s);

  fflush(stdout);
  fprintf(stderr, "nfft: %s", s);
  exit(EXIT_FAILURE);
}
