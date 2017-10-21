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

/**
 * Prints an error message for a failed assertion together with filename and the
 * line where the assertion failed.
 */
void Y(assertion_failed)(const char *s, int line, const char *file)
{
  fflush(stdout);
  fprintf(stderr, "nfft: %s:%d: assertion failed: %s\n", file, line, s);
#ifdef HAVE_ABORT
  /* Use abort function. */
  abort();
#else
  /* Use exit function. */
  exit(EXIT_FAILURE);
#endif
}
