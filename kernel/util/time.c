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

#ifdef HAVE_TIME_H
#include <time.h>
#endif

R Y(elapsed_seconds)(ticks t1, ticks t0)
{
  UNUSED(t1);
  UNUSED(t0);
  return (R)(elapsed(t1,t0)) / (R)(TICKS_PER_SECOND);
}

R Y(clock_gettime_seconds)(void)
{
#if defined(HAVE_CLOCK_GETTIME)
  struct timespec tp;
  if (clock_gettime(CLOCK_REALTIME, &tp) != 0)
  {
    tp.tv_sec = 0;
    tp.tv_nsec = 0;
  }

  return tp.tv_sec+(R)tp.tv_nsec/K(1e9);
#else
  return K(0.0);
#endif
}

