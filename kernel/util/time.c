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

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif

#ifdef HAVE_MACH_MACH_TIME_H
#include <mach/mach_time.h>
#endif

R Y(elapsed_seconds)(ticks t1, ticks t0)
{
  UNUSED(t1);
  UNUSED(t0);
  return (R)(elapsed(t1,t0)) / (R)(TICKS_PER_SECOND);
}

/* Monotonic wall clock, so an NTP step cannot corrupt an interval. The ladder
 * mirrors FFTW's timer.c: prefer a monotonic source, fall back to a
 * wall-clock one only where no monotonic source exists. The origin is
 * arbitrary and differs per branch; only differences are meaningful. */
double Y(clock_gettime_seconds)(void)
{
#if defined(HAVE_CLOCK_GETTIME) && defined(CLOCK_MONOTONIC)
  struct timespec tp;
  if (clock_gettime(CLOCK_MONOTONIC, &tp) != 0)
    return 0.0;
  return (double)tp.tv_sec + (double)tp.tv_nsec / 1e9;
#elif defined(_WIN32) || defined(_WIN64)
  LARGE_INTEGER t, f;
  if (!QueryPerformanceFrequency(&f) || f.QuadPart == 0
      || !QueryPerformanceCounter(&t))
    return 0.0;
  return (double)t.QuadPart / (double)f.QuadPart;
#elif defined(HAVE_MACH_ABSOLUTE_TIME)
  static mach_timebase_info_data_t tb;
  if (tb.denom == 0 && mach_timebase_info(&tb) != 0)
    return 0.0;
  return (double)mach_absolute_time() * (double)tb.numer / (double)tb.denom / 1e9;
#elif defined(HAVE_GETTIMEOFDAY)
  struct timeval tv;
  if (gettimeofday(&tv, 0) != 0)
    return 0.0;
  return (double)tv.tv_sec + (double)tv.tv_usec / 1e6;
#else
  return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

