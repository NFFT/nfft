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

#ifdef HAVE_TIME_H
# include <time.h>
#endif

/* Measurement of pln->adt->apply(pln, p), on two clocks: a fine clock times
 * the executions (raw ticks via ticks.h where a cycle counter exists, else
 * wall seconds), and a budget clock (wall seconds, Y(planner_clock_now))
 * bounds total measurement time.
 *
 * The iteration count n doubles from 1. Each level takes up to
 * PLNR_TIME_REPEAT batches of n back-to-back applies and keeps the minimum
 * batch reading: scheduler interruptions can only inflate a batch, never
 * shrink it, so the minimum is the faithful signal. A level is accepted once
 * that minimum reaches the fine clock's floor, returning minimum / n.
 *
 * Returns a strictly-positive reading or exactly -1.0; zero and negatives
 * never escape. Every degraded exit returns the best positive reading seen,
 * or -1.0. A whole level of zero readings at a sizable n, or n past the hard
 * cap, means a broken clock and returns -1.0.
 *
 * Never touches pln->awake_state; the caller owns wakefulness. */

#define PLNR_TIME_N_CAP                                                       \
  (1L << 30) /* hard cap on n; past this the clock is broken */
#define PLNR_TIME_N_ZERO                                                      \
  (1L << 16) /* all-zero level at this n => frozen clock */

#if defined(HAVE_TICK_COUNTER)
/* cycle.h sets TIME_MIN per counter where short reads are noise (5000 ticks
 * for rdtsc); 100.0 is the fallback where it defines none. */
# ifndef TIME_MIN
#  define TIME_MIN 100.0
# endif
typedef ticks fine_t;
# define FINE_NOW() getticks()
# define FINE_DIFF(t1, t0) elapsed(t1, t0)
# define FINE_FLOOR TIME_MIN
#elif defined(HAVE_CLOCK_GETTIME)
typedef double fine_t;
# define FINE_NOW() Y(planner_clock_now)()
# define FINE_DIFF(t1, t0) ((t1) - (t0))
# define FINE_FLOOR PLNR_TIME_MIN_SLOW_SECONDS
#endif

/* The budget clock. Routing every budget read through this one wrapper keeps
 * plan_measure_cost and the measured-race timelimit check on the same
 * underlying clock. The result is a double whatever the build precision is:
 * in the float build an epoch near 1.7e9 has a 128 s ulp, which would starve
 * every budget check. */
double Y(planner_clock_now)(void)
{
#if defined(HAVE_CLOCK_GETTIME)
# ifdef CLOCK_MONOTONIC
#  define PLNR_CLOCK_ID CLOCK_MONOTONIC
# else
#  define PLNR_CLOCK_ID CLOCK_REALTIME
# endif
  struct timespec ts;
  if (clock_gettime(PLNR_CLOCK_ID, &ts) == 0)
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
  return 0.0;
}

double Y(planner_elapsed_seconds)(double since)
{
  return Y(planner_clock_now)() - since;
}

double Y(plan_measure_cost)(plan *pln, const problem *p)
{
#if defined(FINE_NOW)
  double budget_start;
  long n;
  double best_ratio = -1.0; /* best strictly-positive (min-batch / n) seen */

  budget_start = Y(planner_clock_now)();

  for (n = 1;; n *= 2) {
    double batch_min = -1.0; /* minimum fine-clock batch reading at this n */
    int rep;

    for (rep = 0; rep < PLNR_TIME_REPEAT; rep++) {
      fine_t t_batch_start, t_batch_end;
      double batch;
      long k;

      t_batch_start = FINE_NOW();
      for (k = 0; k < n; k++)
        pln->adt->apply(pln, p);
      t_batch_end = FINE_NOW();

      batch = FINE_DIFF(t_batch_end, t_batch_start);

      if (batch_min < 0.0 || batch < batch_min)
        batch_min = batch;

      if (Y(planner_elapsed_seconds)(budget_start)
          >= (double)PLNR_TIME_LIMIT_SECONDS) {
        if (batch_min > 0.0) {
          double ratio = batch_min / (double)n;
          if (best_ratio < 0.0 || ratio < best_ratio)
            best_ratio = ratio;
        }
        return best_ratio;
      }
    }

    if (batch_min <= 0.0 && n >= PLNR_TIME_N_ZERO)
      return -1.0;

    if (batch_min > 0.0) {
      double ratio = batch_min / (double)n;
      if (best_ratio < 0.0 || ratio < best_ratio)
        best_ratio = ratio;
    }

    if (batch_min >= FINE_FLOOR)
      return batch_min / (double)n;

    /* Check before doubling: n must not overflow a 32-bit long. */
    if (n > (PLNR_TIME_N_CAP / 2L))
      return -1.0;
  }
  /* unreachable */
#else
  UNUSED(pln);
  UNUSED(p);
  return -1.0;
#endif
}
