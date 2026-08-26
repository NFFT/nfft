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

/* Measurement of pln->adt->apply(pln, p), on two clocks: a fine clock times
 * the executions (raw ticks via ticks.h where a cycle counter exists, else
 * wall seconds), and a budget clock (wall seconds, Y(planner_clock_now))
 * bounds total measurement time.
 *
 * The iteration count n doubles from 1. Each level takes up to
 * PLNR_TIME_REPEAT batches of n back-to-back applies and keeps the minimum
 * batch reading: scheduler interruptions can only inflate a batch, never
 * shrink it, so the minimum is the faithful signal. A level is accepted once
 * that minimum reaches the mode's floor, returning minimum / n.
 *
 * Returns a strictly-positive reading or exactly -1.0; zero and negatives
 * never escape. Every degraded exit returns the best positive reading seen,
 * or -1.0. A whole level of zero readings at a sizable n, or n past the hard
 * cap, means a broken clock and returns -1.0.
 *
 * Never touches pln->awake_state; the caller owns wakefulness. */

#define PLNR_TIME_N_CAP (1L << 30) /* hard cap on n; past this the clock is broken */
#define PLNR_TIME_N_ZERO (1L << 16) /* all-zero level at this n => frozen clock */

/* The budget clock. Routing every budget read through this one wrapper keeps
 * plan_measure_cost and the measured-race timelimit check on the same
 * underlying clock by construction. */
double Y(planner_clock_now)(void) {
  return (double)Y(clock_gettime_seconds)();
}

double Y(planner_elapsed_seconds)(double since) {
  return Y(planner_clock_now)() - since;
}

double Y(plan_measure_cost)(plan *pln, const problem *p) {
#if defined(HAVE_TICK_COUNTER)
  /* Tick mode: fine clock = raw ticks, budget clock = wall seconds.
   * getticks()/elapsed() must stay inside this arm. On counter-less platforms
   * ticks.h defines getticks() to 0U, which compiles only as a whole
   * statement, not inside a shared expression. */
  double budget_start;
  long n;
  double best_ratio = -1.0; /* best strictly-positive (min-batch / n) seen */

  budget_start = Y(planner_clock_now)();

  for (n = 1;; n *= 2) {
    double batch_min = -1.0; /* minimum fine-clock batch reading at this n */
    int rep;

    for (rep = 0; rep < PLNR_TIME_REPEAT; rep++) {
      ticks t_batch_start, t_batch_end;
      double batch_ticks;
      long k;

      t_batch_start = getticks();
      for (k = 0; k < n; k++)
        pln->adt->apply(pln, p);
      t_batch_end = getticks();

      batch_ticks = elapsed(t_batch_end, t_batch_start);

      if (batch_min < 0.0 || batch_ticks < batch_min)
        batch_min = batch_ticks;

      if ((Y(planner_clock_now)() - budget_start) >= (double)PLNR_TIME_LIMIT_SECONDS) {
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

    if (batch_min >= PLNR_TIME_MIN_TICKS)
      return batch_min / (double)n;

    /* Check before doubling: n must not overflow a 32-bit long. */
    if (n > (PLNR_TIME_N_CAP / 2L))
      return -1.0;
  }
  /* unreachable */
#elif defined(HAVE_CLOCK_GETTIME)
  /* Slow-timer fallback: wall seconds as both fine and budget clock; the
   * floor is PLNR_TIME_MIN_SLOW_SECONDS. */
  double budget_start;
  long n;
  double best_ratio = -1.0;

  budget_start = Y(planner_clock_now)();

  for (n = 1;; n *= 2) {
    double batch_min = -1.0;
    int rep;

    for (rep = 0; rep < PLNR_TIME_REPEAT; rep++) {
      double t_batch_start, t_batch_end;
      double batch_sec;
      long k;

      t_batch_start = Y(planner_clock_now)();
      for (k = 0; k < n; k++)
        pln->adt->apply(pln, p);
      t_batch_end = Y(planner_clock_now)();

      batch_sec = t_batch_end - t_batch_start;

      if (batch_min < 0.0 || batch_sec < batch_min)
        batch_min = batch_sec;

      if ((t_batch_end - budget_start) >= (double)PLNR_TIME_LIMIT_SECONDS) {
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

    if (batch_min >= PLNR_TIME_MIN_SLOW_SECONDS)
      return batch_min / (double)n;

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
