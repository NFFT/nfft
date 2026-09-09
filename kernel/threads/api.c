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

/* The add-on threading library's entry points, after FFTW's threads/api.c.
 *
 * Asking for threads is what installs them: X(plan_with_nthreads) initialises
 * threading first, and initialising registers the threaded roster. FFTW relies
 * on the same order, which is why its NO_NONTHREADEDP needs no third condition.
 * Registering solvers changes the configuration signature and so invalidates
 * every stored decision, so initialisation destroys the planner first, as
 * FFTW's X(cleanup)() does. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

#include <fftw3.h>

static int threads_inited = 0;

int Y(init_threads)(void)
{
  if (threads_inited)
    return 1;

  /* Nothing to register yet: report failure rather than leave a caller with a
   * raised thread count and no threaded solver to serve it. */
  if (Y(nfft_threads_roster_size)() == 0)
    return 0;

  Y(the_planner_destroy)(); /* the roster is about to change */
  Y(nfft_ensure_registered)();
  Y(nfft_threads_conf_standard)(Y(the_planner)());

  /* Let the main library observe FFTW's thread count for the wisdom key. */
  Y(fftw_nthreads_hook) = FFTW(planner_nthreads);

  threads_inited = 1;
  return 1;
}

void Y(cleanup_threads)(void)
{
  if (!threads_inited)
    return;
  Y(fftw_nthreads_hook) = 0;
  Y(the_planner_destroy)();
  threads_inited = 0;
}

/* The maximum number of threads a plan may use, not a target: a threaded solver
 * may use fewer and passes the rest of the budget to its children. The count is
 * part of the wisdom key, so changing it makes existing entries a clean miss.
 * Leaves the count at 1 while no threaded solver exists. */
void Y(plan_with_nthreads)(int nthreads)
{
  planner *pl;
  if (!Y(init_threads)())
    return;
  pl = Y(the_planner)();
  pl->nthr = nthreads < 1 ? 1 : nthreads;
}

int Y(planner_nthreads)(void)
{
  Y(nfft_ensure_registered)();
  return Y(the_planner)()->nthr;
}
