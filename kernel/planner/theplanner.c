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

/* The process-global planner, created lazily on first use. It registers no
 * solvers itself: problem-kind modules register their own, which keeps
 * kernel/planner/ free of transform dependencies. Planning through this
 * instance is not thread-safe; executing distinct plans is.
 *
 * The generation counter increments on each lazy create. Registration hooks
 * must compare generations, never planner pointers: a freed-and-recreated
 * planner can reuse the same address (malloc ABA), and a pointer-identity
 * check would then silently skip re-registration. */

static planner *the_planner_global = 0;
static unsigned the_planner_gen = 0; /* 0 = never created */

planner *Y(the_planner)(void) {
  if (the_planner_global == 0) {
    the_planner_global = Y(planner_create)();
    ++the_planner_gen;
  }
  return the_planner_global;
}

unsigned Y(the_planner_generation)(void) {
  return the_planner_gen;
}

/* Safe to call when absent. The next the_planner() call recreates and bumps
 * the generation. */
void Y(the_planner_destroy)(void) {
  if (the_planner_global != 0) {
    Y(planner_destroy)
    (the_planner_global);
    the_planner_global = 0;
  }
}
