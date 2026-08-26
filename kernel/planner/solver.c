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

/* Refcounted solver base type: a stateless recipe. Callers embed it as the
 * first member of a larger struct and pass the total size here. */

solver *Y(solver_create)(size_t size, const solver_adt *adt) {
  solver *ego = (solver *)Y(malloc)(size); /* Y(malloc) dies on OOM */
  A(size >= sizeof(solver));
  ego->adt = adt;
  ego->refcnt = 0;
  return ego;
}

void Y(solver_use)(solver *ego) {
  A(ego != 0);
  ego->refcnt++;
}

void Y(solver_destroy)(solver *ego) {
  if (ego == 0)
    return;
  A(ego->refcnt > 0);
  if (--ego->refcnt == 0) {
    if (ego->adt->destroy != 0)
      ego->adt->destroy(ego);
    Y(free)
    (ego);
  }
}
