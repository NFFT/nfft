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

/* Problem base type: a hashable description of work ("what"). Concrete
 * problems embed this struct as their first member and pass the total size
 * here; the adt supplies the hash/print/destroy hooks. */

problem *Y(problem_create)(size_t size, const problem_adt *adt) {
  problem *p = (problem *)Y(malloc)(size);
  A(size >= sizeof(problem));
  p->adt = adt;
  return p;
}

/* adt->destroy frees the problem's owned members, then the base is freed here.
 * It may be NULL for a problem that owns nothing beyond the base. */
void Y(problem_destroy)(problem *p) {
  A(p != 0);
  if (p->adt->destroy != 0)
    p->adt->destroy(p);
  Y(free)
  (p);
}
