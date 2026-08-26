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

/* Plan base type: an executable produced by a solver ("how"). Concrete plans
 * embed this struct as their first member and pass the total size here.
 * Wakefulness is SLEEPY (no tables) or AWAKE (precomputed); the
 * adt awake hook builds and releases those tables. */

plan *Y(plan_create)(size_t size, const plan_adt *adt) {
  plan *ego = (plan *)Y(malloc)(size);
  A(size >= sizeof(plan));
  ego->adt = adt;
  ego->pcost = 0.0;
  ego->awake_state = PLNR_SLEEPY;
  ego->uses_x = 1;
  return ego;
}

/* Idempotent awake hook */
void Y(plan_awake)(plan *ego, int wakefulness) {
  A(ego != 0);
  if (ego->awake_state == wakefulness)
    return;
  if (ego->adt->awake != 0)
    ego->adt->awake(ego, wakefulness);
  ego->awake_state = wakefulness;
}

/* Sleep first so the awake hook releases tables exactly once,
 * then adt->destroy frees the solver-specific state, then the base. */
void Y(plan_destroy)(plan *ego) {
  A(ego != 0);
  Y(plan_awake)
  (ego, PLNR_SLEEPY);
  if (ego->adt->destroy != 0)
    ego->adt->destroy(ego);
  Y(free)
  (ego);
}
