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
#include "deconv.h"

static const solvtab deconv_roster = {SOLVTAB(Y(deconv_solver_1d_register)),
                                      SOLVTAB(Y(deconv_solver_2d_register)),
                                      SOLVTAB(Y(deconv_solver_3d_register)),
                                      SOLVTAB(Y(deconv_solver_nd_register)),
                                      SOLVTAB_END};

void Y(deconv_solvers_register)(planner *pl) {
  Y(solvtab_exec)
  (deconv_roster, pl);
}

void Y(deconv_ensure_registered)(void) {
  static unsigned gen = 0;
  planner *pl = Y(the_planner)();
  unsigned g = Y(the_planner_generation)();
  if (g != gen) {
    Y(deconv_solvers_register)
    (pl);
    gen = g;
  }
}
