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

/* CONV solver registry: the roster tying together the rank-gated leaf solvers
 * (kernel/conv/conv-1d.c for d==1, conv-2d.c for d==2, conv-3d.c for d==3,
 * conv-nd.c for the generic d>=4) and the generation-guarded one-shot
 * registration entry point, mirroring kernel/deconv/solver.c's deconv_roster.
 * The leaf register functions are all declared in kernel/conv/conv.h. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

static const solvtab conv_roster = {SOLVTAB(Y(conv_solver_1d_register)),
                                    SOLVTAB(Y(conv_solver_2d_register)),
                                    SOLVTAB(Y(conv_solver_3d_register)),
                                    SOLVTAB(Y(conv_solver_nd_register)),
                                    SOLVTAB_END};

void Y(conv_solvers_register)(planner *pl) {
  Y(solvtab_exec)
  (conv_roster, pl);
}

void Y(conv_ensure_registered)(void) {
  static unsigned gen = 0;
  planner *pl = Y(the_planner)();
  unsigned g = Y(the_planner_generation)();
  if (g != gen) {
    Y(conv_solvers_register)
    (pl);
    gen = g;
  }
}
