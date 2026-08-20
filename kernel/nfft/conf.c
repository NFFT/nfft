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

/* NFFT solver roster (~ FFTW's dft/conf.c).  Y(nfft_solver_fast_native) is the
 * sole fast solver and claims every d, composed from independently-planned
 * DECONV/CONV sub-problems; the direct NDFT is two 1D solvers (plain per-term,
 * blocked recurrence) plus one generic solver for d >= 2.  Applicability is
 * expressed by mkplan returning NULL.  PRE_PSI is the only fast psi strategy
 * this planner offers.
 *
 * Registration order matters for determinism: iteration is reverse
 * registration order, and on exact pcost ties the earlier-encountered plan is
 * kept. */
static const solvtab the_roster = {
    SOLVTAB(Y(nfft_solver_fast_native_register)),
    SOLVTAB(Y(nfft_solver_ndft_1d_register)),
    SOLVTAB(Y(nfft_solver_ndft_1d_blocked_register)),
    SOLVTAB(Y(nfft_solver_ndft_nd_register)),
    SOLVTAB(Y(nfft_solver_rnk0_register)),
    SOLVTAB_END};

void Y(nfft_solvers_register)(planner *pl) {
  Y(solvtab_exec)
  (the_roster, pl);
}

/* Lazy, idempotent registration into the process-global planner, once per
 * global-planner generation (Y(the_planner_generation), theplanner.c).  Not
 * thread-safe by contract: planning is single-threaded.
 *
 * DECONV/CONV are ensured here, eagerly and before the NFFT roster, because
 * mkplan_native_fast (nfft-nd.c) plans its children by recursing into
 * Y(planner_mkplan)().  Registering them from inside that recursion would
 * realloc pl->slvdescs while FORALL_SOLVERS_OF_KIND still holds a raw pointer
 * into the old array -- a use-after-free that surfaces only on a later
 * planner call, once the freed block has been reused. */
void Y(nfft_ensure_registered)(void) {
  static unsigned registered_gen = 0;
  planner *pl = Y(the_planner)();
  unsigned gen = Y(the_planner_generation)();
  if (gen != registered_gen) {
    Y(deconv_ensure_registered)
    ();
    Y(conv_ensure_registered)
    ();
    Y(nfft_solvers_register)
    (pl);
    registered_gen = gen;
  }
}
