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

/* Execute a solver-registration table. Each entry names a registrar and a
 * function registering one or more solvers, which receive reg_ids 0, 1, ...
 * under that registrar name. Those two identifiers are what wisdom files
 * persist, so renaming a registrar invalidates existing wisdom. */
void Y(solvtab_exec)(const solvtab tbl, planner *pl) {
  int i;
  for (i = 0; tbl[i].reg != 0; i++) {
    pl->cur_reg_nam = tbl[i].reg_nam;
    pl->cur_reg_id = 0;
    tbl[i].reg(pl);
  }
  pl->cur_reg_nam = 0;
}
