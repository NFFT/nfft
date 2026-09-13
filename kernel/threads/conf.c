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

/* The threaded solver roster. This is the one place a threaded solver is
 * registered, and it exists only in the add-on library, so a program that does
 * not link that library has no threaded solvers and cannot raise the planner's
 * thread count. FFTW's threads/conf.c has the same job. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "deconv.h"
#include "conv.h"

static const solvtab s = {SOLVTAB(Y(deconv_solver_threaded_register)),
                          SOLVTAB(Y(conv_solver_threaded_register)),
                          SOLVTAB_END};

void Y(nfft_threads_conf_standard)(planner *pl)
{
  Y(solvtab_exec)(s, pl);
}

int Y(nfft_threads_roster_size)(void)
{
  int n = 0;
  while (s[n].reg != 0)
    n++;
  return n;
}
