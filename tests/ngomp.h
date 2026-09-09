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

#ifndef NGOMP_TEST_H
#define NGOMP_TEST_H

#include "infft.h"

/* The add-on library's contract while its roster is empty: no threaded solver
 * to reach, no way to raise the thread count. */
void Y(check_ngomp_empty_roster)(void);

/* Planning still works at every patience level through the add-on library. */
void Y(check_ngomp_plans_serially)(void);

#endif /* NGOMP_TEST_H */
