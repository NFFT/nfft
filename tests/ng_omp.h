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

#ifndef NG_OMP_TEST_H
#define NG_OMP_TEST_H

#include "infft.h"

/* The add-on library's roster contract now that it carries a threaded solver. */
void Y(check_ng_omp_roster)(void);

/* The threaded DECONV solver agrees with the serial solver of the same rank. */
void Y(check_ng_omp_deconv_1d)(void);

/* A one-thread budget makes the threaded solver decline. */
void Y(check_ng_omp_declines_serial_budget)(void);

#endif /* NG_OMP_TEST_H */
