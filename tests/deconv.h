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

#ifndef DECONV_TEST_H
#define DECONV_TEST_H

#include "infft.h"

/* The DECONV problem type carries exactly the attributes its stage needs and
 * omits the ones that must be absent; the hash is stable and discriminating. */
void Y(check_deconv_problem)(void);
/* The DECONV solver on a 1D KB problem, asserted as values: forward divides a
 * zero-padded grid by phi_hut; the adjoint multiplies by phi_hut_inv and
 * gathers.  Also covers the type-II +1 frequency shift. */
void Y(check_deconv_solver)(void);
/* Index mapping and scaling per rank, swept over the axis-type space (even
 * type-I, even type-II, odd), where the frequency layout differs.  Rank 1 also
 * covers the n < N decline. */
void Y(check_deconv_1d)(void);
void Y(check_deconv_2d)(void);
void Y(check_deconv_3d)(void);
void Y(check_deconv_nd)(void);

#endif
