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

#ifndef CONV_TEST_H
#define CONV_TEST_H

#include "infft.h"

/* The CONV problem type carries exactly the attributes its stage needs, and is
 * data-blind in x. */
void Y(check_conv_problem)(void);
/* The CONV solver on a 1D KB problem, asserted as values: forward gathers the
 * grid against psi with the wrap formula; the adjoint scatters a node spike
 * onto its wrapped neighbours. */
void Y(check_conv_solver)(void);
/* The same gather/scatter identities per rank, swept over the axis-type space.
 * The rank-2/3/4 leaves have their own index arithmetic and were previously
 * exercised only indirectly, through the composed solver. */
void Y(check_conv_1d)(void);
void Y(check_conv_2d)(void);
void Y(check_conv_3d)(void);
void Y(check_conv_nd)(void);

#endif
