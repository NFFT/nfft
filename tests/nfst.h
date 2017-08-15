/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "infft.h"

#undef X
#define X(name) NFST(name)

void X(check_1d_direct_file)(void);
void X(check_1d_fast_file)(void);
void X(check_1d_online)(void);
void X(check_2d_direct_file)(void);
void X(check_2d_fast_file)(void);
void X(check_2d_online)(void);
void X(check_3d_direct_file)(void);
void X(check_3d_fast_file)(void);
void X(check_3d_online)(void);
void X(check_4d_online)(void);

void X(check_adjoint_1d_direct_file)(void);
void X(check_adjoint_1d_fast_file)(void);
void X(check_adjoint_1d_online)(void);
void X(check_adjoint_2d_direct_file)(void);
void X(check_adjoint_2d_fast_file)(void);
void X(check_adjoint_2d_online)(void);
void X(check_adjoint_3d_direct_file)(void);
void X(check_adjoint_3d_fast_file)(void);
void X(check_adjoint_3d_online)(void);
void X(check_adjoint_4d_online)(void);
