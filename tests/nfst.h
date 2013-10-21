/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: simple_test.c 3509 2010-05-25 19:00:59Z keiner $ */

#include "infft.h"

/* Determine precision and name-mangling scheme. */
#define CONCAT(prefix, name) prefix ## name
#if defined(NFFT_SINGLE)
#define X(name) CONCAT(nfftf_,name)
#define Y(name) CONCAT(nfstf_,name)
#elif defined(NFFT_LDOUBLE)
#define X(name) CONCAT(nfftl_,name)
#define Y(name) CONCAT(nfstl_,name)
#else
#define X(name) CONCAT(nfft_,name)
#define Y(name) CONCAT(nfst_,name)
#endif

void X(check_nfst_1d_file)(void);
void X(check_nfst_1d_online)(void);
//void X(check_nfst_2d_file)(void);
//void X(check_nfst_2d_online)(void);
//void X(check_nfst_3d_file)(void);
//void X(check_nfst_3d_online)(void);
//void X(check_nfst_4d_online)(void);
