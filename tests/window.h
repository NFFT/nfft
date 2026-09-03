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

void X(check_bessel_i0_scaling)(void);
void X(check_kaiser_bessel_peak)(void);
void X(check_kaiser_bessel_reference)(void);
void X(check_kaiser_bessel_cancellation)(void);
void X(check_kaiser_bessel_nfft)(void);
void X(check_kaiser_bessel_run)(void);
void X(check_kaiser_bessel_phi)(void);
void X(check_log_sinc)(void);
void X(check_log_sinc_exp)(void);
void X(check_bspline_run)(void);
void X(check_bspline_phi_hut_reference)(void);
void X(check_bspline_phi_reference)(void);
void X(check_sincpow_phi_reference)(void);
void X(check_sincpow_phi_hut_reference)(void);
void X(check_bspline_cheb)(void);
