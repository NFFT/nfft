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

/* Shared infrastructure for Step A (deconvolution) of the fast NFFT decomposition.
 * Include after nfft3.h, infft.h and iplanner.h. */

#ifndef NFFT_DECONV_H
#define NFFT_DECONV_H

/* Analytical cost of DECONV problem. */
double Y(deconv_d_pcost)(const problem *p);

/* DECONV solvers */
void Y(deconv_solver_1d_register)(planner *pl);
void Y(deconv_solver_2d_register)(planner *pl);
void Y(deconv_solver_3d_register)(planner *pl);
void Y(deconv_solver_nd_register)(planner *pl);

#endif /* NFFT_DECONV_H */
