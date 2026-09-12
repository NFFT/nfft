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

/* opcnt arithmetic, after FFTW's kernel/ops.c. */

void Y(ops_zero)(opcnt *dst)
{
  dst->add = dst->mul = dst->fma = dst->other = 0.0;
}

void Y(ops_other)(INT o, opcnt *dst)
{
  Y(ops_zero)(dst);
  dst->other = (double)o;
}

void Y(ops_cpy)(const opcnt *src, opcnt *dst)
{
  *dst = *src;
}

void Y(ops_madd)(INT m, const opcnt *a, const opcnt *b, opcnt *dst)
{
  dst->add = (double)m * a->add + b->add;
  dst->mul = (double)m * a->mul + b->mul;
  dst->fma = (double)m * a->fma + b->fma;
  dst->other = (double)m * a->other + b->other;
}

void Y(ops_add)(const opcnt *a, const opcnt *b, opcnt *dst)
{
  Y(ops_madd)(1, a, b, dst);
}

void Y(ops_add2)(const opcnt *a, opcnt *dst)
{
  Y(ops_add)(a, dst, dst);
}

void Y(ops_madd2)(INT m, const opcnt *a, opcnt *dst)
{
  Y(ops_madd)(m, a, dst, dst);
}
