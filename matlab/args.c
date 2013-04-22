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

/* $Id$ */

#include "config.h"

#ifdef HAVE_UCHAR_H
  #include <uchar.h>
#endif

#include "imex.h"

int nfft_mex_get_int(const mxArray *p, const char *errmsg)
{
  DM(if (!mxIsDouble(p) || mxIsComplex(p) || mxGetM(p) != 1 || mxGetN(p) != 1)
    mexErrMsgTxt(errmsg);)
  return mxGetScalar(p);
}

double nfft_mex_get_double(const mxArray *p, const char *errmsg)
{
  DM(if (!mxIsDouble(p) || mxIsComplex(p) || mxGetM(p) != 1 || mxGetN(p) != 1)
    mexErrMsgTxt(errmsg);)
  return mxGetScalar(p);
}
