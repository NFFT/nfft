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

#ifndef __API_H__
#define __API_H__

#ifndef CALLING_NFFT /* defined when calling internal functions. */
# ifndef COMPILING_NFFT
#  define COMPILING_NFFT /* used for dymbol exporting in nfft3.h */
# endif
#endif

/* When compiling with GNU libtool on Windows, DLL_EXPORT is #defined for the
   shared-library objects. In this case, we'll define NFFT_DLL to add dllexport 
   attributes to the specified functions in nfft3.h; see also api/api.h in FFTW. */
#ifdef DLL_EXPORT
# ifndef NFFT_DLL
#  define NFFT_DLL
# endif
#endif

#include "nfft3.h"
#include "infft.h"

#endif /* __API_H__ */
