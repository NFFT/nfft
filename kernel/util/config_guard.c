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

/* The build reads its own config.h, checked at compile time.
 *
 * NFFT_CONFIG_GUARD_PRECISION is the precision the build system configured;
 * infft.h reflects the precision of the config.h it actually read. They
 * disagree when one config.h shadows another, which is what infft.h's
 * `#include <config.h>` and the generated directory's place at the head of
 * the include path exist to prevent. This is the check that they still do.
 *
 * An in-tree Autotools configure writes include/config.h into the source
 * tree, beside infft.h, so a CMake build of the same tree has two candidates
 * on its include path and the order decides. Builds that leave the macro
 * unset compile this file to nothing.
 */

#include "infft.h"

#ifdef NFFT_CONFIG_GUARD_PRECISION

#define NFFT_CONFIG_GUARD_SINGLE 1
#define NFFT_CONFIG_GUARD_DOUBLE 2
#define NFFT_CONFIG_GUARD_LDOUBLE 3

#if defined(NFFT_SINGLE)
#define NFFT_CONFIG_GUARD_SEEN NFFT_CONFIG_GUARD_SINGLE
#elif defined(NFFT_LDOUBLE)
#define NFFT_CONFIG_GUARD_SEEN NFFT_CONFIG_GUARD_LDOUBLE
#else
#define NFFT_CONFIG_GUARD_SEEN NFFT_CONFIG_GUARD_DOUBLE
#endif

#if NFFT_CONFIG_GUARD_SEEN != NFFT_CONFIG_GUARD_PRECISION
#error "config.h precision does not match the configured one. A generated include/config.h in the source tree shadows this build's own; delete it (or run 'make distclean') and build again."
#endif

#endif /* NFFT_CONFIG_GUARD_PRECISION */

typedef int nfft_config_guard_dummy;
