/*
 * Copyright (c) 2006, 2010, 2016 The NFFT Authors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see <https://www.gnu.org/licenses/>.
 */

#include "api.h"

void Y(get_version)(unsigned *major, unsigned *minor, unsigned *patch)
{
    major[0] = NFFT_VERSION_MAJOR;
    minor[0] = NFFT_VERSION_MINOR;
    patch[0] = NFFT_VERSION_PATCH;
}
