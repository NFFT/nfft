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

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"
#include "reflect.h"

void X(check_get_version)(void)
{
    unsigned major, minor, patch;
    char v1[20], v2[20];

    Y(get_version)(&major, &minor, &patch);

    snprintf(v1, sizeof(v1), "%u.%u.%u", major, minor, patch);
    snprintf(v2, sizeof(v2), "%u.%u.%u", NFFT_VERSION_MAJOR, NFFT_VERSION_MINOR, NFFT_VERSION_PATCH);

    CU_ASSERT(strncmp(v1, v2, sizeof(v1)) == 0);
}

void X(check_get_window_name)(void)
{
    const char *window_name1 = STRINGIZE(WINDOW_NAME);
    const char *window_name2 = Y(get_window_name)();

    /* Just compare what's returned by the method against the literal from config.h. */
    CU_ASSERT(strncmp(window_name1, window_name2, sizeof(window_name1)) == 0);
}
