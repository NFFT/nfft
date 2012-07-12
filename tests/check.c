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

/* Standard headers. */
#include <CUnit/CUnit.h>
#include <CUnit/TestRun.h>

#include "infft.h"
#include "bspline.h"
#include "bessel.h"
#include "nfft.h"

int main(void)
{
  CU_pSuite util, nfft;
  CU_initialize_registry();
  /*CU_set_output_filename("nfft");*/
  nfft = CU_add_suite("nfft", 0, 0);
  CU_add_test(nfft, "nfft_1d_file", X(check_nfft_1d_file));
  CU_add_test(nfft, "nfft_1d_online", X(check_nfft_1d_online));
  CU_add_test(nfft, "nfft_2d_file", X(check_nfft_2d_file));
  CU_add_test(nfft, "nfft_2d_online", X(check_nfft_2d_online));
  CU_add_test(nfft, "nfft_3d_file", X(check_nfft_3d_file));
  CU_add_test(nfft, "nfft_3d_online", X(check_nfft_3d_online));
  CU_add_test(nfft, "nfft_4d_online", X(check_nfft_4d_online));
  util = CU_add_suite("util", 0, 0);
  CU_add_test(util, "bspline", X(check_bspline));
  CU_add_test(util, "bessel_i0", X(check_bessel_i0));
  CU_automated_run_tests();
  /*CU_basic_run_tests();*/
  {
    unsigned int ok = (CU_get_number_of_tests_failed() == 0);
    CU_cleanup_registry();
    return IF(ok, EXIT_SUCCESS, EXIT_FAILURE);
  }
}
