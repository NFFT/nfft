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
#include <CUnit/CUnit.h>
#include <CUnit/TestRun.h>
#include <CUnit/Automated.h>

#include "infft.h"
#include "reflect.h"
#include "util.h"
#include "bspline.h"
#include "bessel.h"
#include "nfft.h"
#include "nfct.h"
#include "nfst.h"

int main(void)
{
  CU_pSuite util, nfft, nfct, nfst;
  CU_initialize_registry();
  /*CU_set_output_filename("nfft");*/
#ifdef _OPENMP
  CU_set_output_filename("CUnitAutomated_threads");
#endif

#undef X
#define X(name) Y(name)
  util = CU_add_suite("util", 0, 0);
  CU_add_test(util, "bspline", X(check_bspline));
  CU_add_test(util, "bessel_i0", X(check_bessel_i0));
  CU_add_test(util, "version", X(check_get_version));
  CU_add_test(util, "window_name", X(check_get_window_name));
  CU_add_test(util, "log2i", X(check_log2i));
  CU_add_test(util, "next_power_of_2", X(check_next_power_of_2));

#undef X
#define X(name) NFFT(name)
  nfft = CU_add_suite("nfft", 0, 0);
  CU_add_test(nfft, "nfft_1d_direct_file", X(check_1d_direct_file));
  CU_add_test(nfft, "nfft_1d_fast_file", X(check_1d_fast_file));
  CU_add_test(nfft, "nfft_adjoint_1d_direct_file", X(check_adjoint_1d_direct_file));
  CU_add_test(nfft, "nfft_adjoint_1d_fast_file", X(check_adjoint_1d_fast_file));
  CU_add_test(nfft, "nfft_1d_online", X(check_1d_online));
  CU_add_test(nfft, "nfft_adjoint_1d_online", X(check_adjoint_1d_online));

  CU_add_test(nfft, "nfft_2d_direct_file", X(check_2d_direct_file));
  CU_add_test(nfft, "nfft_2d_fast_file", X(check_2d_fast_file));
  CU_add_test(nfft, "nfft_adjoint_2d_direct_file", X(check_adjoint_2d_direct_file));
  CU_add_test(nfft, "nfft_adjoint_2d_fast_file", X(check_adjoint_2d_fast_file));
  CU_add_test(nfft, "nfft_2d_online", X(check_2d_online));
  CU_add_test(nfft, "nfft_adjoint_2d_online", X(check_adjoint_2d_online));

  CU_add_test(nfft, "nfft_3d_direct_file", X(check_3d_direct_file));
  CU_add_test(nfft, "nfft_3d_fast_file", X(check_3d_fast_file));
  CU_add_test(nfft, "nfft_adjoint_3d_direct_file", X(check_adjoint_3d_direct_file));
  CU_add_test(nfft, "nfft_adjoint_3d_fast_file", X(check_adjoint_3d_fast_file));
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  CU_add_test(nfft, "nfft_3d_online", X(check_3d_online));
  CU_add_test(nfft, "nfft_adjoint_3d_online", X(check_adjoint_3d_online));

  CU_add_test(nfft, "nfft_4d_online", X(check_4d_online));
  CU_add_test(nfft, "nfft_adjoint_4d_online", X(check_adjoint_4d_online));
#endif
#ifdef HAVE_NFCT
#undef X
#define X(name) NFCT(name)
#ifndef _OPENMP
  nfct = CU_add_suite("nfct", 0, 0);
  CU_add_test(nfct, "nfct_1d_direct_file", X(check_1d_direct_file));
  CU_add_test(nfct, "nfct_1d_fast_file", X(check_1d_fast_file));
  CU_add_test(nfct, "nfct_adjoint_1d_direct_file", X(check_adjoint_1d_direct_file));
  CU_add_test(nfct, "nfct_adjoint_1d_fast_file", X(check_adjoint_1d_fast_file));
  CU_add_test(nfct, "nfct_1d_online", X(check_1d_online));
  CU_add_test(nfct, "nfct_adjoint_1d_online", X(check_adjoint_1d_online));

  CU_add_test(nfct, "nfct_2d_direct_file", X(check_2d_direct_file));
  CU_add_test(nfct, "nfct_2d_fast_file", X(check_2d_fast_file));
  CU_add_test(nfct, "nfct_adjoint_2d_direct_file", X(check_adjoint_2d_direct_file));
  CU_add_test(nfct, "nfct_adjoint_2d_fast_file", X(check_adjoint_2d_fast_file));
  CU_add_test(nfct, "nfct_2d_online", X(check_2d_online));
  CU_add_test(nfct, "nfct_adjoint_2d_online", X(check_adjoint_2d_online));

  CU_add_test(nfct, "nfct_3d_direct_file", X(check_3d_direct_file));
  CU_add_test(nfct, "nfct_3d_fast_file", X(check_3d_fast_file));
  CU_add_test(nfct, "nfct_adjoint_3d_direct_file", X(check_adjoint_3d_direct_file));
  CU_add_test(nfct, "nfct_adjoint_3d_fast_file", X(check_adjoint_3d_fast_file));
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  CU_add_test(nfct, "nfct_3d_online", X(check_3d_online));
  CU_add_test(nfct, "nfct_adjoint_3d_online", X(check_adjoint_3d_online));

  CU_add_test(nfct, "nfct_4d_online", X(check_4d_online));
  CU_add_test(nfct, "nfct_adjoint_4d_online", X(check_adjoint_4d_online));
#endif
#endif
#endif
#ifdef HAVE_NFST
#undef X
#define X(name) NFST(name)
#ifndef _OPENMP
  nfst = CU_add_suite("nfst", 0, 0);
  CU_add_test(nfst, "nfst_1d_direct_file", X(check_1d_direct_file));
  CU_add_test(nfst, "nfst_1d_fast_file", X(check_1d_fast_file));
  CU_add_test(nfst, "nfst_adjoint_1d_direct_file", X(check_adjoint_1d_direct_file));
  CU_add_test(nfst, "nfst_adjoint_1d_fast_file", X(check_adjoint_1d_fast_file));
  CU_add_test(nfst, "nfst_1d_online", X(check_1d_online));
  CU_add_test(nfst, "nfst_adjoint_1d_online", X(check_adjoint_1d_online));

  CU_add_test(nfst, "nfst_2d_direct_file", X(check_2d_direct_file));
  CU_add_test(nfst, "nfst_2d_fast_file", X(check_2d_fast_file));
  CU_add_test(nfst, "nfst_adjoint_2d_direct_file", X(check_adjoint_2d_direct_file));
  CU_add_test(nfst, "nfst_adjoint_2d_fast_file", X(check_adjoint_2d_fast_file));
  CU_add_test(nfst, "nfst_2d_online", X(check_2d_online));
  CU_add_test(nfst, "nfst_adjoint_2d_online", X(check_adjoint_2d_online));

  CU_add_test(nfst, "nfst_3d_direct_file", X(check_3d_direct_file));
  CU_add_test(nfst, "nfst_3d_fast_file", X(check_3d_fast_file));
  CU_add_test(nfst, "nfst_adjoint_3d_direct_file", X(check_adjoint_3d_direct_file));
  CU_add_test(nfst, "nfst_adjoint_3d_fast_file", X(check_adjoint_3d_fast_file));
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  CU_add_test(nfst, "nfst_3d_online", X(check_3d_online));
  CU_add_test(nfst, "nfst_adjoint_3d_online", X(check_adjoint_3d_online));

  CU_add_test(nfst, "nfst_4d_online", X(check_4d_online));
  CU_add_test(nfst, "nfst_adjoint_4d_online", X(check_adjoint_4d_online));
#endif
#endif
#endif
  CU_automated_run_tests();
  //CU_basic_run_tests();
  {
    unsigned int ok = (CU_get_number_of_tests_failed() == 0);
    CU_cleanup_registry();
    return IF(ok, EXIT_SUCCESS, EXIT_FAILURE);
  }
}
