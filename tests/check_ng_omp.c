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

/* Standard headers. */
#include <CUnit/CUnit.h>
#include <CUnit/TestRun.h>
#include <CUnit/Automated.h>

#include "infft.h"
#include "ng_omp.h"

int main(void)
{
  CU_pSuite ng_omp;
  CU_initialize_registry();
  CU_set_output_filename("CUnitAutomated_ng_omp");

  ng_omp = CU_add_suite("ng_omp", 0, 0);
  CU_add_test(ng_omp, "roster", Y(check_ng_omp_roster));
  CU_add_test(ng_omp, "deconv_1d", Y(check_ng_omp_deconv_1d));
  CU_add_test(ng_omp, "deconv_2d", Y(check_ng_omp_deconv_2d));
  CU_add_test(ng_omp, "deconv_3d", Y(check_ng_omp_deconv_3d));
  CU_add_test(ng_omp, "deconv_4d", Y(check_ng_omp_deconv_4d));
  CU_add_test(ng_omp, "deconv_variants", Y(check_ng_omp_deconv_variants));
  CU_add_test(ng_omp, "thread_count_invariance",
              Y(check_ng_omp_thread_count_invariance));
  CU_add_test(ng_omp, "declines_serial_budget",
              Y(check_ng_omp_declines_serial_budget));
  CU_add_test(ng_omp, "conv_forward", Y(check_ng_omp_conv_forward));
  CU_add_test(ng_omp, "conv_adjoint", Y(check_ng_omp_conv_adjoint));
  CU_add_test(ng_omp, "conv_narrow_grid", Y(check_ng_omp_conv_narrow_grid));
  CU_add_test(ng_omp, "conv_clustered", Y(check_ng_omp_conv_clustered));
  CU_add_test(ng_omp, "conv_wrap", Y(check_ng_omp_conv_wrap));
  CU_add_test(ng_omp, "race_memory", Y(check_ng_omp_race_memory));
  CU_add_test(ng_omp, "conv_sorted_declines",
              Y(check_ng_omp_conv_sorted_declines));
  CU_add_test(ng_omp, "child_race", Y(check_ng_omp_child_race));
  CU_add_test(ng_omp, "nfft_all_levels", Y(check_ng_omp_nfft_all_levels));
  CU_add_test(ng_omp, "no_fast_native_refuses",
              Y(check_ng_omp_no_fast_native_refuses));

  CU_automated_run_tests();
  // CU_basic_run_tests();
  {
    unsigned int ok = (CU_get_number_of_tests_failed() == 0);
    CU_cleanup_registry();
    return IF(ok, EXIT_SUCCESS, EXIT_FAILURE);
  }
}
