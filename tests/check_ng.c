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
#include "planner.h"
#include "nplan.h"
#include "nfast.h"

int main(void) {
  CU_pSuite util, planner_suite, nplan_suite, nplan_data_suite, nfast_suite;
  CU_initialize_registry();
#ifdef _OPENMP
  CU_set_output_filename("CUnitAutomated_ng_threads");
#else
  CU_set_output_filename("CUnitAutomated_ng");
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

  planner_suite = CU_add_suite("planner", 0, 0);
  CU_add_test(planner_suite, "md5_rfc_vectors", Y(check_planner_md5_rfc_vectors));
  CU_add_test(planner_suite, "md5_feeders", Y(check_planner_md5_feeders));
  CU_add_test(planner_suite, "printer", Y(check_planner_printer));
  CU_add_test(planner_suite, "scanner", Y(check_planner_scanner));
  CU_add_test(planner_suite, "registry", Y(check_planner_registry));
  CU_add_test(planner_suite, "hashtable", Y(check_planner_hashtable));
  CU_add_test(planner_suite, "subsumption", Y(check_planner_subsumption));
  CU_add_test(planner_suite, "forget", Y(check_planner_forget));
  CU_add_test(planner_suite, "wisdom_roundtrip", Y(check_planner_wisdom_roundtrip));
  CU_add_test(planner_suite, "wisdom_rejects", Y(check_planner_wisdom_rejects));
  CU_add_test(planner_suite, "tensor_basic", Y(check_planner_tensor_basic));
  CU_add_test(planner_suite, "tensor_canonical", Y(check_planner_tensor_canonical));
  CU_add_test(planner_suite, "trinity_mkplan", Y(check_planner_trinity_mkplan));
  CU_add_test(planner_suite, "trinity_wisdom_memo", Y(check_planner_trinity_wisdom_memo));
  CU_add_test(planner_suite, "trinity_print", Y(check_planner_trinity_print));
  CU_add_test(planner_suite, "measure", Y(check_planner_measure));
  CU_add_test(planner_suite, "candidates", Y(check_planner_candidates));
  CU_add_test(planner_suite, "bless", Y(check_planner_bless));
  CU_add_test(planner_suite, "timelimit_default_and_set", Y(check_planner_timelimit_default_and_set));
  CU_add_test(planner_suite, "clock_now_monotonic", Y(check_planner_clock_now_monotonic));

  nplan_suite = CU_add_suite("nplan", 0, 0);
  CU_add_test(nplan_suite, "problem", Y(check_nplan_problem));
  CU_add_test(nplan_suite, "problem_owns_data", Y(check_nplan_problem_owns_data));
  CU_add_test(nplan_suite, "key_is_data_blind", Y(check_nplan_key_is_data_blind));
  CU_add_test(nplan_suite, "elides_unit_axes_geometry", Y(check_nplan_elides_unit_axes_geometry));
  CU_add_test(nplan_suite, "variant_key", Y(check_nplan_variant_key));
  CU_add_test(nplan_suite, "guard_declines", Y(check_nplan_guard_declines));
  CU_add_test(nplan_suite, "solvers", Y(check_nplan_solvers));
  CU_add_test(nplan_suite, "ndft_dispatch", Y(check_nplan_ndft_dispatch));
  CU_add_test(nplan_suite, "ndft_multivariate_dispatch", Y(check_nplan_ndft_multivariate_dispatch));
  CU_add_test(nplan_suite, "ndft_accuracy", Y(check_nplan_ndft_accuracy));
  CU_add_test(nplan_suite, "correct", Y(check_nplan_correct));
  CU_add_test(nplan_suite, "wisdom_memo", Y(check_nplan_wisdom_memo));
  CU_add_test(nplan_suite, "measured", Y(check_nplan_measured));
  CU_add_test(nplan_suite, "destructive_default", Y(check_nplan_destructive_default));
  CU_add_test(nplan_suite, "execute_on", Y(check_nplan_execute_on));
  CU_add_test(nplan_suite, "measured_wisdom", Y(check_nplan_measured_wisdom));
  CU_add_test(nplan_suite, "measured_prunes_by_estimate",
              Y(check_nplan_measured_prunes_by_estimate));
  CU_add_test(nplan_suite, "timelimit_tight_degrades_to_estimate", Y(check_nplan_timelimit_tight_degrades_to_estimate));
  CU_add_test(nplan_suite, "timelimit_unset_measures_and_blesses", Y(check_nplan_timelimit_unset_measures_and_blesses));
  CU_add_test(nplan_suite, "set_timelimit_roundtrip", Y(check_nplan_set_timelimit_roundtrip));
  CU_add_test(nplan_suite, "public_api", Y(check_nplan_public_api));
  CU_add_test(nplan_suite, "print_includes_registrar_names", Y(check_nplan_print_includes_registrar_names));
  CU_add_test(nplan_suite, "forward_only_race", Y(check_nplan_forward_only_race));
  CU_add_test(nplan_suite, "apply_adjoint", Y(check_nplan_apply_adjoint));
  CU_add_test(nplan_suite, "x_copied_not_aliased", Y(check_nplan_x_copied_not_aliased));
  CU_add_test(nplan_suite, "per_plan_core", Y(check_nplan_per_plan_core));
  CU_add_test(nplan_suite, "core_owns_no_data_arrays", Y(check_nplan_core_owns_no_data_arrays));
  CU_add_test(nplan_suite, "core_elision", Y(check_nplan_core_elision));
  CU_add_test(nplan_suite, "variant_guru", Y(check_nplan_variant_guru));
  CU_add_test(nplan_suite, "odd_n", Y(check_nplan_odd_n));
  CU_add_test(nplan_suite, "type_ii_1d", Y(check_nplan_type_ii_1d));
  CU_add_test(nplan_suite, "type_ii_nd", Y(check_nplan_type_ii_nd));
  CU_add_test(nplan_suite, "awake_states", Y(check_nplan_awake_states));
  CU_add_test(nplan_suite, "user_x_pristine", Y(check_nplan_user_x_pristine));
  CU_add_test(nplan_suite, "awake_zero_restore", Y(check_nplan_awake_zero_restore));
  CU_add_test(nplan_suite, "restore_guard_fires", Y(check_nplan_restore_guard_fires));
  CU_add_test(nplan_suite, "in_race_guard_passes", Y(check_nplan_in_race_guard_passes));
  CU_add_test(nplan_suite, "awake_zero_internal", Y(check_nplan_awake_zero_internal));
  CU_add_test(nplan_suite, "elides_unit_axes", Y(check_nplan_elides_unit_axes));
  CU_add_test(nplan_suite, "rank0_solver", Y(check_nplan_rank0_solver));
  CU_add_test(nplan_suite, "unit_axis_correct", Y(check_nplan_unit_axis_correct));
  CU_add_test(nplan_suite, "newarray_native_fast", Y(check_nplan_newarray_native_fast));
  CU_add_test(nplan_suite, "unit_axis_execute_on", Y(check_nplan_unit_axis_execute_on));
  CU_add_test(nplan_suite, "guru rejects NULL args", Y(check_nplan_guru_rejects_null_args));
  CU_add_test(nplan_suite, "guru rejects bad geometry", Y(check_nplan_guru_rejects_bad_geometry));

  nplan_data_suite = CU_add_suite("nplan_data", 0, 0);
  CU_add_test(nplan_data_suite, "native_reference", Y(check_nplan_data));

  nfast_suite = CU_add_suite("nfast", 0, 0);
  CU_add_test(nfast_suite, "window_id", Y(check_nfast_window_id));
  CU_add_test(nfast_suite, "deconv_problem", Y(check_nfast_deconv_problem));
  CU_add_test(nfast_suite, "conv_problem", Y(check_nfast_conv_problem));
  CU_add_test(nfast_suite, "window_vtable", Y(check_nfast_window_vtable));
  CU_add_test(nfast_suite, "window_normalized", Y(check_nfast_window_normalized));
  CU_add_test(nfast_suite, "window_all", Y(check_nfast_window_all));
  CU_add_test(nfast_suite, "bessel_log_scaled", Y(check_nfast_bessel_log_scaled));
  CU_add_test(nfast_suite, "window_key", Y(check_nfast_window_key));
  CU_add_test(nfast_suite, "deconv_solver", Y(check_nfast_deconv_solver));
  CU_add_test(nfast_suite, "deconv_1d_general", Y(check_nfast_deconv_1d_general));
  CU_add_test(nfast_suite, "deconv_2d_general", Y(check_nfast_deconv_2d_general));
  CU_add_test(nfast_suite, "deconv_3d_general", Y(check_nfast_deconv_3d_general));
  CU_add_test(nfast_suite, "deconv_nd_general", Y(check_nfast_deconv_nd_general));
  CU_add_test(nfast_suite, "conv_solver", Y(check_nfast_conv_solver));
  CU_add_test(nfast_suite, "native_fast_accuracy", Y(check_nfast_native_fast_accuracy));
  CU_add_test(nfast_suite, "native_tree", Y(check_nfast_native_tree));
  CU_add_test(nfast_suite, "native_declines_window", Y(check_nfast_native_declines_window));
  CU_add_test(nfast_suite, "flag_selective", Y(check_nfast_flag_selective));
  CU_add_test(nfast_suite, "native_fast_2d", Y(check_nfast_native_fast_2d));
  CU_add_test(nfast_suite, "native_fast_2d_adjoint", Y(check_nfast_native_fast_2d_adjoint));
  CU_add_test(nfast_suite, "native_fast_3d", Y(check_nfast_native_fast_3d));
  CU_add_test(nfast_suite, "native_fast_3d_adjoint", Y(check_nfast_native_fast_3d_adjoint));
  CU_add_test(nfast_suite, "native_fast_4d", Y(check_nfast_native_fast_4d));
  CU_add_test(nfast_suite, "native_fast_4d_adjoint", Y(check_nfast_native_fast_4d_adjoint));
  CU_add_test(nfast_suite, "float_accuracy", Y(check_nfast_float_accuracy));
  CU_add_test(nfast_suite, "window_cancellation", Y(check_nfast_window_cancellation));
  CU_add_test(nfast_suite, "window_apply", Y(check_nfast_window_apply));
  CU_add_test(nfast_suite, "native_window_select", Y(check_nfast_native_window_select));
  CU_add_test(nfast_suite, "window_accuracy", Y(check_nfast_window_accuracy));

  CU_automated_run_tests();
  // CU_basic_run_tests();
  {
    unsigned int ok = (CU_get_number_of_tests_failed() == 0);
    CU_cleanup_registry();
    return IF(ok, EXIT_SUCCESS, EXIT_FAILURE);
  }
}
