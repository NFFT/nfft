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
#include "plan.h"
#include "window.h"
#include "deconv.h"
#include "conv.h"
#include "fast_native.h"
#include "nfft_ng.h"
#include "tune_ng.h"

int main(void) {
  CU_pSuite util, planner_suite, plan_suite, window_suite, deconv_suite,
      conv_suite, fast_native_suite, nfft_ng_suite, tune_suite;
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

  plan_suite = CU_add_suite("plan", 0, 0);
  CU_add_test(plan_suite, "problem", Y(check_plan_problem));
  CU_add_test(plan_suite, "problem_owns_data", Y(check_plan_problem_owns_data));
  CU_add_test(plan_suite, "key_is_data_blind", Y(check_plan_key_is_data_blind));
  CU_add_test(plan_suite, "elides_unit_axes_geometry", Y(check_plan_elides_unit_axes_geometry));
  CU_add_test(plan_suite, "variant_key", Y(check_plan_variant_key));
  CU_add_test(plan_suite, "guard_declines", Y(check_plan_guard_declines));
  CU_add_test(plan_suite, "solvers", Y(check_plan_solvers));
  CU_add_test(plan_suite, "ndft_dispatch", Y(check_plan_ndft_dispatch));
  CU_add_test(plan_suite, "ndft_multivariate_dispatch", Y(check_plan_ndft_multivariate_dispatch));
  CU_add_test(plan_suite, "ndft_accuracy", Y(check_plan_ndft_accuracy));
  CU_add_test(plan_suite, "correct", Y(check_plan_correct));
  CU_add_test(plan_suite, "wisdom_memo", Y(check_plan_wisdom_memo));
  CU_add_test(plan_suite, "measured", Y(check_plan_measured));
  CU_add_test(plan_suite, "destructive_default", Y(check_plan_destructive_default));
  CU_add_test(plan_suite, "execute_on", Y(check_plan_execute_on));
  CU_add_test(plan_suite, "measured_wisdom", Y(check_plan_measured_wisdom));
  CU_add_test(plan_suite, "measured_prunes_by_estimate",
              Y(check_plan_measured_prunes_by_estimate));
  CU_add_test(plan_suite, "timelimit_tight_degrades_to_estimate", Y(check_plan_timelimit_tight_degrades_to_estimate));
  CU_add_test(plan_suite, "timelimit_unset_measures_and_blesses", Y(check_plan_timelimit_unset_measures_and_blesses));
  CU_add_test(plan_suite, "set_timelimit_roundtrip", Y(check_plan_set_timelimit_roundtrip));
  CU_add_test(plan_suite, "public_api", Y(check_plan_public_api));
  CU_add_test(plan_suite, "print_includes_registrar_names", Y(check_plan_print_includes_registrar_names));
  CU_add_test(plan_suite, "forward_only_race", Y(check_plan_forward_only_race));
  CU_add_test(plan_suite, "apply_adjoint", Y(check_plan_apply_adjoint));
  CU_add_test(plan_suite, "x_copied_not_aliased", Y(check_plan_x_copied_not_aliased));
  CU_add_test(plan_suite, "per_plan_core", Y(check_plan_per_plan_core));
  CU_add_test(plan_suite, "core_owns_no_data_arrays", Y(check_plan_core_owns_no_data_arrays));
  CU_add_test(plan_suite, "core_elision", Y(check_plan_core_elision));
  CU_add_test(plan_suite, "variant_guru", Y(check_plan_variant_guru));
  CU_add_test(plan_suite, "odd_n", Y(check_plan_odd_n));
  CU_add_test(plan_suite, "type_ii_1d", Y(check_plan_type_ii_1d));
  CU_add_test(plan_suite, "type_ii_nd", Y(check_plan_type_ii_nd));
  CU_add_test(plan_suite, "awake_states", Y(check_plan_awake_states));
  CU_add_test(plan_suite, "user_x_pristine", Y(check_plan_user_x_pristine));
  CU_add_test(plan_suite, "awake_zero_restore", Y(check_plan_awake_zero_restore));
  CU_add_test(plan_suite, "restore_guard_fires", Y(check_plan_restore_guard_fires));
  CU_add_test(plan_suite, "in_race_guard_passes", Y(check_plan_in_race_guard_passes));
  CU_add_test(plan_suite, "awake_zero_internal", Y(check_plan_awake_zero_internal));
  CU_add_test(plan_suite, "elides_unit_axes", Y(check_plan_elides_unit_axes));
  CU_add_test(plan_suite, "rank0_solver", Y(check_plan_rank0_solver));
  CU_add_test(plan_suite, "unit_axis_correct", Y(check_plan_unit_axis_correct));
  CU_add_test(plan_suite, "newarray_native_fast", Y(check_plan_newarray_native_fast));
  CU_add_test(plan_suite, "unit_axis_execute_on", Y(check_plan_unit_axis_execute_on));
  CU_add_test(plan_suite, "guru rejects NULL args", Y(check_plan_guru_rejects_null_args));
  CU_add_test(plan_suite, "guru rejects bad geometry", Y(check_plan_guru_rejects_bad_geometry));

  window_suite = CU_add_suite("window", 0, 0);
  CU_add_test(window_suite, "id", Y(check_window_id));
  CU_add_test(window_suite, "vtable", Y(check_window_vtable));
  CU_add_test(window_suite, "normalized", Y(check_window_normalized));
  CU_add_test(window_suite, "all", Y(check_window_all));
  CU_add_test(window_suite, "key", Y(check_window_key));
  CU_add_test(window_suite, "bessel_log_scaled", Y(check_window_bessel_log_scaled));
  CU_add_test(window_suite, "cancellation", Y(check_window_cancellation));
  CU_add_test(window_suite, "apply", Y(check_window_apply));

  deconv_suite = CU_add_suite("deconv", 0, 0);
  CU_add_test(deconv_suite, "problem", Y(check_deconv_problem));
  CU_add_test(deconv_suite, "solver", Y(check_deconv_solver));
  CU_add_test(deconv_suite, "1d", Y(check_deconv_1d));
  CU_add_test(deconv_suite, "2d", Y(check_deconv_2d));
  CU_add_test(deconv_suite, "3d", Y(check_deconv_3d));
  CU_add_test(deconv_suite, "nd", Y(check_deconv_nd));

  conv_suite = CU_add_suite("conv", 0, 0);
  CU_add_test(conv_suite, "problem", Y(check_conv_problem));
  CU_add_test(conv_suite, "solver", Y(check_conv_solver));
  CU_add_test(conv_suite, "1d", Y(check_conv_1d));
  CU_add_test(conv_suite, "2d", Y(check_conv_2d));
  CU_add_test(conv_suite, "3d", Y(check_conv_3d));
  CU_add_test(conv_suite, "nd", Y(check_conv_nd));

  fast_native_suite = CU_add_suite("fast_native", 0, 0);
  CU_add_test(fast_native_suite, "tree", Y(check_fast_native_tree));
  CU_add_test(fast_native_suite, "declines_window", Y(check_fast_native_declines_window));
  CU_add_test(fast_native_suite, "flag_selective", Y(check_fast_native_flag_selective));
  CU_add_test(fast_native_suite, "2d", Y(check_fast_native_2d));
  CU_add_test(fast_native_suite, "2d_adjoint", Y(check_fast_native_2d_adjoint));
  CU_add_test(fast_native_suite, "3d", Y(check_fast_native_3d));
  CU_add_test(fast_native_suite, "3d_adjoint", Y(check_fast_native_3d_adjoint));
  CU_add_test(fast_native_suite, "4d", Y(check_fast_native_4d));
  CU_add_test(fast_native_suite, "4d_adjoint", Y(check_fast_native_4d_adjoint));
  CU_add_test(fast_native_suite, "window_select", Y(check_fast_native_window_select));
  CU_add_test(fast_native_suite, "window_accuracy", Y(check_fast_native_window_accuracy));

  tune_suite = CU_add_suite("tune", 0, 0);
  CU_add_test(tune_suite, "meets_goal", Y(check_tune_meets_goal));
  CU_add_test(tune_suite, "unreachable", Y(check_tune_unreachable));
  CU_add_test(tune_suite, "geometries", Y(check_tune_geometries));
  CU_add_test(tune_suite, "bad_args", Y(check_tune_bad_args));
  CU_add_test(tune_suite, "sigma_agrees", Y(check_tune_sigma_agrees));
  CU_add_test(tune_suite, "sigma_limits", Y(check_tune_sigma_limits));
  CU_add_test(tune_suite, "plan", Y(check_tune_plan));
  CU_add_test(tune_suite, "plan_capped", Y(check_tune_plan_capped));

  nfft_ng_suite = CU_add_suite("nfft_ng", 0, 0);
  CU_add_test(nfft_ng_suite, "nfft_ng_1d_file", Y(check_nfft_ng_1d_file));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_1d_file", Y(check_nfft_ng_adjoint_1d_file));
  CU_add_test(nfft_ng_suite, "nfft_ng_1d_online", Y(check_nfft_ng_1d_online));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_1d_online", Y(check_nfft_ng_adjoint_1d_online));

  CU_add_test(nfft_ng_suite, "nfft_ng_2d_file", Y(check_nfft_ng_2d_file));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_2d_file", Y(check_nfft_ng_adjoint_2d_file));
  CU_add_test(nfft_ng_suite, "nfft_ng_2d_online", Y(check_nfft_ng_2d_online));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_2d_online", Y(check_nfft_ng_adjoint_2d_online));

  CU_add_test(nfft_ng_suite, "nfft_ng_3d_file", Y(check_nfft_ng_3d_file));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_3d_file", Y(check_nfft_ng_adjoint_3d_file));

  CU_add_test(nfft_ng_suite, "nfft_ng_4d_file", Y(check_nfft_ng_4d_file));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_4d_file", Y(check_nfft_ng_adjoint_4d_file));
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  CU_add_test(nfft_ng_suite, "nfft_ng_3d_online", Y(check_nfft_ng_3d_online));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_3d_online", Y(check_nfft_ng_adjoint_3d_online));

  CU_add_test(nfft_ng_suite, "nfft_ng_4d_online", Y(check_nfft_ng_4d_online));
  CU_add_test(nfft_ng_suite, "nfft_ng_adjoint_4d_online", Y(check_nfft_ng_adjoint_4d_online));
#endif

  CU_add_test(nfft_ng_suite, "nfft_ng_fast_variants_1d", Y(check_nfft_ng_fast_variants_1d));
  CU_add_test(nfft_ng_suite, "nfft_ng_fast_variants_2d", Y(check_nfft_ng_fast_variants_2d));
  CU_add_test(nfft_ng_suite, "nfft_ng_fast_variants_3d", Y(check_nfft_ng_fast_variants_3d));
  CU_add_test(nfft_ng_suite, "nfft_ng_fast_variants_4d", Y(check_nfft_ng_fast_variants_4d));

  CU_add_test(nfft_ng_suite, "nfft_ng_fast_unit_axes_1d", Y(check_nfft_ng_fast_unit_axes_1d));
  CU_add_test(nfft_ng_suite, "nfft_ng_fast_unit_axes_2d", Y(check_nfft_ng_fast_unit_axes_2d));
  CU_add_test(nfft_ng_suite, "nfft_ng_fast_unit_axes_3d", Y(check_nfft_ng_fast_unit_axes_3d));
  CU_add_test(nfft_ng_suite, "nfft_ng_fast_unit_axes_4d", Y(check_nfft_ng_fast_unit_axes_4d));

  CU_automated_run_tests();
  // CU_basic_run_tests();
  {
    unsigned int ok = (CU_get_number_of_tests_failed() == 0);
    CU_cleanup_registry();
    return IF(ok, EXIT_SUCCESS, EXIT_FAILURE);
  }
}
