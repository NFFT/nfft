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
#include "ngomp.h"

int main(void)
{
  CU_pSuite ngomp;
  CU_initialize_registry();
  CU_set_output_filename("CUnitAutomated_ngomp");

  ngomp = CU_add_suite("ngomp", 0, 0);
  CU_add_test(ngomp, "empty_roster", Y(check_ngomp_empty_roster));
  CU_add_test(ngomp, "plans_serially", Y(check_ngomp_plans_serially));

  CU_automated_run_tests();
  // CU_basic_run_tests();
  {
    unsigned int ok = (CU_get_number_of_tests_failed() == 0);
    CU_cleanup_registry();
    return IF(ok, EXIT_SUCCESS, EXIT_FAILURE);
  }
}
