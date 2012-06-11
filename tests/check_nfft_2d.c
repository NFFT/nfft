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

#include <CUnit/CUnit.h>

#include "check_nfft.h"
#include "infft.h"

static const testcase_delegate_file_t nfft_1d_1_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_1.txt"};
static const testcase_delegate_file_t nfft_1d_1_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_10.txt"};
static const testcase_delegate_file_t nfft_1d_1_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_20.txt"};
static const testcase_delegate_file_t nfft_1d_1_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_50.txt"};
static const testcase_delegate_file_t nfft_1d_2_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_1.txt"};
static const testcase_delegate_file_t nfft_1d_2_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_10.txt"};
static const testcase_delegate_file_t nfft_1d_2_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_20.txt"};
static const testcase_delegate_file_t nfft_1d_2_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_50.txt"};
static const testcase_delegate_file_t nfft_1d_4_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_1.txt"};
static const testcase_delegate_file_t nfft_1d_4_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_10.txt"};
static const testcase_delegate_file_t nfft_1d_4_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_20.txt"};
static const testcase_delegate_file_t nfft_1d_4_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_50.txt"};
static const testcase_delegate_file_t nfft_1d_10_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_1.txt"};
static const testcase_delegate_file_t nfft_1d_10_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_10.txt"};
static const testcase_delegate_file_t nfft_1d_10_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_20.txt"};
static const testcase_delegate_file_t nfft_1d_10_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_50.txt"};
static const testcase_delegate_file_t nfft_1d_20_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_1.txt"};
static const testcase_delegate_file_t nfft_1d_20_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_10.txt"};
static const testcase_delegate_file_t nfft_1d_20_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_20.txt"};
static const testcase_delegate_file_t nfft_1d_20_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_50.txt"};
static const testcase_delegate_file_t nfft_1d_50_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_1.txt"};
static const testcase_delegate_file_t nfft_1d_50_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_10.txt"};
static const testcase_delegate_file_t nfft_1d_50_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_20.txt"};
static const testcase_delegate_file_t nfft_1d_50_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_50.txt"};

static const testcase_delegate_file_t *testcases[] =
{
  &nfft_1d_1_1, &nfft_1d_1_10, &nfft_1d_1_20, &nfft_1d_1_50,
  &nfft_1d_2_1, &nfft_1d_2_10, &nfft_1d_2_20, &nfft_1d_2_50,
  &nfft_1d_4_1, &nfft_1d_4_10, &nfft_1d_4_20, &nfft_1d_4_50,
  &nfft_1d_10_1, &nfft_1d_10_10, &nfft_1d_10_20, &nfft_1d_10_50,
  &nfft_1d_20_1, &nfft_1d_20_10, &nfft_1d_20_20, &nfft_1d_20_50,
  &nfft_1d_50_1, &nfft_1d_50_10, &nfft_1d_50_20, &nfft_1d_50_50,
};

/* Delegates. */
static const init_delegate_t* initializers[] =
{
  &init_2d,
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const trafo_delegate_t* trafos[] = {&trafo_direct, &trafo};

static void check_nfft_2d_file(void)
{
  X(check_many)(SIZE(testcases), SIZE(initializers), SIZE(trafos), testcases,
    initializers, trafos);
}

int main(void)
{
  CU_pSuite s;
  CU_initialize_registry();
  CU_set_output_filename("nfft");
  s = CU_add_suite("nfft_2d", 0, 0);
  CU_add_test(s, "nfft_2d_file", check_nfft_2d_file);
  CU_automated_run_tests();
  /*CU_basic_run_tests();*/
  CU_cleanup_registry();
  return 0;
}
