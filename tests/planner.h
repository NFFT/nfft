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

#ifndef PLANNER_TEST_H
#define PLANNER_TEST_H

#include "infft.h"

void Y(check_planner_md5_rfc_vectors)(void);
void Y(check_planner_md5_feeders)(void);
void Y(check_planner_printer)(void);
void Y(check_planner_scanner)(void);
void Y(check_planner_registry)(void);
void Y(check_planner_hashtable)(void);
void Y(check_planner_subsumption)(void);
void Y(check_planner_forget)(void);
void Y(check_planner_wisdom_roundtrip)(void);
void Y(check_planner_wisdom_rejects)(void);
void Y(check_planner_tensor_basic)(void);
void Y(check_planner_tensor_canonical)(void);
void Y(check_planner_trinity_mkplan)(void);
void Y(check_planner_trinity_wisdom_memo)(void);
void Y(check_planner_trinity_print)(void);
void Y(check_planner_measure)(void);
void Y(check_planner_candidates)(void);
void Y(check_planner_bless)(void);
void Y(check_planner_timelimit_default_and_set)(void);
void Y(check_planner_clock_now_monotonic)(void);

#endif
