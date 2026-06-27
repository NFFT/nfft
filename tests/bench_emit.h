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

#ifndef NFFT_TESTS_BENCH_EMIT_H
#define NFFT_TESTS_BENCH_EMIT_H

/* Append one raw accuracy datapoint as an NDJSON record to the file named by
 * the NFFT_BENCH_OUT environment variable. A no-op (beyond a getenv) when that
 * variable is unset or empty, so ordinary `make check` runs are unaffected.
 *
 * Grouping and aggregation are NOT done here -- the converter
 * (tests/bench/ndjson_to_bmf.py) collapses N/M and builds the metric names.
 * `oracle` is "file" or "online". Values are emitted as long double with ample
 * precision and are valid JSON numbers. */
void bench_emit_accuracy(const char *module, const char *oracle, int d,
                         const int *N, int M, const char *init_name,
                         const char *trafo_name, long double accuracy,
                         long double bound, int ok);

#endif /* NFFT_TESTS_BENCH_EMIT_H */
