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

#include <stdio.h>
#include <stdlib.h>

#include "bench_emit.h"

void bench_emit_accuracy(const char *module, const char *oracle, int d,
                         const int *N, int M, const char *init_name,
                         const char *trafo_name, long double accuracy,
                         long double bound, int ok)
{
  const char *base = getenv("NFFT_BENCH_OUT");
  const char *path;
  FILE *fp;
  int j;
#ifdef _OPENMP
  static char threads_path[4096];
  const int openmp = 1;
#else
  const int openmp = 0;
#endif

  if (base == NULL || base[0] == '\0')
    return;

  /* The serial (checkall) and OpenMP (checkall_threads) binaries both run under
     `make check`; route the OpenMP build to a distinct "<base>.threads" file so
     their records never interleave. The compile-time _OPENMP macro is the only
     thing distinguishing the two binaries (they share this source). */
#ifdef _OPENMP
  {
    int n = snprintf(threads_path, sizeof(threads_path), "%s.threads", base);
    if (n < 0 || (size_t)n >= sizeof(threads_path))
      return;
    path = threads_path;
  }
#else
  path = base;
#endif

  fp = fopen(path, "a");
  if (fp == NULL)
    return;

  fprintf(fp, "{\"module\": \"%s\", \"oracle\": \"%s\", \"openmp\": %d, "
              "\"dim\": %d, \"N\": [", module, oracle, openmp, d);
  for (j = 0; j < d; j++)
    fprintf(fp, "%s%d", j ? "," : "", N[j]);
  /* %.20Le yields a valid JSON number (no leading '+'/spaces; err,bound >= 0). */
  fprintf(fp, "], \"M\": %d, \"init\": \"%s\", \"trafo\": \"%s\", "
              "\"accuracy\": %.20Le, \"bound\": %.20Le, \"ok\": %d}\n",
          M, init_name, trafo_name, accuracy, bound, ok);

  fclose(fp);
}
