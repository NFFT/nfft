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
  const char *path = getenv("NFFT_BENCH_OUT");
  FILE *fp;
  int j;

  if (path == NULL || path[0] == '\0')
    return;

  fp = fopen(path, "a");
  if (fp == NULL)
    return;

  fprintf(fp, "{\"module\": \"%s\", \"oracle\": \"%s\", \"dim\": %d, \"N\": [",
          module, oracle, d);
  for (j = 0; j < d; j++)
    fprintf(fp, "%s%d", j ? "," : "", N[j]);
  /* %.20Le yields a valid JSON number (no leading '+'/spaces; err,bound >= 0). */
  fprintf(fp, "], \"M\": %d, \"init\": \"%s\", \"trafo\": \"%s\", "
              "\"accuracy\": %.20Le, \"bound\": %.20Le, \"ok\": %d}\n",
          M, init_name, trafo_name, accuracy, bound, ok);

  fclose(fp);
}
