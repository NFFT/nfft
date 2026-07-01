#include <stdio.h>
#include <stdlib.h>

#include "accuracy_log.h"

void accuracy_log_append(const char *module, const char *oracle, int d,
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

// Route to serial/threaded file, respectively. 
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
  // %.20Le yields a valid JSON number.
  fprintf(fp, "], \"M\": %d, \"init\": \"%s\", \"trafo\": \"%s\", "
              "\"accuracy\": %.20Le, \"bound\": %.20Le, \"ok\": %d}\n",
          M, init_name, trafo_name, accuracy, bound, ok);

  fclose(fp);
}
