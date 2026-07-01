#ifndef NFFT_TESTS_ACCURACY_LOG_H
#define NFFT_TESTS_ACCURACY_LOG_H

/**
 * Append one raw accuracy datapoint as an NDJSON record to the file named by
 * the NFFT_BENCH_OUT environment variable. A no-op if variable is unset or empty. */
void accuracy_log_append(const char *module, const char *oracle, int d,
                         const int *N, int M, const char *init_name,
                         const char *trafo_name, long double accuracy,
                         long double bound, int ok);

#endif /* NFFT_TESTS_ACCURACY_LOG_H */
