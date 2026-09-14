/* Live-byte accounting through Y(malloc_hook) and Y(free_hook), for tests
 * only. Every block freed while the hooks are installed must have been
 * allocated while they were installed: create and destroy everything the
 * measurement touches between alloc_count_begin() and alloc_count_end(). */
#ifndef NFFT_TESTS_ALLOC_COUNT_H
#define NFFT_TESTS_ALLOC_COUNT_H

#include <string.h>

#define ALLOC_COUNT_PAD 64 /* keeps FFTW's SIMD alignment for the user block */

static size_t alloc_count_live, alloc_count_peak;

static void *alloc_count_malloc(size_t n)
{
  unsigned char *p = (unsigned char *)FFTW(malloc)(n + ALLOC_COUNT_PAD);
  memcpy(p, &n, sizeof(n));
  alloc_count_live += n;
  if (alloc_count_live > alloc_count_peak)
    alloc_count_peak = alloc_count_live;
  return p + ALLOC_COUNT_PAD;
}

static void alloc_count_free(void *q)
{
  unsigned char *p = (unsigned char *)q - ALLOC_COUNT_PAD;
  size_t n;
  memcpy(&n, p, sizeof(n));
  alloc_count_live -= n;
  FFTW(free)(p);
}

static void alloc_count_begin(void)
{
  alloc_count_live = alloc_count_peak = 0;
  Y(malloc_hook) = alloc_count_malloc;
  Y(free_hook) = alloc_count_free;
}

static void alloc_count_end(void)
{
  Y(malloc_hook) = 0;
  Y(free_hook) = 0;
}

#endif
