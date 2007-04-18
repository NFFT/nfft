#include "nfft3.h"
#include <stdio.h>
#include <stdlib.h>

nfft_malloc_type_function nfft_malloc_hook = 0;
nfft_free_type_function nfft_free_hook = 0;
nfft_die_type_function nfft_die_hook = 0;

void *nfft_malloc(size_t n)
{
  void *p;

  if (nfft_malloc_hook)
  {
    return nfft_malloc_hook(n);
  }

  if (n == 0)
  {
    n = 1;
  }

  p = fftw_malloc(n);

  if (!p)
  {
    nfft_die("nfft_malloc: out of memory\n");
  }

  return p;
}

void nfft_free(void *p)
{
  if (p) 
  {
	if (nfft_free_hook) 
	{
	  nfft_free_hook(p);
	  return;
	}

	fftw_free(p);
  }
}

void nfft_die(const char *s)
{
  if (nfft_die_hook)
  {
    nfft_die_hook(s);

    fflush(stdout);
    fprintf(stderr, "nfft: %s", s);
    exit(EXIT_FAILURE);
  }
}
