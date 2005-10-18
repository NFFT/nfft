AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/lowlevel/nfsft/obsolete
AM_LDFLAGS = -L$(top_srcdir)

TESTS = nfsft_test
  
## The program name
noinst_PROGRAMS = nfsft_test

## The source files.
nfsft_test_SOURCES = nfsft_test.c

nfsft_test_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
