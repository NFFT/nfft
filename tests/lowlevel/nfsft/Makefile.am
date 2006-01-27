AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_builddir)/.libs

TESTS = test

## The program name
noinst_PROGRAMS = test

## The source files.
test_SOURCES = test.c

test_LDADD = -lnfft3 -lm -lfftw3 -lcunit