AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

TESTS = test

## The program name
nodist_PROGRAMS = test
nodistdir = .

## The source files.
test_SOURCES = test.c

test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
