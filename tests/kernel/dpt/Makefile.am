AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/kernel/dpt -I$(top_srcdir)/lowlevel/nfsft/obsolete
AM_LDFLAGS = -L$(top_srcdir)

TESTS = test

## The program name
noinst_PROGRAMS = test

## The source files.
test_SOURCES = test.c

test_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3 -lcunit
