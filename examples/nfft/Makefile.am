AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = simple_test ndft_fast taylor_nfft flags nfft_times

## The source files.
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

ndft_fast_SOURCES = ndft_fast.c
ndft_fast_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

taylor_nfft_SOURCES = taylor_nfft.c
taylor_nfft_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

flags_SOURCES = flags.c
flags_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

nfft_times_SOURCES = nfft_times.c
nfft_times_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3