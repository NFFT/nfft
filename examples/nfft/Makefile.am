# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_PROGRAMS = simple_test ndft_fast taylor_nfft flags nfft_times

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la

ndft_fast_SOURCES = ndft_fast.c
ndft_fast_LDADD = $(top_builddir)/libnfft3.la

taylor_nfft_SOURCES = taylor_nfft.c
taylor_nfft_LDADD = $(top_builddir)/libnfft3.la

flags_SOURCES = flags.c
flags_LDADD = $(top_builddir)/libnfft3.la

nfft_times_SOURCES = nfft_times.c
nfft_times_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = flags.m ndft_fast.m taylor_nfft.m README
