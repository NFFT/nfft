AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_PROGRAMS = simple_test accuracy

simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

accuracy_SOURCES = accuracy.c
accuracy_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = accuracy.m README
