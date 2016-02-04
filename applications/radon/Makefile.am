AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_PROGRAMS = radon inverse_radon

radon_SOURCES = radon.c
radon_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

inverse_radon_SOURCES = inverse_radon.c
inverse_radon_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = phantom.m radon.m ridgelet.m README
