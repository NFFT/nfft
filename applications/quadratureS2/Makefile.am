AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_PROGRAMS = quadratureS2

quadratureS2_SOURCES = quadratureS2.c doxygen.h
quadratureS2_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = quadratureS2.m writeTestcase.m writeWeights.m readTestcase.m lgwt.m plotGrid.m example.in example.out
