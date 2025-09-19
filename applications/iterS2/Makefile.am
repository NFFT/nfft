AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/kernel/nfsft

noinst_PROGRAMS = iterS2

iterS2_SOURCES = iterS2.c doxygen.h
iterS2_LDFLAGS = @fftw3_LDFLAGS@
iterS2_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@

EXTRA_DIST = writeTestcase.m writeImageTestcase.m README
#iterS2.m readTestcase.m
