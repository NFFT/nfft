AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_PROGRAMS = fastsumS2

fastsumS2_SOURCES = fastsumS2.c doxygen.h
fastsumS2_LDADD = $(top_builddir)/libnfft3.la -lm

EXTRA_DIST = fastsumS2.m writeTestcase.m readTestcase.m README example.in example.out
