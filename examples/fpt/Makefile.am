AM_CPPFLAGS = -I$(top_srcdir)/include @NFFT_PROGRAM_CPPFLAGS@

# programs
noinst_PROGRAMS = simple_test

# sources
simple_test_SOURCES = simple_test.c
simple_test_LDFLAGS = @fftw3_LDFLAGS@
simple_test_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@
