AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = radon polar_fft_test mpolar_fft_test linogram_fft_test

## The source files.
radon_SOURCES = radon.c
polar_fft_test_SOURCES = polar_fft_test.c
mpolar_fft_test_SOURCES = mpolar_fft_test.c
linogram_fft_test_SOURCES = linogram_fft_test.c

radon_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
polar_fft_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
mpolar_fft_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
linogram_fft_test_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3
