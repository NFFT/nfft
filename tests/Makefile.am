# $Id: Makefile.am 3456 2010-03-17 19:25:47Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_PROGRAMS = check_bessel_i0 check_nfft_1d

check_bessel_i0_SOURCES = check_bessel_i0.c
check_bessel_i0_LDADD = $(top_builddir)/libnfft3.la

check_nfft_1d_SOURCES = check_nfft_1d.c
check_nfft_1d_LDADD = $(top_builddir)/libnfft3.la
