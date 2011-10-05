# $Id: Makefile.am 3456 2010-03-17 19:25:47Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_NFSFT
  CHECK_NFSFT=check_nfsft
  check_nfsft_SOURCES = check_nfsft.c
  check_nfsft_LDADD = $(top_builddir)/libnfft3.la
else
  CHECK_NFSFT=
endif

check_LTLIBRARIES = libcheck_nfft.la

libcheck_nfft_la_SOURCES = check_nfft.c check_nfft.h

check_PROGRAMS = check_nfft_1d
#check_bessel_i0 check_nfft_1d check_nfft_2d $(CHECK_NFSFT) 

TESTS = $(check_PROGRAMS)

check_bessel_i0_SOURCES = check_bessel_i0.c
check_bessel_i0_LDADD = $(top_builddir)/libnfft3.la

check_nfft_1d_SOURCES = check_nfft_1d.c
check_nfft_1d_LDADD = libcheck_nfft.la $(top_builddir)/libnfft3.la

check_nfft_2d_SOURCES = check_nfft_2d.c
check_nfft_2d_LDADD = $(top_builddir)/libnfft3.la

EXTRA_DIST=check_nfft.c check_nfft.h