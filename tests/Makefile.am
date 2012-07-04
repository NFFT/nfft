# $Id: Makefile.am 3456 2010-03-17 19:25:47Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include

# check_nfft
if HAVE_CUNIT
  CHECK_NFFT=check_bspline
else
  CHECK_NFFT=
endif

#  TODO implement nfsft check
#    CHECK_NFSFT=check_nfsft
#    check_nfsft_SOURCES = check_nfsft.c
#    check_nfsft_LDADD = $(top_builddir)/libnfft3.la

check_PROGRAMS = $(CHECK_NFFT)

#check_bessel_i0 check_nfft_1d check_nfft_2d $(CHECK_NFSFT) 

TESTS = $(check_PROGRAMS)

#check_bessel_i0_SOURCES = check_bessel_i0.c
#check_bessel_i0_LDADD = $(top_builddir)/libnfft3.la

#check_nfft_SOURCES = check_nfft.c
#check_nfft_LDADD = $(top_builddir)/libnfft3.la -lcunit -lncurses

check_bspline_SOURCES = check_bspline.c
check_bspline_LDADD = $(top_builddir)/libnfft3.la -lcunit -lncurses
