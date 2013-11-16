# $Id: Makefile.am 3456 2010-03-17 19:25:47Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include -DSRCDIR=@abs_srcdir@

SUBDIRS= data

if HAVE_CUNIT
  CHECK=checkall
else
  CHECK=
endif

check_PROGRAMS =  $(CHECK)

TESTS = $(check_PROGRAMS)

checkall_SOURCES = check.c bspline.c bspline.h bessel.c bessel.h nfft.c nfft.h nfct.c nfct.h nfst.c nfst.h
checkall_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la -lcunit -lncurses

clean-local:
	rm -f CUnitAutomated-Results.xml
