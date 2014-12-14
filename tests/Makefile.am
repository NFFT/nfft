# $Id: Makefile.am 3456 2010-03-17 19:25:47Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include -DSRCDIR=@abs_srcdir@

SUBDIRS= data

if HAVE_CUNIT
  CHECK=checkall
else
  CHECK=
endif

check_PROGRAMS = $(CHECK) debug

TESTS = $(check_PROGRAMS)

if HAVE_NFCT
  NFCT_SOURCES=nfct.c nfct.h
else
  NFCT_SOURCES=
endif

if HAVE_NFST
  NFST_SOURCES=nfst.c nfst.h
else
  NFST_SOURCES=
endif

checkall_SOURCES = check.c bspline.c bspline.h bessel.c bessel.h nfft.c nfft.h $(NFCT_SOURCES) $(NFST_SOURCES)
checkall_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la -lcunit -lncurses

debug_SOURCES = debug.c nfct.c nfst.c nfct.h nfst.h
debug_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la -lcunit -lncurses

clean-local:
	rm -f CUnitAutomated-Results.xml
