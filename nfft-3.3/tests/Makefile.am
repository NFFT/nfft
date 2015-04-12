# $Id: Makefile.am 3456 2010-03-17 19:25:47Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include -DSRCDIR=@abs_srcdir@

SUBDIRS= data

if HAVE_THREADS
if HAVE_OPENMP
  CHECK_THREADS = checkall_threads
else
  CHECK_THREADS =
endif
else
  CHECK_THREADS =
endif

if HAVE_CUNIT
  CHECK=checkall
else
  CHECK=
endif

check_PROGRAMS = $(CHECK) $(CHECK_THREADS)
#debug

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
checkall_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la -lm -lcunit -lncurses

if HAVE_THREADS
if HAVE_OPENMP
  checkall_threads_SOURCES = $(checkall_SOURCES)
  checkall_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@ -lcunit -lncurses
  checkall_threads_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

#debug_SOURCES = debug.c nfct.c nfst.c nfct.h nfst.h
#debug_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la -lcunit -lncurses

clean-local:
	rm -f CUnitAutomated-Results.xml CUnitAutomated_threads-Results.xml checkall.log checkall.trs checkall_threads.log checkall_threads.trs
