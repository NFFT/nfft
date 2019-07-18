AM_CPPFLAGS = -I$(top_srcdir)/include -DSRCDIR=@abs_srcdir@

SUBDIRS= data

if HAVE_CUNIT
if HAVE_THREADS
if HAVE_OPENMP
  CHECK_THREADS = checkall_threads
else
  CHECK_THREADS =
endif
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

checkall_SOURCES = check.c util.c util.h reflect.c reflect.h bspline.c bspline.h bessel.c bessel.h nfft.c nfft.h $(NFCT_SOURCES) $(NFST_SOURCES)
checkall_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la -lm -lcunit

if HAVE_THREADS
if HAVE_OPENMP
  checkall_threads_SOURCES = $(checkall_SOURCES)
  checkall_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@ -lcunit
  checkall_threads_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

clean-local:
	rm -f CUnitAutomated-Results.xml CUnitAutomated_threads-Results.xml checkall.log checkall.trs checkall_threads.log checkall_threads.trs
