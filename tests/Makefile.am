AM_CPPFLAGS = -I$(top_srcdir)/include @cunit_CPPFLAGS@ -DSRCDIR=@abs_srcdir@

SUBDIRS= data

if HAVE_CUNIT
if ENABLE_OPENMP
CHECK_THREADS = checkall_threads checkall_ng_threads
else
CHECK_THREADS =
endif
else
CHECK_THREADS =
endif

if HAVE_CUNIT
  CHECK=checkall checkall_ng
else
  CHECK=
endif

if ENABLE_NG_OMP
CHECK_NG_OMP = checkall_ngomp
else
CHECK_NG_OMP =
endif

check_PROGRAMS = $(CHECK) $(CHECK_THREADS) $(CHECK_NG_OMP)

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

# checkall: legacy nfft/nfct/nfst suites + shared util infrastructure.
checkall_SOURCES = check.c util.c util.h accuracy_log.c accuracy_log.h reflect.c reflect.h bspline.c bspline.h bessel.c bessel.h window.c window.h nfft.c nfft.h $(NFCT_SOURCES) $(NFST_SOURCES)
checkall_LDFLAGS = @fftw3_LDFLAGS@ @cunit_LDFLAGS@
checkall_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @cunit_LIBS@

# checkall_ng: new planner (plan_ng / planner-native) suites + shared util.
checkall_ng_SOURCES = check_ng.c util.c util.h accuracy_log.c accuracy_log.h reflect.c reflect.h bspline.c bspline.h bessel.c bessel.h planner.c planner.h nplan.c nplan.h nplan_data.c nplan_perm.c nfast.c nfast.h
checkall_ng_LDFLAGS = @fftw3_LDFLAGS@ @cunit_LDFLAGS@
checkall_ng_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @cunit_LIBS@

# checkall_ngomp: the add-on threading library's own contract, linked against
# libnfft3<suffix>_ng_omp in addition to libnfft3<suffix>.
checkall_ngomp_SOURCES = check_ngomp.c ngomp.c ngomp.h
checkall_ngomp_CFLAGS = $(OPENMP_CFLAGS)
checkall_ngomp_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@ @cunit_LDFLAGS@
checkall_ngomp_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_ng_omp.la $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @fftw3_LIBS_omp@ @cunit_LIBS@

EXTRA_DIST = data/generated/nfft_testcases.h \
             data/generated/nfft_native_testcases.h \
             data/generated/nfct_testcases.h \
             data/generated/nfst_testcases.h

if ENABLE_OPENMP
checkall_threads_SOURCES = $(checkall_SOURCES)
checkall_threads_CPPFLAGS = -I$(top_srcdir)/include $(cunit_CPPFLAGS)
checkall_threads_CFLAGS = $(OPENMP_CFLAGS)
checkall_threads_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@ @cunit_LDFLAGS@
checkall_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@  $(OPENMP_LIBS) @cunit_LIBS@

checkall_ng_threads_SOURCES = $(checkall_ng_SOURCES)
checkall_ng_threads_CPPFLAGS = -I$(top_srcdir)/include $(cunit_CPPFLAGS)
checkall_ng_threads_CFLAGS = $(OPENMP_CFLAGS)
checkall_ng_threads_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@ @cunit_LDFLAGS@
checkall_ng_threads_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_omp.la @fftw3_LIBS_omp@ @fftw3_LIBS@  $(OPENMP_LIBS) @cunit_LIBS@
endif

clean-local:
	rm -f CUnitAutomated-Results.xml CUnitAutomated_threads-Results.xml checkall.log checkall.trs checkall_threads.log checkall_threads.trs CUnitAutomated_ng-Results.xml CUnitAutomated_ng_threads-Results.xml checkall_ng.log checkall_ng.trs checkall_ng_threads.log checkall_ng_threads.trs CUnitAutomated_ngomp-Results.xml checkall_ngomp.log checkall_ngomp.trs
