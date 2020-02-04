AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  FASTSUM_TEST_THREADS=fastsum_test_threads
  LIBFASTSUM_THREADS=libfastsum_threads.la
else
  FASTSUM_TEST_THREADS=
  LIBFASTSUM_THREADS=
endif

if HAVE_THREADS
if HAVE_OPENMP
  FASTSUM_BENCHOMP_PROGS=fastsum_benchomp fastsum_benchomp_createdataset fastsum_benchomp_detail_single fastsum_benchomp_detail_threads
else
  FASTSUM_BENCHOMP_PROGS=
endif
else
  FASTSUM_BENCHOMP_PROGS=
endif

noinst_LTLIBRARIES = libfastsum.la libkernels.la $(LIBFASTSUM_THREADS)
noinst_PROGRAMS = fastsum_test fastsum_matlab $(FASTSUM_TEST_THREADS) $(FASTSUM_BENCHOMP_PROGS)

libfastsum_la_SOURCES = fastsum.c fastsum.h kernels.h
libfastsum_la_LIBADD = @fftw3_LIBS@
libfastsum_la_LDFLAGS = @fftw3_LDFLAGS@

libkernels_la_SOURCES = kernels.c

if HAVE_THREADS
  libfastsum_threads_la_SOURCES = fastsum.c fastsum.h kernels.h
  libfastsum_threads_la_LIBADD = @fftw3_threads_LIBS@
  libfastsum_threads_la_LDFLAGS = @fftw3_threads_LDFLAGS@
if HAVE_OPENMP
  libfastsum_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

fastsum_test_SOURCES = fastsum_test.c
fastsum_test_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

fastsum_matlab_SOURCES = fastsum_matlab.c
if HAVE_THREADS
fastsum_matlab_LDADD = libfastsum_threads.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
else
fastsum_matlab_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@
endif

if HAVE_THREADS
  fastsum_test_threads_SOURCES = fastsum_test.c
  fastsum_test_threads_LDADD = libfastsum_threads.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
if HAVE_OPENMP
  fastsum_test_threads_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

if HAVE_THREADS
if HAVE_OPENMP
  fastsum_benchomp_SOURCES = fastsum_benchomp.c
  fastsum_benchomp_LDADD = libfastsum_threads.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
  fastsum_benchomp_CFLAGS = $(OPENMP_CFLAGS)

  fastsum_benchomp_createdataset_SOURCES = fastsum_benchomp_createdataset.c
  fastsum_benchomp_createdataset_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

  fastsum_benchomp_detail_single_SOURCES = fastsum_benchomp_detail.c
  fastsum_benchomp_detail_single_LDADD = libfastsum.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LDFLAGS@ @fftw3_LIBS@

  fastsum_benchomp_detail_threads_SOURCES = fastsum_benchomp_detail.c
  fastsum_benchomp_detail_threads_LDADD = libfastsum_threads.la libkernels.la $(top_builddir)/libnfft3@PREC_SUFFIX@_threads.la @fftw3_LDFLAGS@ @fftw3_threads_LIBS@
  fastsum_benchomp_detail_threads_CFLAGS = $(OPENMP_CFLAGS)
endif
endif

EXTRA_DIST = fastsum.m fastsum_test.m README
