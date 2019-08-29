AM_CPPFLAGS = -I$(top_srcdir)/include

if HAVE_THREADS
  LIBUTIL_THREADS_LA = libutil_threads.la
else
  LIBUTIL_THREADS_LA =
endif

noinst_LTLIBRARIES = libutil.la $(LIBUTIL_THREADS_LA)
libutil_la_SOURCES = malloc.c sinc.c lambda.c bessel_i0.c float.c int.c error.c bspline.c assert.c sort.c rand.c vector1.c vector2.c vector3.c print.c damp.c thread.c time.c window.c version.c
# Unused file: voronoi.c

if HAVE_THREADS
  libutil_threads_la_SOURCES = $(libutil_la_SOURCES)
if HAVE_OPENMP
  libutil_threads_la_CFLAGS = $(OPENMP_CFLAGS)
endif
endif
