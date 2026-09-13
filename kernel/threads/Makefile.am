AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/kernel/deconv @fftw3_CPPFLAGS@

noinst_LTLIBRARIES = libthreads_ng.la

libthreads_ng_la_SOURCES = api.c conf.c deconv.c
libthreads_ng_la_CFLAGS = $(OPENMP_CFLAGS)
