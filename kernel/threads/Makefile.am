AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/kernel/deconv \
              -I$(top_srcdir)/kernel/conv @fftw3_CPPFLAGS@

noinst_LTLIBRARIES = libthreads_ng.la

libthreads_ng_la_SOURCES = api.c conf.c deconv.c conv.c
libthreads_ng_la_CFLAGS = $(OPENMP_CFLAGS)
