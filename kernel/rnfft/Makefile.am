AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_OPENMP
  LIBRNFFT_THREADS_LA = librnfft_threads.la
else
  LIBRNFFT_THREADS_LA =
endif

noinst_LTLIBRARIES = librnfft.la $(LIBRNFFT_THREADS_LA)
librnfft_la_SOURCES = rnfft.c

if ENABLE_OPENMP
librnfft_threads_la_SOURCES = rnfft.c
librnfft_threads_la_CFLAGS = $(OPENMP_CFLAGS)
librnfft_threads_la_LDLAGS = $(OPENMP_CFLAGS)
endif
