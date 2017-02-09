# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/applications/fastsum -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
fastsummatlabdir = $(datadir)/nfft/matlab/fastsum

# library
lib_LTLIBRARIES = libfastsum.la
libfastsum_la_SOURCES = fastsummex.c

if HAVE_MATLAB_THREADS
libfastsum_la_LIBADD = $(top_srcdir)/applications/fastsum/libfastsum_threads.la $(top_srcdir)/applications/fastsum/libkernels.la  $(top_builddir)/libnfft3_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)
else
libfastsum_la_LIBADD = $(top_srcdir)/applications/fastsum/libfastsum.la $(top_srcdir)/applications/fastsum/libkernels.la  $(top_builddir)/libnfft3_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)
endif

libfastsum_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libfastsum_la_CFLAGS = $(OPENMP_CFLAGS)
endif


dist_fastsummatlab_DATA = fastsum_get_num_threads.m fastsum_init.m fastsum_set_x.m fastsum_set_y.m fastsum_trafo.m fastsum_trafo_direct.m fastsum_get_f.m fastsum_finalize.m

# target all-am builds .libs/libfastsum@matlab_mexext@
fastsummex@matlab_mexext@: all-am 
	rm -f ./fastsummex@matlab_mexext@
	$(LN_S) .libs/libfastsum@matlab_mexext@ ./fastsummex@matlab_mexext@

all: all-am fastsummex@matlab_mexext@

clean: clean-am
	rm -f ./fastsummex@matlab_mexext@

install-exec-hook:
	rm -f $(libdir)/fastsummex@matlab_mexext@
	$(LN_S) $(libdir)/libfastsum@matlab_mexext@ $(libdir)/fastsummex@matlab_mexext@

uninstall-hook:
	rm -f $(libdir)/fastsummex@matlab_mexext@

