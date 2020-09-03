# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/applications/fastsum -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
fastsummatlabdir = $(datadir)/nfft/matlab/fastsum

# library
lib_LTLIBRARIES = libfastsum.la
libfastsum_la_SOURCES = fastsummex.c

libfastsum_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)

libfastsum_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libfastsum_la_CFLAGS = $(OPENMP_CFLAGS)
endif


dist_fastsummatlab_DATA = Contents.m EXACT_NEARFIELD.m NEARFIELD_BOXES.m fastsum.m fastsum_get_num_threads.m fastsum_init.m fastsum_set_x.m fastsum_set_alpha.m fastsum_set_y.m fastsum_trafo.m fastsum_trafo_direct.m fastsum_get_f.m fastsum_finalize.m simple_test.m test_fastsum.m fastsum_set_num_threads.m

# target all-am builds .libs/libfastsum@matlab_mexext@
fastsummex@matlab_mexext@: all-am
	rm -f ./fastsummex@matlab_mexext@
	$(LN_S) .libs/libfastsum@matlab_mexext@ ./fastsummex@matlab_mexext@

all: all-am fastsummex@matlab_mexext@

clean: clean-am
	rm -f ./fastsummex@matlab_mexext@

install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/fastsummex@matlab_mexext@
	$(LN_S) $(DESTDIR)$(libdir)/libfastsum@matlab_mexext@ $(DESTDIR)$(libdir)/fastsummex@matlab_mexext@

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/fastsummex@matlab_mexext@

