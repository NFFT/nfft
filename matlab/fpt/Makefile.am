# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
fptmatlabdir = $(datadir)/nfft/matlab/fpt

# library
lib_LTLIBRARIES = libfpt.la
libfpt_la_SOURCES = fptmex.c
libfpt_la_LIBADD = $(top_builddir)/libnfft3_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)
libfpt_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libfpt_la_CFLAGS = $(OPENMP_CFLAGS)
endif


dist_fptmatlab_DATA = simple_test.m fpt_init.m fpt_precompute.m fpt_trafo.m fpt_trafo_direct.m fpt_transposed.m fpt_transposed_direct.m fpt_finalize.m eval_clenshaw.m \
	FPT_NO_FAST_ALGORITHM.m FPT_NO_DIRECT_ALGORITHM.m FPT_NO_STABILIZATION.m FPT_FUNCTION_VALUES.m FPT_AL_SYMMETRY.m

# target all-am builds .libs/libfpt@matlab_mexext@
fptmex@matlab_mexext@: all-am
	rm -f ./fptmex@matlab_mexext@
	$(LN_S) .libs/libfpt@matlab_mexext@ ./fptmex@matlab_mexext@

all: all-am fptmex@matlab_mexext@

clean: clean-am
	rm -f ./fptmex@matlab_mexext@

install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/fptmex@matlab_mexext@
	$(LN_S) $(DESTDIR)$(libdir)/libfpt@matlab_mexext@ $(DESTDIR)$(libdir)/fptmex@matlab_mexext@

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/fptmex@matlab_mexext@

