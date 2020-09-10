# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
nfctmatlabdir = $(datadir)/nfft/matlab/nfct

# library
lib_LTLIBRARIES = libnfct.la
libnfct_la_SOURCES = nfctmex.c

libnfct_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)

libnfct_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libnfct_la_CFLAGS = $(OPENMP_CFLAGS)
endif

dist_nfctmatlab_DATA = Contents.m FFTW_ESTIMATE.m FFTW_MEASURE.m \
	FFT_OUT_OF_PLACE.m FG_PSI.m PRE_FG_PSI.m PRE_FULL_PSI.m PRE_LIN_PSI.m \
	PRE_PHI_HUT.m PRE_PSI.m ndct_adjoint.m ndct_trafo.m nfct.m \
	nfct_adjoint.m nfct_finalize.m nfct_get_f.m nfct_get_f_hat.m \
	nfct_get_num_threads.m nfct_get_x.m nfct_init_1d.m nfct_init_2d.m \
	nfct_init_3d.m nfct_init_guru.m nfct_set_f.m nfct_set_f_hat.m \
	nfct_set_x.m nfct_trafo.m simple_test.m test_nfct1d.m test_nfct2d.m \
	test_nfct3d.m nfct_set_num_threads.m

# target all-am builds .libs/libnfct@matlab_mexext@
nfctmex@matlab_mexext@: all-am
	rm -f ./nfctmex@matlab_mexext@
	$(LN_S) .libs/libnfct@matlab_mexext@ ./nfctmex@matlab_mexext@

all: all-am nfctmex@matlab_mexext@

clean: clean-am
	rm -f ./nfctmex@matlab_mexext@

install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/nfctmex@matlab_mexext@
	$(LN_S) $(DESTDIR)$(libdir)/libnfct@matlab_mexext@ $(DESTDIR)$(libdir)/nfctmex@matlab_mexext@

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/nfctmex@matlab_mexext@

