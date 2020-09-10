# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
nfstmatlabdir = $(datadir)/nfft/matlab/nfst

# library
lib_LTLIBRARIES = libnfst.la
libnfst_la_SOURCES = nfstmex.c

libnfst_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)

libnfst_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libnfst_la_CFLAGS = $(OPENMP_CFLAGS)
endif

dist_nfstmatlab_DATA = Contents.m FFTW_ESTIMATE.m FFTW_MEASURE.m \
	FFT_OUT_OF_PLACE.m FG_PSI.m PRE_FG_PSI.m PRE_FULL_PSI.m PRE_LIN_PSI.m \
	PRE_PHI_HUT.m PRE_PSI.m ndst_adjoint.m ndst_trafo.m nfst.m \
	nfst_adjoint.m nfst_finalize.m nfst_get_f.m nfst_get_f_hat.m \
	nfst_get_num_threads.m nfst_get_x.m nfst_init_1d.m nfst_init_2d.m \
	nfst_init_3d.m nfst_init_guru.m nfst_set_f.m nfst_set_f_hat.m \
	nfst_set_x.m nfst_trafo.m simple_test.m test_nfst1d.m test_nfst2d.m \
	test_nfst3d.m nfst_set_num_threads.m

# target all-am builds .libs/libnfst@matlab_mexext@
nfstmex@matlab_mexext@: all-am
	rm -f ./nfstmex@matlab_mexext@
	$(LN_S) .libs/libnfst@matlab_mexext@ ./nfstmex@matlab_mexext@

all: all-am nfstmex@matlab_mexext@

clean: clean-am
	rm -f ./nfstmex@matlab_mexext@

install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/nfstmex@matlab_mexext@
	$(LN_S) $(DESTDIR)$(libdir)/libnfst@matlab_mexext@ $(DESTDIR)$(libdir)/nfstmex@matlab_mexext@

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/nfstmex@matlab_mexext@

