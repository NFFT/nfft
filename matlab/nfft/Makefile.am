# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
nfftmatlabdir = $(datadir)/nfft/matlab/nfft

# library
lib_LTLIBRARIES = libnfft.la
libnfft_la_SOURCES = nfftmex.c

libnfft_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)

libnfft_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libnfft_la_CFLAGS = $(OPENMP_CFLAGS)
endif


dist_nfftmatlab_DATA = FFT_OUT_OF_PLACE.m FFTW_ESTIMATE.m FFTW_MEASURE.m FG_PSI.m Contents.m ndft_adjoint.m ndft_trafo.m \
	nfft_adjoint.m nfft_finalize.m nfft_get_f.m nfft_get_f_hat.m nfft_get_x.m nfft_init.m nfft_init_1d.m nfft_init_2d.m \
        nfft_init_3d.m nfft_init_guru.m nfft_precompute_psi.m nfft_set_f.m nfft_set_f_hat.m nfft_set_x.m nfft_trafo.m \
        PRE_FG_PSI.m PRE_FULL_PSI.m PRE_LIN_PSI.m PRE_PHI_HUT.m PRE_PSI.m simple_test.m \
	nfft_get_num_threads.m nfft.m test_nfft1d.m test_nfft2d.m test_nfft3d.m test_nfft4d.m \
	NFFT_OMP_BLOCKWISE_ADJOINT.m nfft_set_num_threads.m

# target all-am builds .libs/libnfft@matlab_mexext@
nfftmex@matlab_mexext@: all-am
	rm -f ./nfftmex@matlab_mexext@
	$(LN_S) .libs/libnfft@matlab_mexext@ ./nfftmex@matlab_mexext@

all: all-am nfftmex@matlab_mexext@

clean: clean-am
	rm -f ./nfftmex@matlab_mexext@

install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/nfftmex@matlab_mexext@
	$(LN_S) $(DESTDIR)$(libdir)/libnfft@matlab_mexext@ $(DESTDIR)$(libdir)/nfftmex@matlab_mexext@

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/nfftmex@matlab_mexext@

