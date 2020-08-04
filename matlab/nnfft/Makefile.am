# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
nnfftmatlabdir = $(datadir)/nfft/matlab/nnfft

# library
lib_LTLIBRARIES = libnnfft.la
libnnfft_la_SOURCES = nnfftmex.c

libnnfft_la_LIBADD = $(top_builddir)/libnfft3_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)

libnnfft_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libnnfft_la_CFLAGS = $(OPENMP_CFLAGS)
endif


dist_nnfftmatlab_DATA = Contents.m nnfft.m nnfft_display.m nnfft_finalize.m \
	nnfft_get_f.m nnfft_get_f_hat.m nnfft_get_x.m nnfft_init.m nnfft_init_1d.m \
	nnfft_init_2d.m nnfft_init_3d.m nnfft_init_guru.m nnfft_precompute_psi.m \
	nnfft_set_f.m nnfft_set_f_hat.m nnfft_set_v.m nnfft_set_x.m nnfft_trafo.m \
	nnfft_trafo_direct.m PRE_FULL_PSI.m PRE_LIN_PSI.m PRE_PHI_HUT.m PRE_PSI.m \
	simple_test.m test_nnfft1d.m test_nnfft2d.m nnfft_get_num_threads.m \
        nnfft_set_num_threads.m

# target all-am builds .libs/libnnfft@matlab_mexext@
nnfftmex@matlab_mexext@: all-am
	rm -f ./nnfftmex@matlab_mexext@
	$(LN_S) .libs/libnnfft@matlab_mexext@ ./nnfftmex@matlab_mexext@

all: all-am nnfftmex@matlab_mexext@

clean: clean-am
	rm -f ./nnfftmex@matlab_mexext@

install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/nnfftmex@matlab_mexext@
	$(LN_S) $(DESTDIR)$(libdir)/libnnfft@matlab_mexext@ $(DESTDIR)$(libdir)/nnfftmex@matlab_mexext@

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/nnfftmex@matlab_mexext@


