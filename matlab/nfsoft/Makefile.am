# $Id$

#SUBDIRS = @f_hat

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
nfsoftmatlabdir = $(datadir)/nfft/matlab/nfsoft

# library
lib_LTLIBRARIES = libnfsoft.la
libnfsoft_la_SOURCES = nfsoftmex.c
libnfsoft_la_LIBADD = $(top_builddir)/libnfft3_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)
libnfsoft_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
libnfsoft_la_CFLAGS = $(OPENMP_CFLAGS)
endif


dist_nfsoftmatlab_DATA = Contents.m \
	nfsoft_get_num_threads.m

# target all-am builds .libs/libnfsoft@matlab_mexext@
nfsoftmex@matlab_mexext@: all-am 
	rm -f ./nfsoftmex@matlab_mexext@
	$(LN_S) .libs/libnfsoft@matlab_mexext@ ./nfsoftmex@matlab_mexext@

all: all-am nfsoftmex@matlab_mexext@

clean: clean-am
	rm -f ./nfsoftmex@matlab_mexext@

install-exec-hook:
	rm -f $(libdir)/nfsoftmex@matlab_mexext@
	$(LN_S) $(libdir)/libnfsoft@matlab_mexext@ $(libdir)/nfsoftmex@matlab_mexext@

uninstall-hook:
	rm -f $(libdir)/nfsoftmex@matlab_mexext@

