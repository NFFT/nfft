# $Id$

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# application directory
applicationdir = $(pkgdatadir)/matlab

# library
lib_LTLIBRARIES = libnfft.la
libnfft_la_SOURCES = nfftmex.c
libnfft_la_LIBADD = $(top_builddir)/libnfft3noinst.la $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)
libnfft_la_LDFLAGS = -module -shared -shrext $(matlab_mexext) -avoid-version $(matlab_LDFLAGS)

EXTRA_DIST = Contents.m ndft_adjoint.m ndft_trafo.m nfft_adjoint.m nfft_trafo.m \
  nfft_get_f_hat.m nfft_get_f.m nfft_get_x.m \
  nfft_set_f_hat.m nfft_set_f.m nfft_set_x.m \
  nfft_init_1d.m nfft_precompute_psi.m nfft_finalize.m \
  nfftmex.m simple_test.m

application_DATA = $(EXTRA_DIST)

nfftmex@matlab_mexext@: .libs/libnfft@matlab_mexext@
	rm -f ./nfftmex@matlab_mexext@
	$(LN_S) .libs/libnfft@matlab_mexext@ ./nfftmex@matlab_mexext@

all: all-am nfftmex@matlab_mexext@

clean: clean-am
	rm -f ./nfftmex@matlab_mexext@
