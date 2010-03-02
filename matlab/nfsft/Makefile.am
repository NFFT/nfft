# $Id$

SUBDIRS = @f_hat

# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# application directory
applicationdir = $(pkgdatadir)/matlab/nfsft

# library
application_LTLIBRARIES = libnfsft.la
libnfsft_la_SOURCES = nfsftmex.c
libnfsft_la_LIBADD = $(top_builddir)/libnfft3.la $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)
libnfsft_la_LDFLAGS = -module -shared -shrext $(matlab_mexext) -avoid-version $(matlab_LDFLAGS)

EXTRA_DIST = cc.m Contents.m f_hat_index.m ndsft_adjoint.m ndsft_trafo.m \
  nfsft_adjoint.m nfsft_finalize.m nfsft_forget.m nfsft_get_f_hat.m \
  nfsft_get_f.m nfsft_get_x.m nfsft_init_advanced.m nfsft_init_guru.m \
  nfsft_init.m NFSFT_NO_DIRECT_ALGORITHM.m NFSFT_NO_FAST_ALGORITHM.m \
  NFSFT_NORMALIZED.m nfsft_precompute_x.m nfsft_precompute.m \
  NFSFT_PRESERVE_F_HAT.m nfsft_set_f_hat.m nfsft_set_f.m nfsft_set_x.m \
  nfsft_trafo.m NFSFT_USE_DPT.m NFSFT_USE_NDFT.m nfsftmex.m gl.m simple_test.m \
  nfsft_get_f_hat_linear.m nfsft_set_f_hat_linear.m projection.m

application_DATA = $(EXTRA_DIST)

nfsftmex@matlab_mexext@: .libs/libnfsft@matlab_mexext@
	rm -f ./nfsftmex@matlab_mexext@
	$(LN_S) .libs/libnfsft@matlab_mexext@ ./nfsftmex@matlab_mexext@

all: all-am nfsftmex@matlab_mexext@

clean: clean-am
	rm -f ./nfsftmex@matlab_mexext@

install-exec-hook:
	$(LN_S) $(applicationdir)/libnfsft@matlab_mexext@ $(applicationdir)/nfsftmex@matlab_mexext@
