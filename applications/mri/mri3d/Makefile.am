# $Id$

AM_CFLAGS = -I$(top_srcdir)/include

applicationdir = $(pkgdatadir)/applications/mri/mri3d

application_PROGRAMS = \
  construct_data_2d1d \
  construct_data_3d \
  reconstruct_data_gridding \
  reconstruct_data_2d1d \
  reconstruct_data_3d

construct_data_2d1d_SOURCES = construct_data_2d1d.c
construct_data_2d1d_LDADD = $(top_builddir)/libnfft3.la @FFTW3_LIBS@ @FFTW3_LIBCMD@

construct_data_3d_SOURCES = construct_data_3d.c
construct_data_3d_LDADD = $(top_builddir)/libnfft3.la

reconstruct_data_gridding_SOURCES = reconstruct_data_gridding.c
reconstruct_data_gridding_LDADD = $(top_builddir)/libnfft3.la @FFTW3_LIBS@ @FFTW3_LIBCMD@

reconstruct_data_2d1d_SOURCES = reconstruct_data_2d1d.c
reconstruct_data_2d1d_LDADD = $(top_builddir)/libnfft3.la @FFTW3_LIBS@ @FFTW3_LIBCMD@

reconstruct_data_3d_SOURCES = reconstruct_data_3d.c
reconstruct_data_3d_LDADD = $(top_builddir)/libnfft3.la

EXTRA_DIST = mri.m \
			construct_knots_linogram.m \
			construct_knots_radial_3d.m \
			construct_knots_radial.m \
			construct_knots_rose.m \
			construct_knots_spiral.m \
			construct_phantom.m \
			precompute_weights_2d.m \
			rms.m \
			visualize_data.m \
			README

application_DATA = $(EXTRA_DIST)
