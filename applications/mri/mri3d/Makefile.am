AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

applicationdir = $(pkgdatadir)/applications/mri/mri3d

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

## The program name
application_PROGRAMS =	construct_data_2d1d \
			construct_data_3d \
			reconstruct_data_gridding \
			reconstruct_data_2d1d \
			reconstruct_data_3d

## The source files.
construct_data_2d1d_SOURCES = construct_data_2d1d.c
construct_data_2d1d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

construct_data_3d_SOURCES = construct_data_3d.c
construct_data_3d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_gridding_SOURCES = reconstruct_data_gridding.c
reconstruct_data_gridding_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3 

reconstruct_data_2d1d_SOURCES = reconstruct_data_2d1d.c
reconstruct_data_2d1d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_3d_SOURCES = reconstruct_data_3d.c
reconstruct_data_3d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
