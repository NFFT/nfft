AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = construct_data_2d_1d construct_data_3d reconstruct_data_gridding reconstruct_data_2d_1d reconstruct_data_3d

## The source files.
construct_data_2d_1d_SOURCES = construct_data_2d_1d.c

construct_data_2d_1d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

construct_data_3d_SOURCES = construct_data_3d.c

construct_data_3d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_gridding_SOURCES = reconstruct_data_gridding.c

reconstruct_data_gridding_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_2d_1d_SOURCES = reconstruct_data_2d_1d.c

reconstruct_data_2d_1d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_3d_SOURCES = reconstruct_data_3d.c

reconstruct_data_3d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
