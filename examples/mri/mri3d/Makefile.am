AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = construct_data reconstruct_data_gridding reconstruct_data_2_1d reconstruct_data_3d

## The source files.
construct_data_SOURCES = construct_data.c

construct_data_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_gridding_SOURCES = reconstruct_data_gridding.c

reconstruct_data_gridding_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_2_1d_SOURCES = reconstruct_data_2_1d.c

reconstruct_data_2_1d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_3d_SOURCES = reconstruct_data_3d.c

reconstruct_data_3d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
