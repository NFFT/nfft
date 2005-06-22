AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = construct_data construct_data_inhomogen reconstruct_data_gridding \
          reconstruct_data_2d reconstruct_data_special \
          reconstruct_data_special_invers reconstruct_data_inhomog_nnfft

## The source files.
construct_data_SOURCES = construct_data.c

construct_data_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

construct_data_inhomogen_SOURCES = construct_data_inhomogen.c

construct_data_inhomogen_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_gridding_SOURCES = reconstruct_data_gridding.c

reconstruct_data_gridding_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_2d_SOURCES = reconstruct_data_2d.c

reconstruct_data_2d_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_special_SOURCES = reconstruct_data_special.c

reconstruct_data_special_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_special_invers_SOURCES = reconstruct_data_special_invers.c

reconstruct_data_special_invers_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3

reconstruct_data_inhomog_nnfft_SOURCES = reconstruct_data_inhomog_nnfft.c

reconstruct_data_inhomog_nnfft_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3