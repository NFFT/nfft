AM_CFLAGS = "-Wall -I$(top_srcdir)/include -I../texture_util"
AM_CPPFLAGS = "-D DEBUG_RESIDUUM -D DEBUG_CONVERGENCE_QUALITY #-D DEBUG_TRUE_RESIDUUM"

## The program name
noinst_PROGRAMS = calculateOmega calculateX validate

## The source files for the library
calculateOmega_SOURCES = calculateOmega.c
calculateX_SOURCES = calculateX.c
validate_SOURCES = validate.c

## Additional libraries for the library
calculateOmega_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3 ../texture_util/libtexture_util.la
calculateX_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3 ../texture_util/libtexture_util.la
validate_LDADD = -lm ../texture_util/libtexture_util.la
