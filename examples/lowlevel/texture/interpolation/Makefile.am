AM_CFLAGS = -Wall -I$(top_srcdir)/include
AM_CPPFLAGS = -D DEBUG_RESIDUUM -D DEBUG_CONVERGENCE_QUALITY

## The program name
noinst_PROGRAMS = calculateOmega calculateX validate

## The source files for the library
calculateOmega_SOURCES = calculateOmega.c
calculateX_SOURCES = calculateX.c
validate_SOURCES = validate.c

## Additional libraries for the library
calculateOmega_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
calculateX_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
validate_LDADD = -lm
