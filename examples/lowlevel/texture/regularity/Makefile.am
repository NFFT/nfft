AM_CFLAGS = -I$(top_srcdir)/include

## The program name
noinst_PROGRAMS = regularityTest

## The source files for the library
regularityTest_SOURCES = regularityTest.c

## Additional libraries for the library
regularityTest_LDADD = $(top_builddir)/libnfft3.la -lm -lfftw3
