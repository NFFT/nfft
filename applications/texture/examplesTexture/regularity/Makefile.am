AM_CFLAGS = -Wall -I$(top_srcdir)/include -I../texture_util -I../..
AM_CPPFLAGS = -D DEBUG_RESIDUUM -D DEBUG_CONVERGENCE_QUALITY -D DEBUG_N # -D DEBUG_TRUE_RESIDUUM

## The program name
noinst_PROGRAMS = regularityTest

## The source files for the library
regularityTest_SOURCES = regularityTest.c

## Additional libraries for the library
regularityTest_LDADD = ../texture_util/libtexture_util.la
