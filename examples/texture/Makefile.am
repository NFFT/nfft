AM_CFLAGS = -Wall -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The library name
noinst_LTLIBRARIES = libtexture.la

## The source files for the library
libtexture_la_SOURCES = texture.c

## Additional libraries for the library
#libtexture_la_LDADD =

## The program name
noinst_PROGRAMS = textureTest

## The source files for the library
textureTest_SOURCES = textureTest.c

## Additional libraries for the library
textureTest_LDADD = libtexture.la
