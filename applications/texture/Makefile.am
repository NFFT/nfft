# Development code extra directories.
if WITH_DEVCODE
  nodist_SUBDIRS = examplesTexture
endif  

AM_CFLAGS = -Wall -I$(top_srcdir)/include -I$(top_srcdir)/kernel/solver 

## The library name
nodist_LTLIBRARIES = libtexture.la

## The source files for the library
libtexture_la_SOURCES = texture.c texture.h nfft3_texture.h

## Additional libraries for the library
libtexture_la_LIBADD = $(top_srcdir)/libnfft3.la

nodistdir = .