AM_CFLAGS = -Wall -I$(top_srcdir)/include -I../../texture_util -I../../.. \
-I../../data/util 

## The program name
nodist_PROGRAMS = create

## The source files for the library
create_SOURCES = create.c

## Additional libraries for the library
create_LDADD = ../../texture_util/libtexture_util.la \
../../data/util/libdata_util.la

nodistdir = $(abs_srcdir)
