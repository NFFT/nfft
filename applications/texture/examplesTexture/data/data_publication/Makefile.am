AM_CFLAGS = -Wall
AM_CPPFLAGS = -I../util -I$(top_srcdir)/include -I../../.. \
							-I../../texture_util

noinst_PROGRAMS = create compareSamples

create_SOURCES = create.c
compareSamples_SOURCES = compareSamples.c

create_LDADD = ../../texture_util/libtexture_util.la ../util/libdata_util.la
compareSamples_LDADD = ../../texture_util/libtexture_util.la
