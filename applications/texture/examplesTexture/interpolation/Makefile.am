AM_CFLAGS = -Wall -I$(top_srcdir)/include -I../texture_util -I../..
AM_CPPFLAGS = -D DEBUG_RESIDUUM -D DEBUG_CONVERGENCE_QUALITY -D DEBUG_TRUE_RESIDUUM

## The program name
nodist_PROGRAMS = calculateOmega calculateX validate

## The source files for the library
calculateOmega_SOURCES = calculateOmega.c
calculateX_SOURCES = calculateX.c
validate_SOURCES = validate.c

## Additional libraries for the library
calculateOmega_LDADD = ../texture_util/libtexture_util.la
calculateX_LDADD = ../texture_util/libtexture_util.la
validate_LDADD = ../texture_util/libtexture_util.la

nodistdir = .  
