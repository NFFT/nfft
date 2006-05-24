AM_CFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)

## The program name
noinst_PROGRAMS = docking

## The source files.
docking_SOURCES = docking.c
docking_LDADD = $(top_builddir)/.libs/libnfft3.a -lm -lfftw3

EXTRA_DIST = Molekuel_A_komplett.txt molekuel_A_new.txt Molekuel_B_komplett.txt
