# $Id$

# Set compiler flags
AM_CFLAGS = -I$(top_srcdir) $(MATLAB_CFLAGS)
# Set linker flags
AM_LDFLAGS = -L$(top_builddir) $(MATLAB_LDFLAGS)

# Set extension of linker output to extension for mex files
EXEEXT=$(MEXEXT)

# Application directory
applicationdir = $(pkgdatadir)/matlab

# Program name
application_PROGRAMS = nfsft

# Program sources
nfsft_SOURCES = nfsft.c ../malloc.c ../imex.h

# Additional libraries
nfsft_LDADD = $(top_builddir)/.libs/libnfft3.a $(MATLAB_LDADD)

EXTRA_DIST = Contents.m f_hat_index.m ndsft_adjoint.m ndsft_trafo.m nfsft_adjoint.m nfsft_finalize.m nfsft_forget.m nfsft_get_f_hat.m nfsft_get_f.m nfsft_get_x.m nfsft_init_advanced.m nfsft_init_guru.m nfsft_init.m NFSFT_NO_DIRECT_ALGORITHM.m NFSFT_NO_FAST_ALGORITHM.m NFSFT_NORMALIZED.m nfsft_precompute_x.m nfsft_precompute.m NFSFT_PRESERVE_F_HAT.m nfsft_set_f_hat.m nfsft_set_f.m nfsft_set_x.m nfsft_trafo.m NFSFT_USE_DPT.m NFSFT_USE_NDFT.m nfsft.m

application_DATA = $(EXTRA_DIST)

if HAVE_WINDOWS
all: all-am nfsft.dll.hack

nfsft.dll.hack: nfsft.dll.exe
	mv nfsft.dll.exe nfsft.dll

clean: clean-am
	rm -f nfsft.dll
	rm -f nfsft.dll.exe
endif

