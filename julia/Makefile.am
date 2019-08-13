if HAVE_NFCT
  DIR_NFCT=nfct
else
  DIR_NFCT=
endif

if HAVE_NFST
  DIR_NFST=nfst
else
  DIR_NFST=
endif

if HAVE_APPLICATIONS
  DIR_FASTSUM=fastsum
else
  DIR_FASTSUM=
endif


SUBDIRS = nfft $(DIR_NFCT) $(DIR_NFST) $(DIR_FASTSUM)

EXTRA_DIST = README.md
