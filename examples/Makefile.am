if HAVE_NFCT
  DIR_NFCT=nfct
else
  LIB_NFCT=
endif

if HAVE_NFST
  DIR_NFST=nfst
else
  DIR_NFST=
endif

if HAVE_NFSFT
  DIR_NFSFT=nfsft
else
  DIR_NFSFT=
endif

if HAVE_NFSOFT
  DIR_NFSOFT=nfsoft
else
  DIR_NFSOFT=
endif

if HAVE_NNFFT
  DIR_NNFFT=nnfft
else
  DIR_NNFFT=
endif

if HAVE_NSFFT
  DIR_NSFFT=nsfft
else
  DIR_NSFFT=
endif

if HAVE_MRI
  DIR_MRI=mri
else
  DIR_MRI=
endif

if HAVE_FPT
  DIR_FPT=fpt
else
  DIR_FPT=
endif

SUBDIRS = $(DIR_FPT) $(DIR_MRI) $(DIR_NFCT) nfft $(DIR_NFSFT) $(DIR_NFSOFT) $(DIR_NFST) \
  $(DIR_NNFFT) $(DIR_NSFFT) solver
 
EXTRA_DIST = doxygen.c
