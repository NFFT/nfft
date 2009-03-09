if HAVE_NFCT
  LIB_NFCT=nfct/libnfct.la
  DIR_NFCT=nfct
else
  LIB_NFCT=
endif

if HAVE_NFST
  LIB_NFST=nfst/libnfst.la
  DIR_NFST=nfst
else
  LIB_NFST=
  DIR_NFST=
endif

if HAVE_NFSFT
  LIB_NFSFT=nfsft/libnfsft.la
  DIR_NFSFT=nfsft
else
  LIB_NFSFT=
  DIR_NFSFT=
endif

if HAVE_NFSOFT
  LIB_NFSOFT=nfsoft/libnfsoft.la
  DIR_NFSOFT=nfsoft
else
  LIB_NFSOFT=
  DIR_NFSOFT=
endif

if HAVE_NNFFT
  LIB_NNFFT=nnfft/libnnfft.la
  DIR_NNFFT=nnfft
else
  LIB_NNFFT=
  DIR_NNFFT=
endif

if HAVE_NSFFT
  LIB_NSFFT=nsfft/libnsfft.la
  DIR_NSFFT=nsfft
else
  LIB_NSFFT=
  DIR_NSFFT=
endif

if HAVE_MRI
  LIB_MRI=mri/libmri.la
  DIR_MRI=mri
else
  LIB_MRI=
  DIR_MRI=
endif

if HAVE_FPT
  LIB_FPT=fpt/libfpt.la
  DIR_FPT=fpt
else
  LIB_FPT=
  DIR_FPT=
endif

SUBDIRS = $(DIR_FPT) $(DIR_MRI) $(DIR_NFCT) nfft $(DIR_NFSFT) $(DIR_NFSOFT) $(DIR_NFST) \
  $(DIR_NNFFT) $(DIR_NSFFT) solver .

noinst_LTLIBRARIES = libkernel.la

libkernel_la_SOURCES =

libkernel_la_LIBADD = $(LIB_FPT) $(LIB_MRI) $(LIB_NFCT) nfft/libnfft.la \
  $(LIB_NFSFT) $(LIB_NFSOFT) $(LIB_NFST) $(LIB_NNFFT) $(LIB_NSFFT) solver/libsolver.la
