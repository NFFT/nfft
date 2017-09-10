if HAVE_NFSFT
  DIR_FASTSUMS2=fastsumS2
  DIR_QUADRATURES2=quadratureS2
  DIR_ITERS2=
else
  DIR_FASTSUMS2=
  DIR_QUADRATURES2=
  DIR_ITERS2=
endif

if HAVE_MRI
  DIR_MRI=mri
else
  DIR_MRI=
endif

DIST_SUBDIRS = fastgauss fastsum fastsumS2 mri polarFFT \
  quadratureS2 radon 

SUBDIRS= fastgauss fastsum $(DIR_FASTSUMS2) $(DIR_MRI) polarFFT \
  $(DIR_QUADRATURES2) radon $(DIR_ITERS2)

EXTRA_DIST = doxygen.c
