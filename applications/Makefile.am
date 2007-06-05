#DIST_SUBDIRS = fastgauss fastsum fastsumS2 mri polarFFT quadratureS2 radon 
DIST_SUBDIRS = fastgauss fastsum fastsumS2 mri polarFFT radon 

if MAINTAINER_MODE 
  EXTRA_DIRS= texture iterS2
else
  EXTRA_DIRS=   
endif 

SUBDIRS= $(DIST_SUBDIRS) $(EXTRA_DIRS)

EXTRA_DIST = doxygen.c
