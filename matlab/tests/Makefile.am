if HAVE_OCTAVE
MAT_OCT=octave
endif
if HAVE_MATLAB
MAT_OCT=matlab
endif

TEST_NFFT=check_nfft_$(MAT_OCT).sh 
if HAVE_NFSFT
TEST_NFSFT=check_nfsft_$(MAT_OCT).sh
else
TEST_NFSFT=
endif
if HAVE_NFSOFT
TEST_NFSOFT=check_nfsoft_$(MAT_OCT).sh
else
TEST_NFSOFT=
endif

TESTS = $(TEST_NFFT) $(TEST_NFSFT) $(TEST_NFSOFT)
TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = /bin/sh

matlabtestsdir = $(datadir)/nfft/matlab/tests
dist_matlabtests_DATA = nfftTestcaseCheckDelegate.m nfftTestcaseDelegateFile.m \
 nfftTestcaseDelegate.m nfftTestcaseDelegateOnline.m \
 nfftTestcaseInitDelegate.m nfftTestcaseTrafoDelegate.m \
 nfftUnitTestsRunAndExit.m \
 nfsftTestcaseCheckDelegate.m nfsftTestcaseDelegate.m \
 nfsftTestcaseDelegateGLQuadratureOnline.m nfsftTestcaseDelegateOnline.m \
 nfsftTestcaseInitDelegate.m nfsftTestcaseTrafoDelegate.m nfsftUnitTests.m \
 nfsftUnitTestsRunAndExit.m \
 nfsoftTestcaseCheckDelegate.m nfsoftTestcaseDelegate.m nfsoftTestcaseDelegateGLQuadratureOnline.m nfsoftTestcaseDelegateOnline.m nfsoftTestcaseInitDelegate.m nfsoftTestcaseTrafoDelegate.m nfsoftUnitTests.m nfsoftUnitTestsRunAndExit.m 

clean-local:
	rm -f check_nfft_octave.log check_nfft_octave.trs check_nfft_matlab.log check_nfft_matlab.trs check_nfft_matlab.output check_nfsft_octave.log check_nfsft_octave.trs check_nfsft_matlab.log check_nfsft_matlab.trs check_nfsft_matlab.output check_nfsoft_octave.log check_nfsoft_octave.trs check_nfsoft_matlab.log check_nfsoft_matlab.trs check_nfsoft_matlab.output
