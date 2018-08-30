if HAVE_OCTAVE
if HAVE_NFSFT
TESTS = check_nfft_octave.sh check_nfsft_octave.sh
else
TESTS = check_nfft_octave.sh 
endif
else
if HAVE_MATLAB
if HAVE_NFSFT
TESTS = check_nfft_matlab.sh check_nfsft_matlab.sh
else
TESTS = check_nfft_matlab.sh
endif
endif
endif
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
 nfsftUnitTestsRunAndExit.m

clean-local:
	rm -f check_nfft_octave.log check_nfft_octave.trs check_nfft_matlab.log check_nfft_matlab.trs check_nfft_matlab.output check_nfsft_octave.log check_nfsft_octave.trs check_nfsft_matlab.log check_nfsft_matlab.trs check_nfsft_matlab.output
