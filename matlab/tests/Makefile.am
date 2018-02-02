if HAVE_OCTAVE
  TESTS = check_nfft_octave.sh
else
if HAVE_MATLAB
  TESTS = check_nfft_matlab.sh
endif
endif
TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = /bin/sh

clean-local:
	rm -f check_nfft_octave.log check_nfft_octave.trs check_nfft_matlab.log check_nfft_matlab.trs
