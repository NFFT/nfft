if HAVE_OCTAVE
  TESTS = check_nfft_octave.sh
else
if HAVE_MATLAB
  TESTS = check_nfft_matlab.sh
endif
endif

clean-local:
	rm -f check_nfft_octave.sh.log check_nfft_octave.sh.trs check_nfft_matlab.sh.log check_nfft_matlab.sh.trs
