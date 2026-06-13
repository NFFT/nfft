# MATLAB mex host smoke test — run once on a host with a licensed MATLAB.
export MATLABROOT=/path/to/MATLAB/R20XXx

# (A) SERIAL — the must-pass case.
cmake -S . -B build-mat -DNFFT_WITH_MATLAB="$MATLABROOT" -DNFFT_ENABLE_MATLAB_THREADS=OFF
cmake --build build-mat -j --target matlab_nfftmex
cd build-mat/matlab/nfft && "$MATLABROOT/bin/matlab" -batch "simple_test"; cd -
# PASS = nfftmex.<mexext> built AND simple_test prints an error norm ~1e-15 with
#        no 'Undefined function' / 'Invalid MEX-file' error.
#        Also confirm: ldd build-mat/matlab/nfft/nfftmex.<ext> resolves libmwfftw3.so.3
#        against MATLAB's own copy (the versioned-name path on real MATLAB).

# (B) THREADED — the known-hazard case (libgomp vs MATLAB's libiomp5/MKL).
cmake -S . -B build-matT -DNFFT_WITH_MATLAB="$MATLABROOT" -DNFFT_ENABLE_OPENMP=ON -DNFFT_ENABLE_MATLAB_THREADS=ON
cmake --build build-matT -j --target matlab_nfftmex
cd build-matT/matlab/nfft && "$MATLABROOT/bin/matlab" -batch "simple_test"; cd -
# If (B) hangs/crashes on load → the OpenMP-runtime clash is real on this MATLAB;
# record it and recommend NFFT_ENABLE_MATLAB_THREADS=OFF until libiomp5-aware linking lands.
