#!/bin/sh
echo "Running MATLAB nfft tests..."
matlab -nodisplay -r "nfftUnitTestsRunAndExit;"
