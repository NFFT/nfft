#!/bin/sh
echo "Running Octave nfft tests..."
octave-cli --eval 'nfftUnitTestsRunAndExit;'
