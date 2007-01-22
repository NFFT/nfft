#!/bin/sh
############################################################################ NOTE: If you just want to build NFFT 3, do not use this file. Just follow 
# the installation instructions as described in the tutorial found under
# doc/tutorial.
#
# This file is based on the bootstrap.sh script from fftw 3.1.2 by 
# M. Frigo and S. G. Johnson
###########################################################################

touch ChangeLog

echo "PLEASE IGNORE WARNINGS AND ERRORS"

# paranoia: sometimes autoconf doesn't get things right the first time
rm -rf autom4te.cache
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force

rm -f config.cache
