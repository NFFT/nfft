#!/bin/sh
# Bootstrap + configure + build the Autotools tree, to check the Makefile.am
# wiring for tune.c / tune_ng.c.
set -e
cd "$1"
./bootstrap.sh
./configure --enable-all --enable-tests --with-window=kaiserbessel
make -j8
echo "AUTOTOOLS BUILD OK"
