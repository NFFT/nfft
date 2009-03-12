#!/bin/sh
function replace {
  sed -e "$1" -e "$2" < $3 > $3.tmp
  if test "x$(stat -c%s $3.tmp)" = "x0"; then
    echo "Warning: $3 does not seem to have the correct header format."
    rm -f $3.tmp
  else
    mv $3.tmp $3
  fi
}

# C source and header files
for name in $(find .. -wholename "../applications/texture" -prune -o -name "*.[ch]" -print -o -name "nfftconf.h.in" -print); do
  replace "/^ \* Franklin Street/r copyright.txt" "/^ \* Copyright/,/^ \* Franklin Street/d" $name
done

# MATLAB scripts
for name in $(find .. -wholename "../applications/texture" -prune -o -name "*.m" -print); do
  replace "/^% Franklin Street/r copyright_matlab.txt" "/^% Copyright/,/^% Franklin Street/d" $name
done
