#!/bin/sh
function replace {
  sed "$1" < $2 > $2.tmp
  if test "x$(stat -c%s $2.tmp)" = "x0"; then
    echo "Warning: $2 does not seem to have the correct header format."
  else
    cat copyright.txt $2.tmp > $2
  fi
  rm -f $2.tmp
}

# C source and header files
for name in $(find .. -wholename "../applications/texture" -prune -o -name "*.[ch]" -print -o -name "nfftconf.h.in" -print); do
  replace "1,/^\/\* \$Id/d" $name
done

# MATLAB scripts
for name in $(find .. -wholename "../applications/texture" -prune -o -name "*.m" -print); do
  replace "1,/^% \$Id/d" $name
done
