#!/bin/sh
function replace {
  sed "1,/^\/\* \$Id/d" < $1 > $1.tmp
  if test "x$(stat -c%s $1.tmp)" = "x0"; then
    echo "Warning: $1 does not seem to have the correct header format."
  else
    cat copyright.txt $1.tmp > $1
  fi
  rm -f $1.tmp
}

for name in $(find .. -wholename "../applications/texture" -prune -o -name "*.[ch]" -print -o -name "nfftconf.h.in" -print); do
  replace $name
done

