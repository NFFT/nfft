#!/bin/sh
for name in $(find . -name "*.[ch]"); do
  sed "1,/^\/\* \$Id/d" < $name > $name.tmp
  if $($name.tmp) = "0"; then
    echo "Warning: $name does not seem to have the correct header format."
  else
    cat copyright.txt $name.tmp > $name
  fi
  rm -f $name.tmp
done

