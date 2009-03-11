#!/bin/sh
for name in $(find . -name "*.[ch]"); do
  echo $name
done

# sed "1,/^\/\* \$Id/d" < nfsft.c > nfsft_sed.c
