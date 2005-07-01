#!/bin/bash
for file in $(find $1 -name *.c)
do
  strip.sed -i $file
done
