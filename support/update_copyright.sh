#!/bin/sh
function warn {
echo "warn: $1 does not seem to have the correct copyright format"
}

function replace {
  m=$(sed -e "/::Package::/!d" -e "/::Package::/=" -e"/::Package::/d" < $3)
  if test -n "$m"; then
    echo "warn: $3 skipped (Mathematica Package)"
  else
    # first occurence of start of copyright block
    s=$(sed -e "$1!d" -e "$1=" -e"$1d" < $3)
    # first occurence of end of copyright block
    e=$(sed -e "$2!d" -e "$2=" -e"$2d" < $3)
    # Check if copyright block delimters found.
    if test -n "$s" -a -n "$e"; then
      # Check if copyright delimters are ordered properly.
      if test "$s" -le "$e"; then
        # Replace copyright block with template.
        sed -e "${e}r$4" -e "${s},${e}d" < $3 > $3.tmp
        mv $3.tmp $3
        echo "ok: $3"
      else
        warn $3
      fi
    else
      if test "$7" == "yes"; then
        # Add copyright block to beginning of file.
        echo "$5" > head.txt
        echo "$6" > tail.txt
        cat head.txt $4 tail.txt $3 > $3.tmp
        mv $3.tmp $3
        rm head.txt
        rm tail.txt
        echo "add: $3"
      else
        echo "warn: $3 text from $4 not added"
      fi
    fi
  fi
}

# C source and header files
for name in $(find .. -wholename "../applications/texture" -prune -o -name "cycle.h" -prune -o -name "config.h" -prune -o -name "*.[ch]" -print -o -name "nfftconf.h.in" -print); do
  replace "/^ \* Copyright/" "/^ \* Franklin Street/" $name "copyright.txt" '/*' ' */' "yes"
done

# MATLAB scripts
for name in $(find .. -wholename "../applications/texture" -prune -o -name "Contents.m" -prune -o -name "*.m" -print); do
  replace "/^% Copyright/" "/^% Franklin Street/" $name "copyright_matlab.txt" "" "" "yes"
  replace "/^%   Copyright/" "/^%   Copyright/" $name "copyright_matlab_single_line.txt" "" "" "no"
done
