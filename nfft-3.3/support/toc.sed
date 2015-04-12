#!/usr/bin/sed -f
s/^ \* \\section \([a-zA-Z0-9_]*\) \(.*\)/<ul><li><a href="#\1">\2<\/a><\/li><\/ul>/w "${srcdir}"/doc/api/html/toc.txt
s/^ \* \\subsection \([a-zA-Z0-9_]*\) \(.*\)/<ul><ul><li><a href="#\1">\2<\/a><\/li><\/ul><\/ul>/w "${srcdir}"/doc/api/html/toc.txt
s/^ \* \\subsubsection \([a-zA-Z0-9_]*\) \(.*\)/<ul><ul><ul><li><a href="#\1">\2<\/a><\/li><\/ul><\/ul><\/ul>/w "${srcdir}"/doc/api/html/toc.txt

