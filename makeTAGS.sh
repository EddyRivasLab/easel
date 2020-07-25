#! /bin/sh

# `-x`: Higher-level ~/src/makeTAGS.sh calls us with
# this option, for indexing entire stack of lab code (easel, hmmer3,
# hmmer4, infernal) without redundancy. Output TAGS file
# to TAGS.part, and the higher level script assembles.
#

while [[ "$#" -gt 0 ]]; do case $1 in
  -x) optx=1;;
  *) echo "Unknown option: $1"; exit 1;;
esac; shift; done

if [ $optx ]; then
    out="-o TAGS.part"
else
    out=""
fi    

etags    $out configure.ac
etags -a $out LICENSE

find . -name "*.c"       -print -or -name "*.h"  -print | xargs etags -a $out
find . -name "*.pl"      -print -or -name "*.pm" -print | xargs etags -a $out
find . -name "*.py"      -print                         | xargs etags -a $out
find . -name "*.sh"      -print                         | xargs etags -a $out
find . -name "*.md"      -print                         | xargs etags -a $out
find . -name "*.tex"     -print                         | xargs etags -a $out
find . -name "*.man"     -print                         | xargs etags -a $out
find . -name "*.in"      -print                         | xargs etags -a $out
find . -name "*.sqc"     -print                         | xargs etags -a $out
find . -name "*README"   -print                         | xargs etags -a $out
