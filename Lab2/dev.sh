#!/bin/bash
set -ex
DEBUG=1 make -j
in=${1-./3d.txt}
out=${in}.out
t=${2-6}
./router.debug "$in" "$out" "$t"
perl ./eval2008.pl "$in" "$out"
