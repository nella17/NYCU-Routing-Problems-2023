#!/bin/bash
set -ex
DEBUG=1 make -j
in=${1-./3d.txt}
out=${in}.out
./router.debug "$in" "$out"
perl ./eval2008.pl "$in" "$out"
