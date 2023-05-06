#!/bin/bash
set -ex
make -j
in=${1-'ISPD 2008 Benchmarks/newblue2.gr'}
out=${in}.out
/usr/bin/time -p ./router "$in" "$out"
time perl ./eval2008.pl "$in" "$out"
rm "$out"
