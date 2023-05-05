#!/bin/bash
set -ex
make -j
in=${1-'ISPD 2008 Benchmarks/adaptec1.gr'}
out=${in}.out
/usr/bin/time -lp ./router "$in" "$out"
time perl ./eval2008.pl "$in" "$out"
rm "$out"
