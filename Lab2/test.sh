#!/bin/bash
set -ex
make -j
in=${1-'ISPD-2008-Benchmarks/newblue2.fastplace90.3d.50.20.100.gr'}
out=${in}.out
t=${2-600}
/usr/bin/time -p ./router "$in" "$out" "$t"
time perl ./eval2008.pl "$in" "$out"
rm "$out"
