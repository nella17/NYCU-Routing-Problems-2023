#!/bin/bash
set -ex
CXX=g++-12 DEBUG=1 make -j
i=${1-0}
c="./case/case$i.txt"
r="./case/result$i.txt"
./Lab3.debug "$c" "$r"
python plotter.py "$r"
./verifier "$r" "$c"
