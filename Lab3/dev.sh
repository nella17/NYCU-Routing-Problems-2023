#!/bin/bash
set -ex
DEBUG=1 make -j
i=${1-0}
c="./case/case$i.txt"
r="./case/result$i.txt"
time ./Lab3.debug "$c" "$r"
./verifier "$r" "$c"
python plotter.py "$r" $2
