#!/bin/bash
set -eux
make -j

for i in {1..5}; do
  c="./case/case$i.txt"
  r="./case/result$i.txt"
  ./Lab3 "$c" "$r"
  python plotter.py "$r"
  ./verifier "$r" "$c"
done
