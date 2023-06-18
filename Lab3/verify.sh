#!/bin/bash
set -u

for i in {0..4}; do
  c="./case/case$i.txt"
  r="./case/result$i.txt"
  echo "case$i"
  ./verifier "$r" "$c"
  echo
done
