#!/bin/bash
set -u

for i in {0..4}; do
  c="./case/case$i.txt"
  r="./case/result$i.txt"
  l="./case/test$i.log"
  echo "case$i"
  grep real "$l"
  ./verifier "$r" "$c"
  echo
done
