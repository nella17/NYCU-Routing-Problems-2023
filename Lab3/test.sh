#!/bin/bash
set -ux
make -j

for i in {0..5}; do
  c="./case/case$i.txt"
  r="./case/result$i.txt"
  p="./case/result$i.png"
  time ./Lab3 "$c" "$r" 2>/dev/null
  if [ $? = 0 ]; then
    python plotter.py "$r"
    mv path.png "$p"
    ./verifier "$r" "$c"
  else
    rm "$r" "$p"
  fi
done
