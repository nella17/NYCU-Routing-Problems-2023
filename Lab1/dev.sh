#!/bin/bash
set -eux
case=${1:-case1}
CXX=g++-12 make -j DEBUG=1
./Lab1 "./testbench/$case.txt" "./output/$case.txt"
python ./plotter.py \
  --in_file "./testbench/$case.txt" \
  --out_file "./output/$case.txt" \
  --img_name "./output/$case.png"
