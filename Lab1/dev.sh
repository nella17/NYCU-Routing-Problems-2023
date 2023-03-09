#!/bin/bash
set -eux
CXX=g++-12 make -j DEBUG=1
./Lab1 ./testbench/case1.txt ./output/case1.txt
python ./plotter.py \
  --in_file ./testbench/case1.txt \
  --out_file ./output/case1.txt \
  --img_name ./output/case1.png
