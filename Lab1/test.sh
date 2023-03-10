#!/bin/bash
set -eux
CXX=g++-12 make -j
fd . ./testbench -x bash -c '
./Lab1 {} ./output/{/} &&
python ./plotter.py \
  --in_file {} \
  --out_file ./output/{/} \
  --img_name ./output/{/.}.png
'
