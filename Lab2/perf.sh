#!/bin/bash
set -ex

task='newblue5*'

FREQ=1000

# make clean
PERF=1 PF_FLAGS='-Ofast -g' \
  make -j

rm perf.data || true
time perf record -F $FREQ -g --call-graph dwarf -- \
  ./router.perf ./ISPD-2008-Benchmarks/$task.gr /dev/null
perf script > out.perf

# https://github.com/brendangregg/FlameGraph.git
./FlameGraph/stackcollapse-perf.pl out.perf > out.folded
./FlameGraph/flamegraph.pl out.folded > out.svg
