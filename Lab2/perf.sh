#!/bin/bash
set -ex
PERF=1 CXX=g++-12 PF_FLAGS='-Ofast -g -L/opt/homebrew/Cellar/gperftools/2.10/lib -lprofiler' \
  make -j
CPUPROFILE=./prof.out time ./router.perf 'ISPD 2008 Benchmarks/adaptec3.gr' /dev/null
# pprof ./router.perf ./prof.out
pprof --pdf ./router.perf ./prof.out > out.pdf
