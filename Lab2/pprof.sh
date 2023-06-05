#!/bin/bash
set -ex

task='newblue5*'

export CPUPROFILE=./prof.out
export CPUPROFILE_FREQUENCY=10000

# make clean
if [ "$(uname -s)" != "Darwin" ]; then
  PERF=1 CXX=g++ PF_FLAGS='-Ofast -g -lprofiler' \
    make -j
else
	PERF=1 CXX=g++-12 PF_FLAGS='-Ofast -g -L/opt/homebrew/Cellar/gperftools/2.10/lib -lprofiler' \
		make -j
fi

time ./router.perf ./ISPD-2008-Benchmarks/$task.gr /dev/null
# pprof ./router.perf ./prof.out
pprof --pdf ./router.perf ./prof.out > out.pdf
