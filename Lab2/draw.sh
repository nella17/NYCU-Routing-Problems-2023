#!/bin/bash
set -eux

mkdir draw || true
g++ -std=c++17 -Ofast -o draw.exe draw.cpp

function task() {
  f="$1"
  in=$f
  out=$f.all.out
  fn=$(basename "$f")
  cmap="draw/$fn"
  img="draw/$fn.png"
  time ./draw.exe "$in" "$out" "$cmap"
  python draw.py "$cmap" "$img"
}

for f in './ISPD-2008-Benchmarks'/*.gr; do
  echo "$f"
  task "$f" &
done
wait
