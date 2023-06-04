#!/bin/bash
set -ex

cases=(
  "http://www.ispd.cc/contests/07/rcontest/benchmark/adaptec1.capo70.3d.35.50.90.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/adaptec2.mpl60.3d.35.20.100.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/adaptec3.dragon70.3d.30.50.90.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/adaptec4.aplace60.3d.30.50.90.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/adaptec1.capo70.3d.35.50.90.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/adaptec5.mfar50.3d.50.20.100.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/newblue1.ntup50.3d.30.50.90.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/newblue2.fastplace90.3d.50.20.100.gr.gz"
  "http://www.ispd.cc/contests/07/rcontest/benchmark/newblue3.kraftwerk80.3d.40.50.90.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/bigblue1.capo60.3d.50.10.100.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/bigblue2.mpl60.3d.40.60.60.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/bigblue3.aplace70.3d.50.10.90.m8.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/bigblue4.fastplace70.3d.80.20.80.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/newblue4.mpl50.3d.40.10.95.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/newblue5.ntup50.3d.40.10.100.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/newblue6.mfar80.3d.60.10.100.gr.gz"
  "http://www.ispd.cc/contests/08/benchmark/newblue7.kraftwerk70.3d.80.20.82.m8.gr.gz"
)
d="ISPD-2008-Benchmarks"

task() {
  url=${1?'url'}
  wget "$url"
  gunzip "$(basename "$url")"
}

mkdir "$d"
cd "$d"

for url in "${cases[@]}"; do
  task "$url" &
done
wait
