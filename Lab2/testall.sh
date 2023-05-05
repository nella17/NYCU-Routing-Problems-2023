#!/bin/bash
result=Benchmarks

route() {
  in=$1
  out=$1.out
  time=$1.txt
  /usr/bin/time -lp ./router "$in" "$out" 2>&1 | tee "$time"
}
test() {
  time perl ./eval2008.pl "$in" "$out" 2>&1 | tee -a "$result"
  rm "$out"
}

make -j || exit
rm "$result"
touch "$result"
for f in './ISPD 2008 Benchmarks'/*.gr; do
  echo "$f"
  route "$f"
  if [ $? ]; then
    if [ "$1" ]; then
      test "$f"
    else
      rm "$out"
    fi
  else
    echo "route fail"
  fi
done
wait
