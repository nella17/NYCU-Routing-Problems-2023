#!/bin/bash
result=Benchmarks

function prepend() { while read -r line; do echo "${1}${line}"; done; }

route() {
  in=$1
  out=$1.out
  time=$1.txt
  /usr/bin/time -p ./router "$in" "$out" $((1)) 2>&1 | tee "$time" | prepend "$(basename "$in") "
  if [ $? ]; then
    echo "$1 done"
    time perl ./eval2008.pl "$in" "$out" 2>&1 | tee -a "$result"
    tail "$time" | tee -a "$result" > /dev/null
    rm "$out"
  else
    echo "route fail"
  fi
}

make -j || exit
rm "$result"
touch "$result"
for f in './ISPD 2008 Benchmarks'/*.gr; do
  echo "$f"
  route "$f" &
done
wait
