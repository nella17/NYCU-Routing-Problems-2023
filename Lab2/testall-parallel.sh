#!/bin/bash
result=Benchmarks
tl=${1-$((30 * 60))}

function prepend() { while read -r line; do echo "${1}${line}"; done; }

route() {
  in=$1
  out=$1.all.out
  time=$1.txt
  timeout "$tl" /usr/bin/time -p ./router "$in" "$out" "$tl" 2>&1 | tee "$time" | prepend "$(basename "$in") "
  if [ $? ]; then
    echo "$1 done"
    case=$(basename "$in" | cut -d'.' -f1)
    rt=$(tail "$time" | grep real | awk '{ print $2 }')
    time perl ./eval2008.pl "$in" "$out" 2>&1 | tee -a "$time" | \
      tail -n1 | tail "-c+40" | \
      awk "{ print \"$case\" \"\t\" \$1 \"\t\" \$2 \"\t\" \$3 \"\t\" \"$rt\" }" | tee -a "$result"
    rm "$out"
  else
    echo "route fail"
  fi
}

make -j || exit
echo -e "case\tTOF\tMOF\tWL\tRT" > "$result"
for f in './ISPD-2008-Benchmarks'/*.gr; do
  echo "$f"
  route "$f" &
done
wait
