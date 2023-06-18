#!/bin/bash
set -ux
make -j

task() {
  c="./case/case$1.txt"
  r="./case/result$1.txt"
  p="./case/result$1.png"
  s="./case/result$1.sat"
  st="./case/result$1.sat.txt"

  d=$(mktemp -d)
  cd "$d" || exit
  ln -s "$OLDPWD/"{Lab3,'case',open-wbo,plotter.py,verifier} .

  if time ./Lab3 "$c" "$r" 2>/dev/null; then
    cp clause.sat "$s"
    cp sat_result.txt "$st"
    python plotter.py "$r"
    mv path.png "$p"
    ./verifier "$r" "$c"
  else
    rm "$r" "$p"
  fi

  cd "$OLDPWD" || exit
  rm -fr "$d"
}

for i in {0..5}; do
  l="./case/test$i.log"
  task "$i" >"$l" 2>&1 &
done
wait
