#!/bin/bash


## warning:
(>&2 echo "WARNING: indexing starts with 1")

## default delimiter ' ':
delim=' '

if [ $# -ne 3 ]; then
  if [ $# -ne 4 ]; then
    echo "usage: 1.number of elements (to take combinations of); 2.combo-length lower-bound; 3.combo-length upper-bound; 4.(optional) delimited character (default is a single space)"
    exit
  fi

  ## must be 3 arguments:
  delim="$4"
fi

## bounds to lengths of combinations:
Ne="$1"
N_lw="$2"
N_up="$3"


## string to pipe into python for combo gen:
cg_script="
from itertools import combinations
for N in range(${N_lw},$((N_up+1))):
  for i in combinations(range(1, $((Ne+1))), N):
    print('${delim}'.join(list(map(str, i))))
"


## pipe into python, printing to stdout:
echo "$cg_script" | python




