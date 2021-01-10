#!/bin/bash

if [ $# -ne 5 ]; then
  echo "usage: 1.input xyz file; 2.output sorted xyz file (to overwrite); 3.field number of comment lines to sort by; 4.histogram file to write; 5.number to subtract from each value (use 0.0 if none)"
  exit
fi

inxyz="$1"
outxyz="$2"
fnum="$3"
outhist="$4"
subval="$5"

if [ -f "$outxyz" ]; then
  echo "ERROR: $outxyz already exists... exiting..."
  exit
fi
if [ -f "$outhist" ]; then
  echo "ERROR: $outhist already exists... exiting..."
  exit
fi


## get number of atoms from first line of xyz:
nats="$(head -1 "$inxyz")"
natsp2=$((nats + 2))

tmpfile="$(mktemp tmpvalues.XXXXXXXX)"
## get values:
awk "(NR - 1) % ${natsp2} == 1" "$inxyz" | awk -v n="$fnum" '{print $n}'  > "$tmpfile"

## get sorted order:
pyscript="
import numpy as np
invals = np.genfromtxt('${tmpfile}')
invals -= $subval
with open('order_${tmpfile}', 'w') as f:
  for i in np.argsort(invals):
    f.write(str(i) + ' ' + str('{:0.8f}'.format(invals[i])) + '\\n')
"
## run sort script:
echo "$pyscript" | python3

## clean up:
rm "$tmpfile"

keep_comment_line_bool=1

## now sort the xyz file:
indat="order_${tmpfile}"
for i in $(seq 1 "$(wc -l "$indat" | awk '{print $1}')"); do
  cfi="$(head -$i "$indat" | tail -1 | awk '{print $1}')"
  ## change to 1-index:
  cfi=$((cfi + 1))
  vali="$(head -$i "$indat" | tail -1 | awk '{print $2}')"
  if [ $keep_comment_line_bool == 1 ]; then
    head -$((cfi * natsp2)) "$inxyz" | tail -$natsp2 >> "$outxyz"
  else
    echo "$nats" >> "$outxyz"
    echo "$vali" >> "$outxyz"
    head -$((cfi * natsp2)) "$inxyz" | tail -$nats >> "$outxyz"
  fi
done

## clean up:
rm "$indat"

## finally, write a histogram file for good measure:
awk "(NR - 1) % ${natsp2} == 1" "$outxyz" | grep -n '^' | sed 's/:/ /g' > "$outhist"



