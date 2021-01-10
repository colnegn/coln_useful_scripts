#!/bin/bash


if [ $# -lt 1 ]; then
  echo 'usage: 1.xyz file (containing configurations with same number of atoms!!!)'
  exit
fi

sum_count=0

for i in `seq 1 $#`; do

  ## xyz file, $i:
  inxyz=$1
  
  ## number of atoms in $inxyz:
  nats=$(head -1 $inxyz | sed -e 's/\r//')

  
  ## make sure same number of atoms in each configuration (and catch simple defects):
  if [ $(awk -v n=$((nats + 2)) '(NR - 1) % n == 0' $inxyz | sed 's/^ * //g;s/ * $//g' | sort -g | uniq | wc -l) -ne 1 ]; then
    echo "ERROR: $inxyz has multiple systems, or a defect... exiting..."
    exit
  fi
  
  
  ## if we made it here, xyz file seems fine, so count the configurations:
  echo -e "\n===>  $inxyz  <==="
  ncfs="$(echo "$(wc -l $inxyz | awk '{print $1}') / $((nats + 2))" | bc)"
  echo -e "${ncfs}\n"


  ## sum number of configs for all xyzs:
  sum_count=$((sum_count + ncfs))

  ## prepare for next xyz file:
  shift
done

echo -e "\ntotal sum: $sum_count"


