#!/bin/bash

ion='fluoride'
spec='Isomer1/sp_012'
num_occ_sporbs='15'
outdat='f_012_wf_overlaps.dat'

ion_list=('fluoride' 'chloride' 'bromide')
nocc_list=('15' '19' '28')

if [ $# -ne 1 ]; then
  echo "usage: 1.list of density functionals (ordered)"
  exit
fi

flist="$1"
## number of functionals:
numfn="$(wc -l "$flist" | awk '{print $1}')"

## data file header:
header='#! ROWS/COLUMNS:\n#!'
header="${header} $(grep -n '^' "$flist" | tr '\n' ' ')"

for i in $(seq 0 2); do
  ion="${ion_list[$i]}"
  num_occ_sporbs="${nocc_list[$i]}"
#for ion in fluoride chloride bromide; do #  iodide

  outdat="${ion}_012_wf_overlaps.dat"

  
  ## bash array to hold wf overlap values:
  wf_overlap_vals=()
  
  ## loop over pairs of functionals:
  for i in $(seq 1 "$numfn"); do 
    for j in $(seq "$i" "$numfn"); do 
      
      ## get directory names for functionals $i, and $j:
      Fni="$(head -$i "$flist" | tail -1)"
      Fnj="$(head -$j "$flist" | tail -1)"
  
      ## get paths to AO-to-MO coef matrix for $ion/$spec, for each functional:
      Pathi="$ion/${Fni}_wat2/$spec/${ion}_${Fni}_wat2_$(dirname "$spec")_$(basename "$spec").Cmat"
      Pathj="$ion/${Fnj}_wat2/$spec/${ion}_${Fnj}_wat2_$(dirname "$spec")_$(basename "$spec").Cmat"
      
      ## get path to AO overlap matrix for $ion/$spec:
      Patho="$ion/overlap/$spec/${ion}_overlap_$(dirname "$spec")_$(basename "$spec").Smat"
  
      
      ## add value of overlap ij to array:
      wf_overlap_vals+=("$(printf "%0.4f" $(./calc_wf_overlap.py "$Pathi" "$Pathj" "$Patho" "$num_occ_sporbs"))")
      
    done
  done


  ## print header:
  echo -e "$header" >> "$outdat"
  
  ## now write matrix:
  for i in $(seq 0 "$((numfn-1))"); do 
    for j in $(seq 0 "$((numfn-1))"); do 
      ## calculate element of $wf_overlap_vals for element ij:
      el="$(echo "n=float($i);m=float($j);k=float($numfn);print(int(min(n,m)*(k+((1-min(n,m))/2))+abs(n-m)))" | python)"
  
      printf "${wf_overlap_vals[$el]}  " >> "$outdat"
    done
    printf '\n' >> "$outdat"
  done


done


