#!/bin/bash

if [ $# -ne 3 ]; then
  echo "usage: 1.Q-Chem output file; 2.prefix for (matrix) output file to write; 3.prefix for MO energy file to write"
  exit
fi

here=`pwd`

infile="$1"
outMOfile="$2"
energyfile="$3"


#!# HARD-CODED grep strings:
grpstr1='RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS'
grpstr2='eigenvalues: '


## line numbers for $grpstr1:
read -ra li1s <<< "$(grep -n "$grpstr1" "$infile" | cut -d':' -f1 | tr '\n' ' ')"
## number of systems:
nsys="${#li1s[@]}"

#!# add number of lines in file to array for tail (4 commands down):
li1s+=("$(wc -l "$infile" | awk '{print $1}')")

#!# LOOP OVER DIFFERENT systems (their AO-to-MO matrices):
for syi in $(seq 0 $((nsys - 1))); do 

  li1=${li1s[$syi]}
  tohead="$(echo "${li1s[$((syi + 1))]} - $li1" | bc)"

  ##!#!#!#!#
  ##!# TEST:
  #echo "${li1s[@]}"
  #echo "$syi"
  #echo "$li1"
  #echo "$tohead"
  #echo
  ##!#!#!#!#


  ## array for line numbers of tops of 'blocks' of matrix columns:
  read -ra blocklines <<< "$(tail -n +$li1 "$infile" | head -$tohead | grep -n "$grpstr2" | cut -d':' -f1 | tr '\n' ' ')"
  #read -ra blocklines <<< "$(tail -n +$li1 "$infile" | grep -n "$grpstr2" | cut -d':' -f1 | tr '\n' ' ')"

  ## number of AOs from the number of rows of matrix:
  #!# if only one block:
  if [ "${#blocklines[@]}" -eq 1 ]; then
    li_2_below_AOS="$(tail -n +$li1 "$infile" | head -$tohead | grep -n 'Ground-State Mulliken Net Atomic Charges' | cut -d':' -f1 | head -1)"
    numAOs="$(( $li_2_below_AOS - ${blocklines[0]} - 2 ))"
  else
    numAOs="$(( ${blocklines[1]} - ${blocklines[0]} - 2 ))"
    #                                                 ^ 2 lines between blocks
  fi

  ##!#!#!#!#
  ##!# TEST:
  #echo "${blocklines[@]}"
  #echo "$numAOs"
  #echo
  ##!#!#!#!#
  
  
  ### now write blocks to separate temp files:
  
  ## array for temp block file names:
  declare -a tmpblocks
  ## loop over blocks:
  for bl in ${blocklines[@]}; do
    ## make temp file:
    tmp_bl="$(mktemp tmpblock.XXXXXXXX)"
    ## add file name to array:
    tmpblocks+=("$tmp_bl")
  
    #!# count number of rows in block, $bl (at most 6):
    nrows_bl="$(( $(tail -n +$(( $li1 + $bl -1 )) "$infile" | head -1 | wc -w) - 1 ))"
  
    #!# write block to $tmp_bl:
    tail -n +$(( $li1 + $bl )) "$infile" | head -$numAOs | tr -s ' ' | rev | cut -d' ' -f-$nrows_bl | rev > "$tmp_bl"
  
    #!# write MO energies for block to $energyfile:
    head -$(( $li1 + $bl - 1)) "$infile" | tail -1 | tr -s ' ' | rev | cut -d' ' -f-$nrows_bl | rev | tr ' ' '\n' >> "${energyfile}_${syi}"
  done
  
  
  ## now combine all blocks into one file (ugh):
  paste ${tmpblocks[@]} > "${outMOfile}_${syi}"
  
  
  ## clean up:
  rm "${tmpblocks[@]}"
  unset blocklines
  unset tmpblocks

done


