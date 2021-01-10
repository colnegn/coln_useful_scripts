#!/bin/bash

if [ $# -ne 3 ]; then
  echo "usage: 1.Q-Chem output file; 2.(matrix) output file to write; 3.MO energy file to write"
  exit
fi

here=`pwd`

infile="$1"
outMOfile="$2"
energyfile="$3"


#!# HARD-CODED grep strings:
grpstr1='RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS'
grpstr2='eigenvalues: '


## line number for $grpstr1:
li1="$(grep -n "$grpstr1" "$infile" | cut -d':' -f1)"

## array for line numbers of tops of 'blocks' of matrix columns:
read -ra blocklines <<< "$(tail -n +$li1 "$infile" | grep -n "$grpstr2" | cut -d':' -f1 | tr '\n' ' ')"

## number of AOs from the number of rows of matrix:
numAOs="$(( ${blocklines[1]} - ${blocklines[0]} - 2 ))"
#                                                 ^ 2 lines between blocks


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
  head -$(( $li1 + $bl - 1)) "$infile" | tail -1 | tr -s ' ' | rev | cut -d' ' -f-$nrows_bl | rev | tr ' ' '\n' >> "$energyfile"
done


## now combine all blocks into one file (ugh):
paste ${tmpblocks[@]} > "$outMOfile"


## clean up:
rm "${tmpblocks[@]}"



