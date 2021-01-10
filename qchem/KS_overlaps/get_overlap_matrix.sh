#!/bin/bash

if [ $# -ne 2 ]; then
  echo "usage: 1.Q-Chem output file; 2.(matrix) output file to write"
  exit
fi

here=`pwd`

infile="$1"
outfile="$2"

#!# HARD-CODED grep strings:
topgrp='^ Overlap Matrix$'
botgrp='^ Core Hamiltonian$'


## get line number for each grep string:
top_li="$(grep -n "$topgrp" "$infile" | cut -d':' -f1)"
bot_li="$(grep -n "$botgrp" "$infile" | cut -d':' -f1)"

## get the number of AOs from the line above $bot_li:
numAOs="$(head -$((bot_li - 1)) "$infile" | tail -1 | awk '{print $1}')"
#    the rows are numbered in the 1st field of the overlap matrix ^^

## get the number of blocks:
numblocks="$(echo "($bot_li - $top_li - 1) / $numAOs" | bc)"


### now write blocks to separate temp files:

## array for temp block file names:
declare -a tmpblocks
## loop over blocks:
for bl in $(seq 1 $numblocks); do
  ## make temp file:
  tmp_bl="$(mktemp tmpblock.XXXXXXXX)"
  ## add file name to array:
  tmpblocks+=("$tmp_bl")

  ## number of lines to head:
  tohead="$(echo "$top_li + (($numAOs + 1) * $bl)" | bc)"

  #!# write block to $tmp_bl:
  head -$tohead "$infile" | tail -$numAOs | tr -s ' ' | cut -d' ' -f3- > "$tmp_bl"
done


## now combine all blocks into one file (ugh):
paste ${tmpblocks[@]} > "$outfile"


## clean up:
rm "${tmpblocks[@]}"


