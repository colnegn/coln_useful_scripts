#!/bin/bash

if [ $# -ne 3 ]; then
  echo "usage: 1.number of atoms (including frozens); 2.directory list of dirs containing force output files; 3.output file to write for all gradients (pasted together)"
  exit
fi

nats="$1"
dlist="$2"
full_outfile="$3"

if [ -f "$outfile" ]; then
  echo "ERROR: $outfile already exists... exiting..."
  exit
fi


## qchem prints out gradient components in blocks of 6 columns (one for each atom, with 3 rows for x,y,z):
n_blocks="$(echo "$nats / 6.0" | bc -l | cut -d'.' -f1)"
## get remainder block:
n_rem="$(echo "$nats % 6" | bc)"


## number of columns for each block:
n_cols_per_block=()
for bl in $(seq 1 $n_blocks); do
  n_cols_per_block+=(6)
done

## number of rows to grep -A for (4 rows per block):
li_per_block=4
n_grep="$((n_blocks * li_per_block))"
if [ $n_rem -gt 0 ]; then
  n_blocks=$((n_blocks + 1))
  n_grep="$((n_grep + li_per_block))"
  n_cols_per_block+=($n_rem)
fi


## keep list of single gradient files to paste together:
all_grad_files=''

grep_str='Gradient of SCF Energy'
for d in $(cat "$dlist"); do
  gradfile="$(mktemp ./tmp_grad.XXXXXXXX)"
  grep -A "$n_grep" "$grep_str" "$d/output" | tail -$n_grep > $gradfile


  #!# output file to write to:
  outfile="$d/gradient.gvec"
  #!# make sure not to append (want to overwrite):
  if [ -f "$outfile" ]; then
    echo "removing ${outfile} ..."
    rm "$outfile"
  fi

  ## add $outfile to $all_grad_files:
  all_grad_files="${all_grad_files} ${outfile}"


  ## loop over blocks in qchem output file:
  for bl in $(seq 1 $n_blocks); do
    blockfile="$(mktemp ./tmp_block.XXXXXXXX)"
    head -$((bl * li_per_block)) "$gradfile" | tail -$((li_per_block - 1)) > "$blockfile"
#                                                    ^^^ only want gradient components
    for col in $(seq 1 ${n_cols_per_block[$((bl - 1))]}); do
      cat "$blockfile" | awk -v n="$((col + 1))" '{print $n}' >> "$outfile"
    done

    ## clean up:
    rm "$blockfile"
  done

  ## clean up:
  rm "$gradfile"
done

      
## finally, paste all gradient files together
paste $all_grad_files > "$full_outfile"




