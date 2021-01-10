#!/bin/bash


if [ $# -lt 1 ]; then
  printf '\n\n ~/alias_scripts/alias_tail.sh requires at least one argument\n\n\n'
  exit
fi


## if there is a line-number argument
if [ ${1:0:1} == '-' ]; then            # use 'getopts' builtin, instead?
  n_li_arg=$1

  ## shift input argument index:
  shift

## default: tail 10 lines
else
  n_li_arg='-10'
fi



## only want to print file names if we are tailing more than one:
if [ $# -gt 1 ]; then
  print_file_names=1
else
  print_file_names=0
fi


## loop over file arguments:
for i in `seq 1 $#`; do

  ## print file name:
  if [ $print_file_names -eq 1 ]; then
    echo ""
    echo "==> $1 <=="
  fi

  ## run tail:
  tail ${n_li_arg} $1
  ## prepare for next file:
  shift

done

exit



