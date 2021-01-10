#!/bin/bash

usage_str='
usage:
       1.list of directories (dlist) containing file to grep inside;
       2.name of file to grep (name must be the same for all directories in the dlist);
       3.string to grep for (grep string) in quotes (must be the same for all files;
       4.(optional) number of occurrences of grep string expected in each file
'

## expect 1 occurrence of grep string by default:
numgrep=1

if [ $# -ne 3 ]; then
  if [ $# -ne 4 ]; then
    echo "$usage_str"
    exit
  fi
  ## must be 4 arguments:
  numgrep="$4"
fi

dlist="$1"
fname="$2"
grepstr="$3"


echo "calculations to resubmit:"

## loop over directories in $dlist:
for d in $(cat "$dlist"); do
  ## make sure that $fname exists in $d:
  if ! [ -f "$d/$fname" ]; then
    echo "$d"
    continue
  fi

  ## count occurrences of $grepstr:
  if [ "$(grep -c "$grepstr" "$d/$fname")" -ne "$numgrep" ]; then
    ## print failing directories to stdout:
    echo "$d"
  fi
done

echo "^^^^^^^^^"


