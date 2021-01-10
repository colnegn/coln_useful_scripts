#!/bin/bash


if [ $# -lt 1 ]; then
  printf '\n\n ~/alias_scripts/alias_open.sh requires at least one argument\n\n\n'
  exit
fi



for i in `seq 1 $#`; do

  xdg-open $1 &> /dev/null

  shift

done

exit



