#!/bin/bash

outstr=''
for i in $(seq 1 $#); do
  outstr="${outstr}$1"
  shift
done

printf "$outstr"

